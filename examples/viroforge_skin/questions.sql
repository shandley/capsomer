-- Setup: the true source of each contig, from its primary alignment to the simulated genomes.
CREATE VIEW truth AS
SELECT a.query_sequence_id AS sequence_id, a.target_sequence_id AS true_sequence_id,
  a.query_coverage, a.target_coverage, a.identity,
  (SELECT min(origin) FROM x_sequence_origin o WHERE o.sequence_id = a.target_sequence_id) AS origin,
  (SELECT bool_and(is_known) FROM x_sequence_origin o WHERE o.sequence_id = a.target_sequence_id) AS has_family
FROM alignment a WHERE a.run_id = 'minimap2.truth' AND a.is_primary
QUALIFY row_number() OVER (PARTITION BY a.query_sequence_id ORDER BY a.alignment_length DESC) = 1;

-- Do the viral calls agree with the truth? Contigs by assembler, geNomad call and true origin.
SELECT s.run_id AS assembly, v.call, coalesce(t.origin, 'no alignment to a simulated genome') AS true_origin,
  count(*) AS contigs
FROM sequence_source s JOIN viral_call v ON v.sequence_id = s.sequence_id AND v.run_id = 'genomad.' || s.run_id
LEFT JOIN truth t ON t.sequence_id = s.sequence_id
WHERE s.run_id IN ('megahit', 'metaspades') GROUP BY ALL ORDER BY 1, 2, 4 DESC;

-- Does estimated completeness track true completeness? Viral contigs by CheckV tier.
SELECT q.run_id, q.tier, count(*) AS contigs, round(avg(q.completeness), 3) AS mean_checkv_completeness,
  round(avg(t.target_coverage), 3) AS mean_true_completeness
FROM quality q JOIN truth t USING (sequence_id)
WHERE q.run_id LIKE 'checkv.%' GROUP BY ALL ORDER BY 1, 4 DESC;

-- Viral genomes with no family in the simulator's taxonomy: were they assembled, called viral, and named?
SELECT t.has_family, count(DISTINCT t.true_sequence_id) AS genomes_with_a_contig,
  count(*) AS contigs, count(*) FILTER (WHERE v.call = 'virus') AS called_virus,
  count(*) FILTER (WHERE a.outcome = 'assigned') AS given_a_taxon
FROM truth t JOIN viral_call v USING (sequence_id) JOIN taxon_assignment a USING (sequence_id, run_id)
WHERE t.origin = 'viral' GROUP BY 1 ORDER BY 1;

-- The same contig from two assemblers is one sequence with two sources.
SELECT s.sequence_id, q.length, string_agg(s.run_id || ': ' || s.name, '; ') AS produced_by
FROM sequence_source s JOIN sequence q USING (sequence_id)
WHERE s.run_id IN ('megahit', 'metaspades') GROUP BY ALL HAVING count(*) > 1;

-- Species clusters: how many hold contigs from both assemblers?
SELECT n_assemblies, count(*) AS clusters, sum(members) AS member_contigs FROM (
  SELECT c.cluster_id, count(DISTINCT s.run_id) AS n_assemblies, count(DISTINCT m.sequence_id) AS members
  FROM cluster c JOIN cluster_member m USING (cluster_id)
  JOIN sequence_source s ON s.sequence_id = m.sequence_id AND s.run_id IN ('megahit', 'metaspades')
  WHERE c.run_id = 'vclust' GROUP BY 1) GROUP BY 1 ORDER BY 1;

-- Connectivity: contigs whose mapped reads were simulated from more than one genome.
-- Each read's true genome comes from the simulator; its position on the contig from the mapper.
WITH placed AS (
  SELECT m.target_sequence_id AS contig_id, tr.target_sequence_id AS true_sequence_id, count(*) AS reads,
    min(m.target_start) AS first_base, max(m.target_end) AS last_base
  FROM alignment m JOIN alignment tr ON tr.query_read_id = m.query_read_id AND tr.run_id = 'viroforge'
  WHERE m.run_id = 'minimap2.reads' AND m.is_primary GROUP BY 1, 2)
SELECT (SELECT min(name) FROM sequence_source x WHERE x.sequence_id = p.contig_id) AS contig, q.length,
  count(*) AS true_genomes, sum(reads) AS reads,
  string_agg(reads || ' reads at ' || first_base || '-' || last_base, '; ' ORDER BY reads DESC) AS by_true_genome
FROM placed p JOIN sequence q ON q.sequence_id = p.contig_id
WHERE reads >= 20 GROUP BY p.contig_id, q.length HAVING count(*) > 1 ORDER BY 4 DESC LIMIT 8;

-- Abundance: reads counted on each species cluster against the fragments simulated from its true genome.
SELECT (SELECT min(name) FROM sequence_source x WHERE x.sequence_id = c.sequence_id AND x.run_id <> 'viroforge') AS representative,
  c.mapped_count AS reads_mapped, 2 * tc.mapped_count AS reads_simulated,
  round(c.mean_depth, 1) AS mean_depth, round(t.target_coverage, 2) AS share_of_true_genome
FROM coverage c JOIN truth t USING (sequence_id)
JOIN coverage tc ON tc.sequence_id = t.true_sequence_id AND tc.run_id = 'viroforge'
WHERE c.run_id = 'coverm' ORDER BY 2 DESC LIMIT 8;
