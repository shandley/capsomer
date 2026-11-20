-- Setup: read clusters that no search assigned to a taxon.
CREATE TEMP VIEW dark AS
SELECT s.sequence_id FROM sequence s
WHERE s.kind = 'read_cluster' AND NOT EXISTS (
  SELECT 1 FROM taxon_assignment t WHERE t.sequence_id = s.sequence_id AND t.outcome = 'assigned');

-- Setup: read clusters and the contig each maps to.
CREATE TEMP VIEW on_contig AS
SELECT query_sequence_id AS sequence_id, target_sequence_id AS contig_id, target_start
FROM alignment WHERE target_sequence_id IS NOT NULL AND is_primary;

-- What the bundle holds.
SELECT kind, count(*) AS sequences FROM sequence GROUP BY 1 ORDER BY 1;

-- Question 8. What fraction of each library's fragments matched nothing in any search,
-- and what fraction is unmatched but mapped to a contig?
SELECT c.library_id, sum(c.count) AS fragments,
  round(sum(c.count) FILTER (WHERE d.sequence_id IS NOT NULL) / sum(c.count), 3) AS unmatched,
  round(sum(c.count) FILTER (WHERE d.sequence_id IS NOT NULL AND o.sequence_id IS NOT NULL) / sum(c.count), 3) AS unmatched_on_a_contig
FROM read_count c LEFT JOIN dark d USING (sequence_id) LEFT JOIN on_contig o USING (sequence_id)
GROUP BY 1 ORDER BY 1;

-- Question 9. Which unmatched sequences recur in more than one sample?
-- By identical sequence, and by mapping to a contig that unmatched reads from several samples share.
SELECT
  (SELECT count(*) FROM (SELECT sequence_id FROM read_count JOIN dark USING (sequence_id)
     GROUP BY 1 HAVING count(DISTINCT library_id) > 1)) AS identical_in_several_samples,
  (SELECT count(*) FROM (SELECT o.contig_id FROM on_contig o JOIN dark USING (sequence_id)
     JOIN read_count c USING (sequence_id) GROUP BY 1 HAVING count(DISTINCT c.library_id) > 1)) AS contigs_shared_by_unmatched_reads;

-- Question 3. Which read clusters have only translated hits below 40% identity and map to no contig?
SELECT count(*) AS read_clusters, sum(c.count) AS fragments, round(min(w.best_identity), 3) AS lowest_identity
FROM (SELECT query_sequence_id AS sequence_id, max(identity) AS best_identity FROM alignment
      WHERE target_reference_entry_id IS NOT NULL GROUP BY 1
      HAVING bool_and(search_type = 'translated_nucleotide_to_protein') AND max(identity) < 0.40) w
JOIN (SELECT sequence_id, sum(count) AS count FROM read_count GROUP BY 1) c USING (sequence_id)
WHERE w.sequence_id NOT IN (SELECT sequence_id FROM on_contig);

-- Weak hits are kept: translated hits by identity band, with the poorest e-value in each.
SELECT CASE WHEN identity < 0.3 THEN 'under 30%' WHEN identity < 0.5 THEN '30 to 50%'
            WHEN identity < 0.7 THEN '50 to 70%' ELSE '70% and over' END AS amino_acid_identity,
  count(*) AS hits, max(evalue) AS poorest_evalue
FROM alignment WHERE search_type = 'translated_nucleotide_to_protein' GROUP BY 1 ORDER BY 1;

-- Question 5. Viral fragments per library by family, summed from read clusters through their final taxon.
SELECT t.family_name, sum(c.count) AS fragments, count(DISTINCT c.library_id) AS libraries
FROM read_count c JOIN taxon_assignment a ON a.sequence_id = c.sequence_id AND a.is_preferred
JOIN taxon t USING (taxon_id) WHERE t.domain_name = 'Viruses'
GROUP BY 1 ORDER BY 2 DESC LIMIT 10;

-- Question 13. Which unmatched read clusters map to a contig that has a taxon?
SELECT t.domain_name, t.family_name, count(*) AS unmatched_read_clusters
FROM dark d JOIN on_contig o USING (sequence_id)
JOIN taxon_assignment a ON a.sequence_id = o.contig_id AND a.is_preferred
JOIN taxon t USING (taxon_id) GROUP BY 1, 2 ORDER BY 3 DESC LIMIT 10;

-- Question 12. Taxa along a contig: contigs whose two halves carry different families of reads.
-- A difference may mean recombination or a chimeric assembly.
WITH hits AS (
  SELECT o.contig_id, o.target_start > s.length / 2 AS right_half, t.family_name
  FROM on_contig o JOIN sequence s ON s.sequence_id = o.contig_id
  JOIN taxon_assignment a ON a.sequence_id = o.sequence_id AND a.is_preferred
  JOIN taxon t USING (taxon_id) WHERE t.family_name IS NOT NULL),
halves AS (SELECT contig_id, right_half, family_name, count(*) AS n FROM hits GROUP BY ALL),
top AS (SELECT contig_id, right_half, arg_max(family_name, n) AS family, sum(n) AS reads,
          max(n) / sum(n) AS share FROM halves GROUP BY 1, 2)
SELECT s.name AS contig, s.length, l.family AS first_half_family, l.reads AS first_half_reads,
  r.family AS second_half_family, r.reads AS second_half_reads
FROM top l JOIN top r ON l.contig_id = r.contig_id AND NOT l.right_half AND r.right_half
JOIN sequence s ON s.sequence_id = l.contig_id
WHERE l.family <> r.family AND l.reads >= 5 AND r.reads >= 5 AND l.share >= 0.8 AND r.share >= 0.8
ORDER BY least(l.reads, r.reads) DESC LIMIT 10;

-- Question 6. What does the bundle record about how it was made?
SELECT r.run_id, r.step, p.name AS setting, p.value
FROM run r LEFT JOIN run_parameter p USING (run_id) ORDER BY 1, 3;

-- Question 1 needs a host category for each taxon. How many taxa have one?
SELECT count(*) AS taxa, count(host_category) AS with_host_category, count(native_id) AS with_release_identifier
FROM taxon;
