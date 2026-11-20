-- Setup: one row per piece per classifier, with what is known about the piece.
CREATE VIEW scored AS
SELECT t.*, split_part(v.run_id, '.', 1) AS tool, v.call = 'virus' AS called_viral,
  CASE WHEN t.refseq_close_coverage >= 0.5 THEN '1 close to RefSeq'
       WHEN t.refseq_protein_coverage >= 0.5 THEN '2 distant'
       WHEN t.refseq_protein_coverage IS NOT NULL THEN '3 very distant' END AS novelty
FROM x_truth t JOIN viral_call v USING (sequence_id) WHERE v.run_id <> 'pieces';

-- Genomes deposited after every tool's database date: share of pieces called viral, by novelty.
SELECT tool, novelty, count(*) AS pieces, round(avg(called_viral::INT), 3) AS found
FROM scored WHERE test_set LIKE 'recent%' GROUP BY ALL ORDER BY 2, 4 DESC;

-- The same, split into cultured phages and uncultivated genomes.
SELECT tool, test_set, novelty, count(*) AS pieces, round(avg(called_viral::INT), 3) AS found
FROM scored WHERE test_set LIKE 'recent%' GROUP BY ALL ORDER BY 2, 3, 5 DESC;

-- How many genomes and studies stand behind each set of pieces?
SELECT test_set, count(*) AS pieces, count(DISTINCT source) AS genomes, count(DISTINCT study) AS studies
FROM x_truth GROUP BY 1 ORDER BY 1;

-- Very distant pieces of the recent set, by deposit. One deposit supplies most of them.
SELECT CASE WHEN study LIKE 'KCR_vOTU%' THEN 'KCR (cucumber rhizosphere)' WHEN study LIKE 'WSV%' THEN 'WSV'
  WHEN test_set = 'recent_cultured' THEN 'cultured phages' ELSE 'other uncultivated' END AS deposit,
  tool, count(*) AS pieces, round(avg(called_viral::INT), 3) AS found
FROM scored WHERE test_set LIKE 'recent%' AND novelty = '3 very distant' GROUP BY ALL ORDER BY 1, 4 DESC;

-- A second set with many studies: uncultivated genomes deposited April to September 2026,
-- at most 25 from one study. Share of pieces called viral, by novelty.
SELECT tool, novelty, count(*) AS pieces, round(avg(called_viral::INT), 3) AS found
FROM scored WHERE test_set = 'earlier_uncultivated' GROUP BY ALL ORDER BY 2, 4 DESC;

-- In that set, studies with at least five very distant pieces: which of geNomad and Jaeger found more?
WITH wide AS (SELECT sequence_id, any_value(study) AS study,
    bool_or(called_viral) FILTER (WHERE tool = 'genomad') AS genomad,
    bool_or(called_viral) FILTER (WHERE tool = 'jaeger') AS jaeger
  FROM scored WHERE test_set = 'earlier_uncultivated' AND novelty = '3 very distant' GROUP BY 1),
per_study AS (SELECT study, avg(genomad::INT) AS genomad, avg(jaeger::INT) AS jaeger, count(*) AS pieces
  FROM wide GROUP BY 1 HAVING count(*) >= 5)
SELECT count(*) AS studies, sum(pieces) AS pieces, count(*) FILTER (WHERE jaeger > genomad) AS jaeger_ahead,
  count(*) FILTER (WHERE genomad > jaeger) AS genomad_ahead, count(*) FILTER (WHERE genomad = jaeger) AS tied,
  round(avg(genomad), 3) AS mean_genomad, round(avg(jaeger), 3) AS mean_jaeger
FROM per_study;

-- Known RefSeq viruses from the simulation.
SELECT tool, count(*) AS pieces, round(avg(called_viral::INT), 3) AS found
FROM scored WHERE test_set = 'simulated_refseq' GROUP BY ALL ORDER BY 3 DESC;

-- Non-viral pieces called viral. A piece several tools call is probably a prophage,
-- so calls a tool makes alone are the better estimate of its errors.
WITH agree AS (SELECT sequence_id, count(*) FILTER (WHERE called_viral) AS n FROM scored WHERE NOT is_viral GROUP BY 1)
SELECT s.tool, count(*) AS non_viral_pieces, count(*) FILTER (WHERE s.called_viral) AS called_viral,
  count(*) FILTER (WHERE s.called_viral AND a.n = 1) AS called_by_this_tool_alone,
  count(*) FILTER (WHERE s.called_viral AND a.n >= 4) AS shared_with_three_or_more_tools
FROM scored s JOIN agree a USING (sequence_id) WHERE NOT s.is_viral GROUP BY 1 ORDER BY 4;

-- Does running a second tool beside geNomad help on the very distant pieces?
WITH wide AS (SELECT sequence_id, bool_or(called_viral) FILTER (WHERE tool = 'genomad') AS genomad,
  bool_or(called_viral) FILTER (WHERE tool = 'jaeger') AS jaeger,
  bool_or(called_viral) FILTER (WHERE tool = 'phabox2') AS phabox2,
  bool_or(called_viral) FILTER (WHERE tool = 'virsorter2') AS virsorter2,
  bool_or(called_viral) FILTER (WHERE tool = 'cenotetaker3') AS cenotetaker3
  FROM scored WHERE novelty = '3 very distant' GROUP BY 1)
SELECT count(*) AS very_distant_pieces, round(avg(genomad::INT), 3) AS genomad_alone,
  round(avg((genomad OR jaeger)::INT), 3) AS with_jaeger,
  round(avg((genomad OR phabox2)::INT), 3) AS with_phabox2,
  round(avg((genomad OR virsorter2)::INT), 3) AS with_virsorter2,
  round(avg((genomad OR cenotetaker3)::INT), 3) AS with_cenotetaker3,
  round(avg((genomad OR jaeger OR phabox2 OR virsorter2 OR cenotetaker3)::INT), 3) AS any_tool
FROM wide;
