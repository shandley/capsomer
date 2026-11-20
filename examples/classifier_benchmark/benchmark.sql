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

-- Uncultivated genomes by study, since the study's own discovery tool is favoured.
SELECT tool, left(study, 45) AS study, count(*) AS pieces, round(avg(called_viral::INT), 3) AS found
FROM scored WHERE test_set = 'recent_uncultivated' AND study IN (
  SELECT study FROM x_truth WHERE test_set = 'recent_uncultivated' GROUP BY 1 HAVING count(*) >= 60)
GROUP BY ALL ORDER BY 2, 4 DESC;

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
