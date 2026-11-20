"""Write the scored benchmark as plain tables, so results can be re-analysed without rerunning any tool.

Usage: export_results.py <bundle folder> <results folder>
Writes calls.tsv: one row per test piece, with what is known about it and each classifier's
call and viral score.
"""

import sys
from pathlib import Path

from capsomer.bundle import connect, literal

TOOLS = ("genomad", "jaeger", "virsorter2", "phabox2", "cenotetaker3")


def main() -> None:
    bundle, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    con = connect(bundle)
    per_tool = ", ".join(
        f"max(v.call) FILTER (WHERE split_part(v.run_id, '.', 1) = '{t}') AS {t}_call, "
        f"max(v.virus_score) FILTER (WHERE split_part(v.run_id, '.', 1) = '{t}') AS {t}_virus_score"
        for t in TOOLS
    )
    con.execute(
        "COPY (SELECT s.name AS piece, t.test_set, t.is_viral, t.origin, t.source, t.study, t.has_family, "
        "round(t.refseq_protein_coverage, 4) AS refseq_protein_coverage, "
        f"round(t.refseq_close_coverage, 4) AS refseq_close_coverage, {per_tool} "
        "FROM x_truth t JOIN sequence_source s ON s.sequence_id = t.sequence_id AND s.run_id = 'pieces' "
        "LEFT JOIN viral_call v ON v.sequence_id = t.sequence_id "
        "GROUP BY ALL ORDER BY t.test_set, s.name) "
        f"TO {literal(str(out / 'calls.tsv'))} (HEADER, DELIMITER '\t')"
    )
    rows = con.execute("SELECT count(*) FROM x_truth").fetchone()
    print(f"wrote {out / 'calls.tsv'}: {rows[0] if rows else 0} pieces")


if __name__ == "__main__":
    main()
