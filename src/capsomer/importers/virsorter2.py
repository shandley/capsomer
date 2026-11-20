"""Import VirSorter2 results: viral calls with scores, and proviruses.

VirSorter2 reports by sequence name. The names are joined to the bundle through
the FASTA file VirSorter2 was given.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Any, cast

import yaml

from capsomer import ids
from capsomer.builder import Builder
from capsomer.bundle import literal, quote
from capsomer.importers.common import fasta_records, load_fasta, tsv

SETTINGS = (
    "GROUPS",
    "MIN_LENGTH",
    "PROBA_CUTOFF",
    "HALLMARK_REQUIRED",
    "HALLMARK_REQUIRED_ON_SHORT",
    "VIRAL_GENE_REQUIRED",
    "VIRAL_GENE_ENRICH_OFF",
    "HIGH_CONFIDENCE_ONLY",
    "PROVIRUS",
    "KEEP_ORIGINAL_SEQ",
)
FIXED = ("seqname", "max_score", "max_score_group", "length", "hallmark", "viral", "cellular")
EXTRA = {
    "best_group": "The VirSorter2 viral group with the highest score.",
    "viral_gene_percent": "Percentage of the sequence's genes that VirSorter2 annotated as viral.",
    "cellular_gene_percent": "Percentage of the sequence's genes that VirSorter2 annotated as cellular.",
    "fewer_than_two_genes": "True when VirSorter2 kept the sequence although it has fewer than two genes.",
}
NOTE = (
    "VirSorter2 lists the sequences it calls viral. A sequence it was given and did not list has "
    "call none in viral_call, with a score when VirSorter2 kept one. VirSorter2 trims the ends of "
    "the sequences it writes out; calls are attached to the untrimmed input sequence."
)


def _number(column: str) -> str:
    """SQL for a column as a number. VirSorter2 writes nan where it has no score."""
    value = f"TRY_CAST({column} AS DOUBLE)"
    return f"CASE WHEN isnan({value}) THEN NULL ELSE {value} END"


def _first(folder: Path, pattern: str) -> Path | None:
    found = sorted(folder.glob(pattern))
    return found[0] if found and found[0].stat().st_size > 0 else None


def _settings(results: Path) -> dict[str, object]:
    path = results / "config.yaml"
    if not path.exists():
        return {}
    config = cast(dict[str, Any], yaml.safe_load(path.read_text()) or {})
    found: dict[str, object] = {}
    for key in SETTINGS:
        value = config.get(key)
        if value is not None:
            listed = cast(list[object], value) if isinstance(value, list) else None
            found[key.lower()] = ",".join(str(v) for v in listed) if listed else value
    return found


def add_virsorter2(
    builder: Builder,
    results: Path,
    fasta: Path,
    *,
    run_id: str,
    tool_version: str | None = None,
    database_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    preferred: bool = False,
) -> dict[str, int]:
    """Add one VirSorter2 run. Returns the number of rows added per table."""
    con = builder.con
    load_fasta(con, fasta)
    missing = builder.scalar(
        "SELECT count(*) FROM _fasta WHERE sequence_id NOT IN (SELECT sequence_id FROM sequence)"
    )
    if missing:
        raise ValueError(f"{missing} sequences of {fasta} are not in the bundle; add them first")
    scores = _first(results, "*final-viral-score.tsv")
    if scores is None:
        raise FileNotFoundError(f"no final-viral-score.tsv in {results}")
    con.execute(f"CREATE OR REPLACE TEMP TABLE _called AS SELECT * FROM {tsv(scores)}")
    columns = [str(r[0]) for r in con.execute("DESCRIBE _called").fetchall()]
    groups = [c for c in columns if c not in FIXED]
    every = _first(results, "*/all-fullseq-proba.tsv")
    if every:
        con.execute(f"CREATE OR REPLACE TEMP TABLE _every AS SELECT * FROM {tsv(every)}")
        found = [str(r[0]) for r in con.execute("DESCRIBE _every").fetchall() if r[0] != "seqname"]
        best = "greatest(" + ", ".join(_number(quote(g)) for g in found) + ")"
        con.execute(
            f"CREATE OR REPLACE TEMP TABLE _every AS SELECT seqname, {best} AS score FROM _every"
        )
    else:
        con.execute("CREATE OR REPLACE TEMP TABLE _every (seqname VARCHAR, score DOUBLE)")

    release = f"virsorter2_db_{database_version}" if database_version else "virsorter2_db"
    builder.add_rows(
        "reference_release",
        ("reference_release_id", "name", "version", "kind"),
        [(release, "VirSorter2 database", database_version, "sequence_database")],
    )
    builder.add_run(
        run_id,
        tool="VirSorter2",
        tool_version=tool_version,
        step="viral_identification",
        description="Classification of sequences as viral, by viral group.",
        parent=parent,
        parameters=_settings(results),
        inputs=inputs,
        references=[release],
    )
    run = literal(run_id)
    added: dict[str, int] = {}

    # ---- partial calls: each is a stretch of an input sequence, taken from VirSorter2's own output
    partial = [
        str(r[0])
        for r in con.execute(
            "SELECT seqname FROM _called WHERE seqname LIKE '%||%_partial'"
        ).fetchall()
    ]
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _region (name VARCHAR, sequence_id VARCHAR, source_id VARCHAR, "
        "region_start BIGINT, region_end BIGINT, length BIGINT, gc_content DOUBLE, bases VARCHAR)"
    )
    written = _first(results, "*final-viral-combined.fa")
    if partial and written:
        wanted = set(partial)
        for header, bases in fasta_records(written):
            name = header.split()[0]
            if name not in wanted:
                continue
            upper = bases.upper()
            source = con.execute(
                "SELECT sequence_id, bases FROM _fasta WHERE name = ?", [name.split("||")[0]]
            ).fetchone()
            if source is None:
                continue
            at = str(source[1]).find(upper)
            con.execute(
                "INSERT INTO _region VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    name,
                    ids.sequence_id(upper),
                    source[0],
                    at + 1 if at >= 0 else None,
                    at + len(upper) if at >= 0 else None,
                    len(upper),
                    ids.gc_content(upper),
                    upper,
                ],
            )
        added["sequence"] = builder.add(
            "sequence",
            "SELECT sequence_id, 'region' AS kind, length, bases, gc_content FROM _region",
        )
        builder.add("sequence_source", f"SELECT sequence_id, {run} AS run_id, name FROM _region")
        added["sequence_relation"] = builder.add(
            "sequence_relation",
            "SELECT DISTINCT sequence_id AS from_sequence_id, source_id AS to_sequence_id, "
            f"'part_of' AS relation, region_start, region_end, "
            f"CASE WHEN region_start IS NOT NULL THEN 'plus' END AS region_strand, {run} AS run_id "
            "FROM _region WHERE sequence_id <> source_id",
        )

    # ---- one call per sequence: listed sequences are viral, the rest of the input is not
    detail = (
        "json_object("
        + ", ".join(f"{literal(g)}, {_number('c.' + quote(g))}" for g in groups)
        + ")"
    )
    flag = "true" if preferred else "NULL"
    added["viral_call"] = builder.add(
        "viral_call",
        "WITH listed AS (SELECT c.*, split_part(c.seqname, '||', 1) AS base, "
        "split_part(c.seqname, '||', 2) AS suffix FROM _called c), "
        "subject AS ("
        " SELECT l.seqname, coalesce(r.sequence_id, f.sequence_id) AS sequence_id FROM listed l "
        " LEFT JOIN _region r ON r.name = l.seqname LEFT JOIN _fasta f ON f.name = l.base "
        " WHERE r.sequence_id IS NOT NULL OR l.suffix NOT LIKE '%_partial') "
        f"SELECT s.sequence_id, {run} AS run_id, 'virus' AS call, {_number('c.max_score')} AS virus_score, "
        f"CAST(c.hallmark AS BIGINT) AS n_hallmark_genes, {flag} AS is_preferred, {detail} AS detail, "
        f"c.max_score_group AS best_group, {_number('c.viral')} AS viral_gene_percent, "
        f"{_number('c.cellular')} AS cellular_gene_percent, "
        "split_part(c.seqname, '||', 2) = 'lt2gene' AS fewer_than_two_genes "
        "FROM subject s JOIN _called c USING (seqname) "
        f"QUALIFY row_number() OVER (PARTITION BY s.sequence_id ORDER BY {_number('c.max_score')} DESC NULLS LAST) = 1",
    )
    added["viral_call"] += builder.add(
        "viral_call",
        f"SELECT f.sequence_id, {run} AS run_id, 'none' AS call, max(e.score) AS virus_score, "
        f"{flag} AS is_preferred FROM _fasta f LEFT JOIN _every e ON e.seqname = f.name "
        f"WHERE f.sequence_id NOT IN (SELECT sequence_id FROM viral_call WHERE run_id = {run}) "
        "GROUP BY f.sequence_id",
    )
    builder.describe("viral_call", EXTRA)
    builder.note(NOTE)
    return added
