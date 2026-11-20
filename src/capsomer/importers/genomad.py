"""Import geNomad results: viral and plasmid calls, taxonomy, terminal repeats, proviruses.

geNomad reports by sequence name. The names are joined to the bundle through
the FASTA file geNomad was given, so the sequences must already be in the
bundle (see ``add_contigs``).
"""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any, cast

from capsomer import ids
from capsomer.builder import Builder
from capsomer.bundle import literal, quote
from capsomer.importers.common import load_fasta, tsv

RANKS = ("domain", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species")
EXTRA = {
    "viral_call": {
        "n_genes": "Genes geNomad predicted on the sequence.",
        "genetic_code": "Genetic code geNomad chose for gene prediction.",
        "marker_enrichment": "geNomad's marker enrichment for the class it called.",
        "fdr": "geNomad's estimated false discovery rate at this score. Empty unless score calibration was on.",
    },
}
NOTE = (
    "geNomad lists viruses and plasmids only. A sequence it scored and listed as neither has "
    "call none in viral_call, with its three scores."
)


def _find(folder: Path, pattern: str) -> Path | None:
    found = sorted(folder.glob(pattern))
    return found[0] if found and found[0].stat().st_size > 0 else None


def _read(builder: Builder, name: str, path: Path | None, columns: Sequence[str]) -> None:
    """A temp table from a geNomad file, or an empty one with the given columns."""
    if path is None:
        spec = ", ".join(f"{quote(c)} VARCHAR" for c in columns)
        builder.con.execute(f"CREATE OR REPLACE TEMP TABLE {name} ({spec})")
    else:
        builder.con.execute(f"CREATE OR REPLACE TEMP TABLE {name} AS SELECT * FROM {tsv(path)}")


def _parameters(genomad_dir: Path) -> dict[str, object]:
    found: dict[str, object] = {}
    for path in sorted(genomad_dir.glob("*/*.json")):
        data = cast(dict[str, Any], json.loads(path.read_text()))
        module = str(data.get("module", path.stem))
        for key, value in cast(dict[str, Any], data.get("parameters", {})).items():
            found[f"{module}.{key}"] = value
    return found


def add_genomad(
    builder: Builder,
    genomad_dir: Path,
    fasta: Path,
    *,
    run_id: str,
    tool_version: str | None = None,
    database_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
    preferred: bool = False,
) -> dict[str, int]:
    """Add one geNomad run. Returns the number of rows added per table.

    preferred marks the run's calls and taxa as the final answer of its
    pipeline. Leave it off when another step makes the final decision.
    """
    con = builder.con
    load_fasta(con, fasta)
    _read(
        builder,
        "_scores",
        _find(genomad_dir, "*_aggregated_classification/*_aggregated_classification.tsv"),
        ("seq_name", "chromosome_score", "plasmid_score", "virus_score"),
    )
    listed = ("seq_name", "topology", "n_genes", "genetic_code", "fdr", "n_hallmarks")
    _read(
        builder,
        "_virus",
        _find(genomad_dir, "*_summary/*_virus_summary.tsv"),
        (*listed, "coordinates", "virus_score", "marker_enrichment"),
    )
    _read(
        builder,
        "_plasmid",
        _find(genomad_dir, "*_summary/*_plasmid_summary.tsv"),
        (*listed, "plasmid_score", "marker_enrichment"),
    )
    taxonomy = ("seq_name", "n_genes_with_taxonomy", "agreement", "taxid", "lineage")
    _read(builder, "_taxa", _find(genomad_dir, "*_annotate/*_taxonomy.tsv"), taxonomy)
    _read(
        builder,
        "_provirus_taxa",
        _find(genomad_dir, "*_find_proviruses/*_provirus_taxonomy.tsv"),
        taxonomy,
    )
    _read(
        builder,
        "_provirus",
        _find(genomad_dir, "*_find_proviruses/*_provirus.tsv"),
        ("seq_name", "source_seq", "start", "end"),
    )

    release = f"genomad_db_{database_version}" if database_version else "genomad_db"
    builder.add_rows(
        "reference_release",
        ("reference_release_id", "name", "version", "kind"),
        [(release, "geNomad database", database_version, "sequence_database")],
    )
    builder.add_run(
        run_id,
        tool="geNomad",
        tool_version=tool_version,
        step="viral_identification",
        description="Classification of sequences as virus or plasmid, with taxonomy of the viruses.",
        parent=parent,
        parameters=_parameters(genomad_dir),
        inputs=inputs,
        references=[release],
    )
    run = literal(run_id)
    flag = "true" if preferred else "NULL"
    added: dict[str, int] = {}

    # ---- proviruses: each is a new sequence, cut from the sequence it was found on
    regions = con.execute(
        'SELECT p.seq_name, f.sequence_id, CAST(p.start AS BIGINT), CAST(p."end" AS BIGINT), '
        'substr(f.bases, CAST(p.start AS BIGINT), CAST(p."end" AS BIGINT) - CAST(p.start AS BIGINT) + 1) '
        "FROM _provirus p JOIN _fasta f ON f.name = p.source_seq"
    ).fetchall()
    cut = [
        (name, ids.sequence_id(bases), source, start, end, len(bases), ids.gc_content(bases), bases)
        for name, source, start, end, bases in regions
    ]
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _region (name VARCHAR, sequence_id VARCHAR, source_id VARCHAR, "
        "region_start BIGINT, region_end BIGINT, length BIGINT, gc_content DOUBLE, bases VARCHAR)"
    )
    if cut:
        con.executemany("INSERT INTO _region VALUES (?, ?, ?, ?, ?, ?, ?, ?)", cut)
        added["sequence"] = builder.add(
            "sequence",
            "SELECT sequence_id, 'region' AS kind, length, bases, gc_content FROM _region",
        )
        builder.add("sequence_source", f"SELECT sequence_id, {run} AS run_id, name FROM _region")
        added["sequence_relation"] = builder.add(
            "sequence_relation",
            "SELECT DISTINCT sequence_id AS from_sequence_id, source_id AS to_sequence_id, "
            f"'part_of' AS relation, region_start, region_end, 'plus' AS region_strand, {run} AS run_id "
            "FROM _region",
        )
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _subject AS SELECT name, sequence_id FROM _fasta "
        "UNION ALL SELECT name, sequence_id FROM _region"
    )
    missing = builder.scalar(
        "SELECT count(*) FROM _fasta WHERE sequence_id NOT IN (SELECT sequence_id FROM sequence)"
    )
    if missing:
        raise ValueError(f"{missing} sequences of {fasta} are not in the bundle; add them first")

    # ---- viral calls: every sequence geNomad scored or listed
    added["viral_call"] = builder.add(
        "viral_call",
        "WITH named AS (SELECT seq_name FROM _scores UNION SELECT seq_name FROM _virus "
        "UNION SELECT seq_name FROM _plasmid) "
        f"SELECT s.sequence_id, {run} AS run_id, "
        "CASE WHEN v.seq_name IS NOT NULL THEN 'virus' WHEN p.seq_name IS NOT NULL THEN 'plasmid' "
        "ELSE 'none' END AS call, "
        "coalesce(a.virus_score, v.virus_score) AS virus_score, "
        "coalesce(a.plasmid_score, p.plasmid_score) AS plasmid_score, "
        "a.chromosome_score AS cellular_score, "
        "coalesce(v.n_hallmarks, p.n_hallmarks) AS n_hallmark_genes, "
        f"{flag} AS is_preferred, "
        "CAST(coalesce(v.n_genes, p.n_genes) AS BIGINT) AS n_genes, "
        "CAST(coalesce(v.genetic_code, p.genetic_code) AS BIGINT) AS genetic_code, "
        "CAST(coalesce(v.marker_enrichment, p.marker_enrichment) AS DOUBLE) AS marker_enrichment, "
        "CAST(coalesce(v.fdr, p.fdr) AS DOUBLE) AS fdr "
        "FROM named n JOIN _subject s ON s.name = n.seq_name "
        "LEFT JOIN _scores a USING (seq_name) LEFT JOIN _virus v USING (seq_name) "
        "LEFT JOIN _plasmid p USING (seq_name) "
        "QUALIFY row_number() OVER (PARTITION BY s.sequence_id ORDER BY n.seq_name) = 1",
    )
    builder.describe("viral_call", EXTRA["viral_call"])
    builder.note(NOTE)

    # ---- terminal repeats, which geNomad reports for the sequences it lists
    added["quality"] = builder.add(
        "quality",
        f"SELECT s.sequence_id, {run} AS run_id, CASE t.topology WHEN 'No terminal repeats' THEN 'none' "
        "WHEN 'DTR' THEN 'direct' WHEN 'ITR' THEN 'inverted' END AS terminal_repeat "
        "FROM (SELECT seq_name, topology FROM _virus UNION ALL SELECT seq_name, topology FROM _plasmid) t "
        "JOIN _subject s ON s.name = t.seq_name "
        "WHERE t.topology IN ('No terminal repeats', 'DTR', 'ITR') "
        "QUALIFY row_number() OVER (PARTITION BY s.sequence_id ORDER BY t.seq_name) = 1",
    )

    # ---- taxonomy
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _lineage AS SELECT * FROM _taxa "
        "UNION ALL SELECT * FROM _provirus_taxa"
    )
    names = ", ".join(
        f"nullif(string_split(lineage, ';')[{i + 1}], '') AS {rank}_name"
        for i, rank in enumerate(RANKS)
    )
    deepest = "string_split(rtrim(lineage, ';'), ';')"
    rank = (
        "CASE len("
        + deepest
        + ") "
        + " ".join(f"WHEN {i + 1} THEN {literal(r)}" for i, r in enumerate(RANKS) if i)
        + " END"
    )
    added["taxon"] = builder.add(
        "taxon",
        f"SELECT DISTINCT {literal(release)} || ':' || taxid AS taxon_id, {literal(release)} AS reference_release_id, "
        f"taxid AS native_id, {deepest}[-1] AS name, {rank} AS rank, {names} "
        "FROM _lineage WHERE taxid IS NOT NULL AND lineage IS NOT NULL",
    )
    added["taxon_assignment"] = builder.add(
        "taxon_assignment",
        f"SELECT s.sequence_id, {run} AS run_id, "
        "CASE WHEN t.taxid IS NULL THEN 'none' ELSE 'assigned' END AS outcome, "
        f"{literal(release)} || ':' || t.taxid AS taxon_id, 'protein_to_profile' AS search_type, "
        "CASE WHEN t.taxid IS NOT NULL THEN 'majority' END AS rule, "
        "CAST(t.agreement AS DOUBLE) AS score, CASE WHEN t.taxid IS NOT NULL THEN 'agreement' END AS score_type, "
        f"{flag} AS is_preferred, "
        "CASE WHEN t.taxid IS NOT NULL THEN json_object('n_genes_with_taxonomy', CAST(t.n_genes_with_taxonomy AS BIGINT)) END AS detail "
        "FROM (SELECT seq_name FROM _scores UNION SELECT seq_name FROM _lineage) n "
        "JOIN _subject s ON s.name = n.seq_name LEFT JOIN _lineage t USING (seq_name) "
        "QUALIFY row_number() OVER (PARTITION BY s.sequence_id ORDER BY n.seq_name) = 1",
    )
    return added
