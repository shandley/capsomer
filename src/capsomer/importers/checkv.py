"""Import CheckV results: completeness, contamination and quality tiers.

CheckV reports by sequence name. The names are joined to the bundle through
the FASTA file CheckV was given.
"""

from collections.abc import Sequence
from pathlib import Path

from capsomer.builder import Builder
from capsomer.bundle import literal
from capsomer.importers.common import load_fasta, tsv

EXTRA = {
    "expected_length": "Genome length CheckV expects from its closest references, in bases.",
    "n_genes": "Genes CheckV predicted on the sequence.",
    "n_viral_genes": "Genes CheckV matched to viral marker families.",
    "n_host_genes": "Genes CheckV matched to microbial marker families.",
    "is_provirus": "True when CheckV found a host region on the sequence.",
    "proviral_length": "Length of the viral region when CheckV found a host region, in bases.",
    "kmer_freq": "CheckV's k-mer frequency. Values well above 1 suggest the genome is repeated on the sequence.",
    "warnings": "Warnings CheckV attached to the sequence.",
}
NOTE = (
    "CheckV completeness and contamination are stored as fractions. CheckV reports them as percentages. "
    "Host regions CheckV found are recorded by is_provirus and proviral_length only; the trimmed "
    "sequences are not imported."
)


def add_checkv(
    builder: Builder,
    checkv_dir: Path,
    fasta: Path,
    *,
    run_id: str,
    tool_version: str | None = None,
    database_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
) -> int:
    """Add one CheckV run and return the number of quality rows added."""
    con = builder.con
    load_fasta(con, fasta)
    con.execute(
        f"CREATE OR REPLACE TEMP TABLE _quality AS SELECT * FROM {tsv(checkv_dir / 'quality_summary.tsv')}"
    )
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _expected (contig_id VARCHAR, aai_expected_length VARCHAR)"
    )
    con.execute("CREATE OR REPLACE TEMP TABLE _repeat (contig_id VARCHAR, prediction_type VARCHAR)")
    for table, name, columns in (
        ("_expected", "completeness.tsv", "contig_id, aai_expected_length"),
        ("_repeat", "complete_genomes.tsv", "contig_id, prediction_type"),
    ):
        if (checkv_dir / name).exists():
            con.execute(f"INSERT INTO {table} SELECT {columns} FROM {tsv(checkv_dir / name)}")
    missing = builder.scalar(
        "SELECT count(*) FROM _fasta WHERE sequence_id NOT IN (SELECT sequence_id FROM sequence)"
    )
    if missing:
        raise ValueError(f"{missing} sequences of {fasta} are not in the bundle; add them first")

    release = f"checkv_db_{database_version}" if database_version else "checkv_db"
    builder.add_rows(
        "reference_release",
        ("reference_release_id", "name", "version", "kind"),
        [(release, "CheckV database", database_version, "sequence_database")],
    )
    builder.add_run(
        run_id,
        tool="CheckV",
        tool_version=tool_version,
        step="quality_assessment",
        description="Completeness and contamination of viral sequences.",
        parent=parent,
        parameters={
            "completeness_unit": "fraction of the expected genome length",
            "contamination_unit": "fraction of the sequence length that is host",
        },
        inputs=inputs,
        references=[release],
    )
    added = builder.add(
        "quality",
        f"SELECT f.sequence_id, {literal(run_id)} AS run_id, "
        "CAST(q.completeness AS DOUBLE) / 100 AS completeness, q.completeness_method, "
        "CAST(q.contamination AS DOUBLE) / 100 AS contamination, q.checkv_quality AS tier, "
        "CASE q.miuvig_quality WHEN 'High-quality' THEN 'high_quality' "
        "WHEN 'Genome-fragment' THEN 'genome_fragment' END AS miuvig_tier, "
        "CASE r.prediction_type WHEN 'DTR' THEN 'direct' WHEN 'ITR' THEN 'inverted' ELSE 'none' END AS terminal_repeat, "
        "CAST(e.aai_expected_length AS DOUBLE) AS expected_length, "
        "CAST(q.gene_count AS BIGINT) AS n_genes, CAST(q.viral_genes AS BIGINT) AS n_viral_genes, "
        "CAST(q.host_genes AS BIGINT) AS n_host_genes, q.provirus = 'Yes' AS is_provirus, "
        "CAST(q.proviral_length AS BIGINT) AS proviral_length, CAST(q.kmer_freq AS DOUBLE) AS kmer_freq, "
        "nullif(q.warnings, '') AS warnings "
        "FROM _quality q JOIN _fasta f ON f.name = q.contig_id "
        "LEFT JOIN _expected e USING (contig_id) LEFT JOIN _repeat r USING (contig_id) "
        "QUALIFY row_number() OVER (PARTITION BY f.sequence_id ORDER BY q.contig_id) = 1",
    )
    builder.describe("quality", EXTRA)
    builder.note(NOTE)
    return added
