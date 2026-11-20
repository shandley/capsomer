"""Import CheckM2 results: completeness and contamination of genomes.

CheckM2 scores one genome per FASTA file and reports by file name. Each file
becomes a genome in the bundle, made of the sequences in the file, so the
sequences must already be in the bundle.
"""

from collections.abc import Sequence
from pathlib import Path

from capsomer import ids
from capsomer.builder import Builder
from capsomer.bundle import literal
from capsomer.importers.common import load_fasta, tsv

EXTRA = {
    "coding_density": "Fraction of the genome that CheckM2 found to be coding.",
    "translation_table": "Genetic code CheckM2 used for gene prediction.",
    "notes": "Notes CheckM2 attached to the genome.",
}
NOTE = (
    "CheckM2 completeness and contamination are stored as fractions. CheckM2 reports them as "
    "percentages, and its contamination can exceed 100."
)


def add_checkm2(
    builder: Builder,
    report: Path,
    genome_dir: Path,
    *,
    run_id: str,
    extension: str = "fna",
    binning_run: str | None = None,
    tool_version: str | None = None,
    database_version: str | None = None,
    inputs: Sequence[str] = (),
    parent: str | None = None,
) -> int:
    """Add one CheckM2 run and return the number of quality rows added.

    binning_run is the run that formed the genomes, if the bundle records one.
    """
    con = builder.con
    con.execute(f"CREATE OR REPLACE TEMP TABLE _report AS SELECT * FROM {tsv(report)}")
    names = [str(r[0]) for r in con.execute("SELECT Name FROM _report").fetchall()]
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _bin (name VARCHAR, genome_id VARCHAR, sequence_id VARCHAR, length BIGINT)"
    )
    for name in names:
        path = genome_dir / f"{name}.{extension}"
        if not path.exists():
            raise FileNotFoundError(f"no genome file {path} for {name} in {report}")
        load_fasta(con, path, "_one")
        missing = builder.scalar(
            "SELECT count(*) FROM _one WHERE sequence_id NOT IN (SELECT sequence_id FROM sequence)"
        )
        if missing:
            raise ValueError(f"{missing} sequences of {path} are not in the bundle; add them first")
        members = con.execute("SELECT DISTINCT sequence_id, length FROM _one").fetchall()
        genome = ids.genome_id(str(m[0]) for m in members)
        con.executemany(
            "INSERT INTO _bin VALUES (?, ?, ?, ?)", [(name, genome, *m) for m in members]
        )
    known = (
        "SELECT genome_id FROM genome" if "genome" in builder.tables else "SELECT NULL WHERE false"
    )
    con.execute(
        "CREATE OR REPLACE TEMP TABLE _new_genome AS SELECT DISTINCT genome_id FROM _bin "
        f"WHERE genome_id NOT IN ({known})"
    )

    release = f"checkm2_db_{database_version}" if database_version else "checkm2_db"
    builder.add_rows(
        "reference_release",
        ("reference_release_id", "name", "version", "kind"),
        [(release, "CheckM2 database", database_version, "sequence_database")],
    )
    builder.add_run(
        run_id,
        tool="CheckM2",
        tool_version=tool_version,
        step="quality_assessment",
        description="Completeness and contamination of genomes.",
        parent=parent,
        parameters={
            "completeness_unit": "fraction of the expected gene content",
            "contamination_unit": "fraction, from gene content that should not be present",
        },
        inputs=inputs,
        references=[release],
    )
    former = literal(binning_run) if binning_run else "NULL"
    builder.add(
        "genome",
        "SELECT genome_id, CASE WHEN count(*) > 1 THEN 'bin' ELSE 'single_sequence' END AS kind, "
        f"CASE WHEN count(*) > 1 THEN {former} END AS run_id, count(*) AS n_sequences, "
        "sum(length) AS total_length FROM (SELECT DISTINCT genome_id, sequence_id, length FROM _bin) "
        "WHERE genome_id IN (SELECT genome_id FROM _new_genome) GROUP BY genome_id",
    )
    builder.add(
        "genome_member",
        f"SELECT DISTINCT genome_id, sequence_id, {former} AS run_id FROM _bin "
        "WHERE genome_id IN (SELECT genome_id FROM _new_genome)",
    )
    added = builder.add(
        "quality",
        f"SELECT b.genome_id, {literal(run_id)} AS run_id, round(CAST(r.Completeness AS DOUBLE) / 100, 6) AS completeness, "
        "r.Completeness_Model_Used AS completeness_method, round(CAST(r.Contamination AS DOUBLE) / 100, 6) AS contamination, "
        "CAST(r.Coding_Density AS DOUBLE) AS coding_density, CAST(r.Translation_Table_Used AS BIGINT) AS translation_table, "
        "nullif(r.Additional_Notes, 'None') AS notes "
        "FROM _report r JOIN (SELECT DISTINCT name, genome_id FROM _bin) b ON b.name = r.Name "
        "QUALIFY row_number() OVER (PARTITION BY b.genome_id ORDER BY r.Name) = 1",
    )
    builder.describe("quality", EXTRA)
    builder.note(NOTE)
    return added
