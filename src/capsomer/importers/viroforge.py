"""Import the ground truth of a ViroForge synthetic dataset into a Capsomer bundle.

ViroForge simulates virome sequencing from known genomes and records what is
true: which genome every read came from, the taxonomy of each genome, the
contaminants, the injected artifacts, and which viral genomes have no family
in its taxonomy (the stratum its benchmark calls dark). The bundle built here
is that truth, so it does not depend on any analysis pipeline.

Truth with no catalog table yet goes in two tables outside the catalog,
``x_sequence_origin`` and ``x_read_artifact``.
"""

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

import duckdb

from capsomer import ids
from capsomer.bundle import build_manifest, literal, write_table
from capsomer.importers.common import fasta_records, values

RUN = "viroforge"
RANKS = ("realm", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus", "species")
RANK_COLUMNS = {
    "realm": "realm_name",
    "kingdom": "kingdom_name",
    "phylum": "phylum_name",
    "class": "class_name",
    "order": "order_name",
    "family": "family_name",
    "genus": "genus_name",
    "species": "species_name",
}
READ_NAME = re.compile(r"^(?P<genome>.+)_\d+_\d+/(?P<end>[12])(?P<dup>_dup\d+)?$")
SETTINGS = ("coverage", "read_length", "insert_size", "contamination_level", "platform")
NOTES = [
    "Ground truth of a ViroForge synthetic dataset. Nothing here was inferred by an analysis pipeline.",
    "Each alignment row states the genome a read was simulated from. It carries no coordinates.",
    "Ancestor taxa have identifiers derived from their lineage path. Only the genome-level taxon carries an NCBI taxonomy identifier.",
    "x_sequence_origin and x_read_artifact hold truth that the version 1 catalog has no table for.",
]
CUSTOM_TABLES: dict[str, dict[str, Any]] = {
    "x_sequence_origin": {
        "shape": "label",
        "one_row_is": "the true origin of one simulated record of a reference sequence",
        "description": "Whether each simulated sequence is a virus or a contaminant, and whether a viral genome has a family in the simulator's taxonomy.",
        "key": ["sequence_id", "run_id", "simulator_name"],
    },
    "x_read_artifact": {
        "shape": "label",
        "one_row_is": "one artifact injected into one read",
        "description": "PCR duplicates, adapter read-through and low-complexity reads that the simulator injected.",
        "key": ["read_id", "run_id", "artifact"],
    },
}
EXTRA_COLUMNS = {
    "sample": {"collection_description": "ViroForge's description of the simulated collection."},
    "coverage": {
        "expected_depth": "Depth the simulator aimed for on this sequence.",
        "designed_relative_abundance": "Relative abundance the simulator assigned, over all sequences.",
    },
    "x_sequence_origin": {
        "sequence_id": "The reference sequence.",
        "run_id": "The simulation run.",
        "origin": "viral, or the kind of contaminant.",
        "is_known": "For viral genomes: false when the simulator's taxonomy has no family for the genome. ViroForge's benchmark scores those as its dark stratum.",
        "organism": "Organism or source of a contaminant.",
        "simulator_name": "The name the simulator gave this record. Several records can share one sequence.",
    },
    "x_read_artifact": {
        "read_id": "The read.",
        "run_id": "The simulation run.",
        "artifact": "pcr_duplicate, adapter or low_complexity.",
        "detail": "The adapter injected, the kind of low-complexity sequence, or the read a duplicate copies.",
    },
}


@dataclass
class Summary:
    sample: str
    rows: dict[str, int] = field(default_factory=dict[str, int])
    reads_without_source: int = 0


def one_file(folder: Path, pattern: str) -> Path:
    found = sorted(folder.glob(pattern))
    if len(found) != 1:
        raise FileNotFoundError(
            f"expected one file matching {pattern} in {folder}, found {len(found)}"
        )
    return found[0]


def write_reads(fastq: Path, out: Any) -> None:
    """One line per read: name, end, source genome, and any artifact labels."""
    with fastq.open() as handle:
        for number, line in enumerate(handle):
            if number % 4:
                continue
            name, *labels = line[1:].split()
            attributes = dict(label.split("=", 1) for label in labels if "=" in label)
            match = READ_NAME.match(name)
            genome = match["genome"] if match else ""
            end = "first" if name.rsplit("/", 1)[-1].startswith("1") else "second"
            out.write(
                "\t".join(
                    [
                        name,
                        end,
                        genome,
                        "1" if match and match["dup"] else "0",
                        attributes.get("adapter_injected", ""),
                        attributes.get("artifact_type", ""),
                    ]
                )
                + "\n"
            )


def taxon_rows(
    taxonomy: dict[str, dict[str, Any]], names: dict[str, str]
) -> tuple[list[dict[str, str | None]], dict[str, str]]:
    """Taxon rows for every genome and its ancestors, and genome -> taxon identifier."""
    rows: dict[str, dict[str, str | None]] = {}
    leaf: dict[str, str] = {}
    for genome, entry in taxonomy.items():
        parent: str | None = None
        path: list[str] = []
        lineage: dict[str, str | None] = {"domain_name": "Viruses"}
        for rank in RANKS:
            name = entry.get(rank)
            if not name or name == "Unknown":
                continue
            path.append(f"{rank}:{name}")
            taxon_id = "viroforge_taxonomy:" + hashlib.md5("|".join(path).encode()).hexdigest()[:16]
            if rank in RANK_COLUMNS:
                lineage[RANK_COLUMNS[rank]] = str(name)
            rows.setdefault(
                taxon_id,
                {
                    "taxon_id": taxon_id,
                    "name": str(name),
                    "rank": rank,
                    "parent_taxon_id": parent,
                    **lineage,
                },
            )
            parent = taxon_id
        taxid = entry.get("ncbi_taxid")
        taxon_id = (
            f"viroforge_taxonomy:ncbi:{taxid}" if taxid else f"viroforge_taxonomy:genome:{genome}"
        )
        rows.setdefault(
            taxon_id,
            {
                "taxon_id": taxon_id,
                "native_id": str(taxid) if taxid else None,
                "name": names.get(genome, genome),
                "parent_taxon_id": parent,
                **lineage,
            },
        )
        leaf[genome] = taxon_id
    return list(rows.values()), leaf


def dict_rows(rows: list[dict[str, str | None]], columns: tuple[str, ...]) -> str:
    return values([tuple(row.get(c) for c in columns) for row in rows], columns)


def import_viroforge(dataset_dir: Path, out_dir: Path, *, sample_id: str | None = None) -> Summary:
    """Build a truth bundle in out_dir from a ViroForge dataset folder."""
    metadata = cast(
        dict[str, Any],
        json.loads(one_file(dataset_dir / "metadata", "*_metadata.json").read_text()),
    )
    config = cast(dict[str, Any], metadata["configuration"])
    bench = cast(dict[str, Any], metadata["benchmarking"])
    sequences = cast(list[dict[str, Any]], metadata["sequences"])
    sample = sample_id or dataset_dir.name
    summary = Summary(sample=sample)

    out_dir.mkdir(parents=True, exist_ok=True)
    work = out_dir.parent / "work"
    work.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()

    # ---- reference sequences
    by_genome: dict[str, str] = {}
    con.execute(
        "CREATE TABLE seq (sequence_id VARCHAR, name VARCHAR, length BIGINT, gc_content DOUBLE, bases VARCHAR)"
    )
    for header, bases in fasta_records(one_file(dataset_dir / "fasta", "*.fasta")):
        genome, upper = header.split()[0], bases.upper()
        by_genome[genome] = ids.sequence_id(upper)
        con.execute(
            "INSERT INTO seq VALUES (?, ?, ?, ?, ?)",
            [by_genome[genome], genome, len(upper), ids.gc_content(upper), upper],
        )

    # ---- reads and their true source
    with (work / "reads.tsv").open("w") as handle:
        handle.write("read_id\tread_end\tgenome\tis_duplicate\tadapter\tartifact_type\n")
        write_reads(one_file(dataset_dir / "fastq", "*_R1.fastq"), handle)
        write_reads(one_file(dataset_dir / "fastq", "*_R2.fastq"), handle)
    con.execute(
        f"CREATE TABLE reads AS SELECT r.*, s.sequence_id FROM read_csv({literal(str(work / 'reads.tsv'))}, "
        "delim='\\t', header=true, quote='', all_varchar=true) r LEFT JOIN seq s ON r.genome = s.name"
    )
    lost = con.execute("SELECT count(*) FROM reads WHERE sequence_id IS NULL").fetchone()
    summary.reads_without_source = int(lost[0]) if lost else 0
    totals = con.execute(
        "SELECT count(*), count(*) FILTER (WHERE read_end = 'first') FROM reads"
    ).fetchone()
    n_reads, n_fragments = (int(totals[0]), int(totals[1])) if totals else (0, 0)

    def write(name: str, query: str) -> None:
        summary.rows[name] = write_table(con, name, query, out_dir)

    # ---- context
    collection = cast(dict[str, Any], metadata["collection"])
    write(
        "sample",
        values(
            [
                (
                    sample,
                    str(collection["name"]),
                    "not_control",
                    str(collection.get("description", "")),
                )
            ],
            ("sample_id", "name", "control_type", "collection_description"),
        ),
    )
    write(
        "library",
        f"SELECT {literal(sample)} AS library_id, {literal(sample)} AS sample_id, "
        f"{literal(str(config.get('molecule_type', '')).upper())} AS nucleic_acid, "
        f"{literal(str(config.get('vlp_protocol')))} AS enrichment, "
        f"{literal(str(config.get('amplification')))} AS amplification, "
        f"{literal(str(config.get('platform')))} AS platform, 'paired' AS layout, "
        f"{n_reads} AS reads_raw, {n_fragments} AS fragments_raw",
    )
    generation = cast(dict[str, Any], metadata["generation_info"])
    write(
        "run",
        values(
            [
                (
                    RUN,
                    "viroforge",
                    str(generation.get("viroforge_version")),
                    "simulation",
                    "Simulation of a virome sequencing library with known truth.",
                )
            ],
            ("run_id", "tool", "tool_version", "step", "description"),
        ),
    )
    settings: list[tuple[str | None, ...]] = [
        (RUN, "random_seed", str(generation.get("random_seed")))
    ]
    settings += [(RUN, key, str(config[key])) for key in SETTINGS if config.get(key) is not None]
    dark = cast(dict[str, Any], config.get("dark_matter") or {})
    if dark:
        settings.append((RUN, "dark_matter_fraction", str(dark.get("dark_matter_fraction"))))
    for block, key in (
        ("duplicate_stats", "duplicate_rate"),
        ("adapter_stats", "adapter_rate"),
        ("low_complexity_stats", "rate"),
    ):
        if key in metadata.get(block, {}):
            settings.append(
                (RUN, f"{block.removesuffix('_stats')}_rate", str(metadata[block][key]))
            )
    settings.append(
        (RUN, "relative_abundance_denominator", "all simulated sequences, viral and contaminant")
    )
    write("run_parameter", values(settings, ("run_id", "name", "value")))
    write(
        "reference_release",
        values(
            [
                (
                    "viroforge_genomes",
                    "ViroForge genome collection: NCBI RefSeq viral genomes and contaminant references",
                    "sequence_database",
                ),
                (
                    "viroforge_taxonomy",
                    "ICTV ranks with NCBI taxonomy identifiers, as held in the ViroForge database",
                    "formal_taxonomy",
                ),
            ],
            ("reference_release_id", "name", "kind"),
        ),
    )
    write(
        "run_reference",
        values(
            [(RUN, "viroforge_genomes"), (RUN, "viroforge_taxonomy")],
            ("run_id", "reference_release_id"),
        ),
    )

    # ---- sequences, genomes, external identifiers
    # Records with identical bases are one sequence. The simulator's names for it stay in x_sequence_origin.
    write(
        "sequence",
        "SELECT sequence_id, 'reference' AS kind, any_value(length) AS length, any_value(bases) AS bases, "
        "any_value(gc_content) AS gc_content, min(name) AS name FROM seq GROUP BY sequence_id",
    )
    genomes = [(ids.genome_id([seq_id]), seq_id, name) for name, seq_id in by_genome.items()]
    con.execute("CREATE TABLE gen (genome_id VARCHAR, sequence_id VARCHAR, name VARCHAR)")
    con.executemany("INSERT INTO gen VALUES (?, ?, ?)", genomes)
    write(
        "genome",
        "SELECT DISTINCT genome_id, 'single_sequence' AS kind, 1 AS n_sequences, s.length AS total_length FROM gen JOIN seq s USING (sequence_id)",
    )
    write("genome_member", "SELECT DISTINCT genome_id, sequence_id FROM gen")
    write(
        "identifier_map",
        "SELECT DISTINCT sequence_id, 'NCBI assembly' AS authority, name AS external_id FROM gen "
        "WHERE name LIKE 'GCF_%' OR name LIKE 'GCA_%'",
    )

    # ---- taxonomy and truth labels
    names = {str(s["genome_id"]): str(s["genome_name"]) for s in sequences}
    taxa, leaf = taxon_rows(cast(dict[str, dict[str, Any]], bench.get("taxonomy", {})), names)
    columns = (
        "taxon_id",
        "native_id",
        "name",
        "rank",
        "parent_taxon_id",
        "domain_name",
        *RANK_COLUMNS.values(),
    )
    write(
        "taxon",
        f"SELECT *, 'viroforge_taxonomy' AS reference_release_id FROM ({dict_rows(taxa, columns)})",
    )
    con.execute("CREATE TABLE leaf (name VARCHAR, taxon_id VARCHAR)")
    con.executemany("INSERT INTO leaf VALUES (?, ?)", list(leaf.items()))
    write(
        "taxon_assignment",
        f"SELECT DISTINCT g.genome_id, {literal(RUN)} AS run_id, 'assigned' AS outcome, l.taxon_id, "
        "'other' AS rule, true AS is_preferred FROM gen g JOIN leaf l USING (name)",
    )

    # ---- counts
    con.execute(
        "CREATE TABLE per_seq AS SELECT sequence_id, count(*) AS fragments FROM reads "
        "WHERE read_end = 'first' AND sequence_id IS NOT NULL GROUP BY 1"
    )
    expected = cast(list[dict[str, Any]], bench.get("expected_coverage", {}).get("per_genome", []))
    con.execute(
        "CREATE TABLE design (name VARCHAR, expected_depth DOUBLE, designed_relative_abundance DOUBLE)"
    )
    con.executemany(
        "INSERT INTO design VALUES (?, ?, ?)",
        [
            (str(e["genome_id"]), e.get("expected_coverage"), e.get("relative_abundance"))
            for e in expected
        ],
    )
    write(
        "coverage",
        f"SELECT p.sequence_id, {literal(sample)} AS library_id, {literal(RUN)} AS run_id, p.fragments AS mapped_count, "
        "'fragments' AS count_unit, sum(d.expected_depth) AS expected_depth, "
        "sum(d.designed_relative_abundance) AS designed_relative_abundance "
        "FROM per_seq p JOIN gen g USING (sequence_id) LEFT JOIN design d USING (name) GROUP BY ALL",
    )
    write(
        "taxon_profile",
        f"SELECT l.taxon_id, {literal(sample)} AS library_id, {literal(RUN)} AS run_id, sum(p.fragments) AS direct_count, "
        "'fragments' AS count_unit, sum(d.designed_relative_abundance) AS relative_abundance "
        "FROM (SELECT DISTINCT sequence_id, name FROM gen) g JOIN leaf l USING (name) JOIN per_seq p USING (sequence_id) "
        "LEFT JOIN design d USING (name) GROUP BY ALL",
    )

    # ---- reads and their true origin
    write(
        "read",
        f"SELECT read_id, {literal(sample)} AS library_id, read_id AS name, read_end FROM reads",
    )
    write(
        "alignment",
        f"SELECT read_id AS query_read_id, sequence_id AS target_sequence_id, {literal(RUN)} AS run_id, "
        "'other' AS search_type, true AS is_primary FROM reads WHERE sequence_id IS NOT NULL",
    )

    # ---- truth outside the catalog
    contaminants = cast(
        list[dict[str, Any]], bench.get("contamination_manifest", {}).get("contaminants", [])
    )
    known = {
        g: bool(e.get("is_known", True))
        for g, e in cast(dict[str, dict[str, Any]], bench.get("taxonomy", {})).items()
    }
    origin: list[tuple[str, str, bool | None, str | None]] = [
        (
            str(c["genome_id"]),
            str(c.get("contaminant_type", "contaminant")),
            None,
            str(c.get("organism", "")) or None,
        )
        for c in contaminants
    ]
    origin += [
        (str(s["genome_id"]), "viral", known.get(str(s["genome_id"])), None)
        for s in sequences
        if s.get("sequence_type") == "viral"
    ]
    con.execute(
        "CREATE TABLE origin (name VARCHAR, origin VARCHAR, is_known BOOLEAN, organism VARCHAR)"
    )
    con.executemany("INSERT INTO origin VALUES (?, ?, ?, ?)", origin)
    con.execute(
        f"COPY (SELECT DISTINCT g.sequence_id, {literal(RUN)} AS run_id, o.origin, o.is_known, o.organism, o.name AS simulator_name "
        f"FROM origin o JOIN gen g USING (name) ORDER BY ALL) TO {literal(str(out_dir / 'x_sequence_origin.parquet'))} "
        "(FORMAT parquet, COMPRESSION zstd)"
    )
    duplicates = sorted((dataset_dir / "metadata").glob("*_duplicate_manifest.tsv"))
    template = (
        f"(SELECT copy_read_id, template_read_id FROM read_csv({literal(str(duplicates[0]))}, delim='\\t', header=true, all_varchar=true))"
        if duplicates
        else "(SELECT NULL AS copy_read_id, NULL AS template_read_id WHERE false)"
    )
    con.execute(
        f"COPY (SELECT read_id, {literal(RUN)} AS run_id, 'pcr_duplicate' AS artifact, 'copy of ' || t.template_read_id AS detail "
        f"FROM reads r LEFT JOIN {template} t ON r.read_id = t.copy_read_id WHERE r.is_duplicate = '1' "
        f"UNION ALL SELECT read_id, {literal(RUN)}, 'adapter', adapter FROM reads WHERE adapter <> '' "
        f"UNION ALL SELECT read_id, {literal(RUN)}, 'low_complexity', artifact_type FROM reads WHERE artifact_type <> '' "
        f"ORDER BY ALL) TO {literal(str(out_dir / 'x_read_artifact.parquet'))} (FORMAT parquet, COMPRESSION zstd)"
    )
    con.close()
    build_manifest(
        out_dir,
        bundle_id=f"viroforge-{sample}",
        title=f"ViroForge ground truth: {collection['name']}",
        extra_columns=EXTRA_COLUMNS,
        notes=NOTES,
        custom_tables=CUSTOM_TABLES,
    )
    return summary
