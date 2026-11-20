"""Import Hecatomb version 1 results into a Capsomer bundle.

Written against the file layout of a 2022 Hecatomb run. Required files are
``bigtable.tsv`` and ``seqtable.fasta``. Used when present: ``assembly.fasta``,
``contigSeqTable.tsv``, ``contig_count_table.tsv`` and
``contig_with_taxonomy.tsv``.

What Hecatomb's output cannot supply is left empty and reported in the
manifest notes, not guessed.
"""

from dataclasses import dataclass, field
from pathlib import Path

import duckdb

from capsomer import ids
from capsomer.bundle import build_manifest, literal, quote, relation_columns, write_table
from capsomer.importers.common import fasta_records, snake, tsv, values

PIPELINE = "hecatomb"
RANKS = (
    ("kingdom", "superkingdom", "domain_name"),
    ("phylum", "phylum", "phylum_name"),
    ("class", "class", "class_name"),
    ("order", "order", "order_name"),
    ("family", "family", "family_name"),
    ("genus", "genus", "genus_name"),
    ("species", "species", "species_name"),
)
BALTIMORE = ("I", "II", "III", "IV", "V", "VI", "VII")
NOTES = [
    "Imported from Hecatomb version 1 output files.",
    "Read clusters are first reads clustered within each sample. The clustering thresholds are not in the output files, so run_parameter does not list them.",
    "The output holds one final hit and taxon per read cluster. The calls of the separate primary and secondary searches are not in it, so each search type is one run.",
    "Read clusters absent from bigtable.tsv are recorded as searched with no result by both search runs.",
    "Taxa come from lineage names only. The output gives no taxonomy identifiers or release, so taxon identifiers are derived from the lineage path and native_id is empty.",
    "Libraries carry no read totals, so counts cannot be checked against them.",
    "Tool versions and database versions are not in the output files.",
]


@dataclass
class Summary:
    samples: list[str]
    rows: dict[str, int] = field(default_factory=dict[str, int])
    identical_clusters_merged: int = 0


def write_read_clusters(seqtable: Path, samples: set[str], out: Path) -> None:
    """Stream seqtable.fasta, keep the chosen samples, and add identifiers."""
    with out.open("w") as handle:
        handle.write("sequence_id\tname\tsample\tcount\tlength\tgc\tbases\n")
        for name, bases in fasta_records(seqtable):
            sample, count, _ = name.rsplit(":", 2)
            if sample in samples:
                upper = bases.upper()
                handle.write(
                    f"{ids.sequence_id(upper)}\t{name}\t{sample}\t{count}\t{len(upper)}\t"
                    f"{ids.gc_content(upper)}\t{upper}\n"
                )


def write_contigs(assembly: Path, out: Path) -> None:
    with out.open("w") as handle:
        handle.write("sequence_id\tname\tlength\tgc\tbases\n")
        for name, bases in fasta_records(assembly):
            upper = bases.upper()
            handle.write(
                f"{ids.sequence_id(upper)}\t{name.split()[0]}\t{len(upper)}\t"
                f"{ids.gc_content(upper)}\t{upper}\n"
            )


def choose_samples(
    con: duckdb.DuckDBPyConnection, bigtable: Path, samples: list[str] | None, n_samples: int | None
) -> list[str]:
    if samples:
        return sorted(samples)
    found = con.execute(f"SELECT DISTINCT sampleID FROM {tsv(bigtable)} ORDER BY 1").fetchall()
    names = [str(row[0]) for row in found]
    return names[:n_samples] if n_samples else names


def lineage_key(alias: str) -> str:
    parts = ", ".join(f"coalesce({alias}.{quote(col)}, '')" for col, _, _ in RANKS)
    return f"concat_ws('|', {parts})"


def build_taxonomy(con: duckdb.DuckDBPyConnection, sources: list[str]) -> None:
    """Turn lineage name columns into taxon rows with parents."""
    cols = ", ".join(quote(col) for col, _, _ in RANKS)
    union = " UNION ALL ".join(f"SELECT {cols} FROM {source}" for source in sources)
    con.execute(
        f"CREATE TABLE lineage AS SELECT {lineage_key('u')} AS lin_key, * FROM "
        f"(SELECT DISTINCT {cols} FROM ({union})) u "
        f"WHERE coalesce({cols}) IS NOT NULL"
    )
    long = " UNION ALL ".join(
        f"SELECT lin_key, {i} AS idx, {literal(rank)} AS rank_name, {quote(col)} AS name FROM lineage"
        for i, (col, rank, _) in enumerate(RANKS, start=1)
    )
    con.execute(
        "CREATE TABLE lin_path AS SELECT *, 'hecatomb_taxonomy:' || substr(md5(path), 1, 16) AS taxon_id "
        "FROM (SELECT *, string_agg(rank_name || ':' || name, '|') OVER "
        "(PARTITION BY lin_key ORDER BY idx ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) AS path "
        f"FROM ({long}) WHERE name IS NOT NULL)"
    )
    con.execute(
        "CREATE TABLE lin_tax AS SELECT *, "
        "lag(taxon_id) OVER (PARTITION BY lin_key ORDER BY idx) AS parent_taxon_id, "
        "idx = max(idx) OVER (PARTITION BY lin_key) AS is_leaf FROM lin_path"
    )
    con.execute("CREATE TABLE leaf AS SELECT lin_key, taxon_id FROM lin_tax WHERE is_leaf")
    rank_names = ", ".join(
        f"any_value(CASE WHEN t.idx >= {i} THEN l.{quote(col)} END) AS {out}"
        for i, (col, _, out) in enumerate(RANKS, start=1)
    )
    con.execute(
        "CREATE TABLE taxon_rows AS SELECT t.taxon_id, 'hecatomb_taxonomy' AS reference_release_id, "
        "any_value(t.name) AS name, any_value(t.rank_name) AS rank, "
        f"any_value(t.parent_taxon_id) AS parent_taxon_id, any_value(t.idx) AS idx, {rank_names} "
        "FROM lin_tax t JOIN lineage l USING (lin_key) GROUP BY t.taxon_id"
    )


def import_hecatomb(
    hecatomb_dir: Path,
    out_dir: Path,
    *,
    samples: list[str] | None = None,
    n_samples: int | None = None,
    sample_table: Path | None = None,
    sample_id_column: str = "sampleID",
    bundle_id: str | None = None,
    title: str | None = None,
) -> Summary:
    """Build a bundle in out_dir from the Hecatomb results in hecatomb_dir."""
    bigtable, seqtable = hecatomb_dir / "bigtable.tsv", hecatomb_dir / "seqtable.fasta"
    assembly = hecatomb_dir / "assembly.fasta"
    contig_reads = hecatomb_dir / "contigSeqTable.tsv"
    contig_counts = hecatomb_dir / "contig_count_table.tsv"
    contig_taxa = hecatomb_dir / "contig_with_taxonomy.tsv"
    has_contigs = assembly.exists()
    has_contig_reads = has_contigs and contig_reads.exists()
    has_contig_counts = has_contigs and contig_counts.exists()
    has_contig_taxa = has_contigs and contig_taxa.exists()

    out_dir.mkdir(parents=True, exist_ok=True)
    work = out_dir.parent / "work"
    work.mkdir(parents=True, exist_ok=True)
    database = work / "import.duckdb"
    database.unlink(missing_ok=True)
    con = duckdb.connect(str(database))

    chosen = choose_samples(con, bigtable, samples, n_samples)
    summary = Summary(samples=chosen)
    sample_list = ", ".join(literal(s) for s in chosen)

    # ---- read clusters
    write_read_clusters(seqtable, set(chosen), work / "read_clusters.tsv")
    con.execute(
        f"CREATE TABLE rc AS SELECT * FROM read_csv({literal(str(work / 'read_clusters.tsv'))}, "
        "delim='\\t', header=true, quote='', columns={'sequence_id': 'VARCHAR', 'name': 'VARCHAR', "
        "'sample': 'VARCHAR', 'count': 'BIGINT', 'length': 'BIGINT', 'gc': 'DOUBLE', 'bases': 'VARCHAR'})"
    )
    merged = con.execute(
        "SELECT count(*) - count(DISTINCT (sequence_id, sample)) FROM rc"
    ).fetchone()
    summary.identical_clusters_merged = int(merged[0]) if merged else 0
    con.execute(
        "CREATE TABLE seq_reads AS SELECT sequence_id, 'read_cluster' AS kind, any_value(length) AS length, "
        "any_value(bases) AS bases, any_value(gc) AS gc_content FROM rc GROUP BY sequence_id"
    )

    # ---- read-level hits
    con.execute(
        f"CREATE TABLE big AS SELECT b.*, r.sequence_id, {lineage_key('b')} AS lin_key, "
        "CASE b.alnType WHEN 'aa' THEN 'translated_nucleotide_to_protein' ELSE 'nucleotide_to_nucleotide' END AS search_type, "
        f"'{PIPELINE}.search_' || b.alnType AS run_id, b.alnType || ':' || b.targetID AS entry_id "
        f"FROM {tsv(bigtable)} b JOIN (SELECT name, sequence_id FROM rc) r ON b.seqID = r.name"
    )
    lineage_sources = ["big"]

    # ---- contigs
    if has_contigs:
        write_contigs(assembly, work / "contigs.tsv")
        con.execute(
            f"CREATE TABLE ctg AS SELECT * FROM read_csv({literal(str(work / 'contigs.tsv'))}, "
            "delim='\\t', header=true, quote='', max_line_size=1000000000, columns={'sequence_id': 'VARCHAR', 'name': 'VARCHAR', "
            "'length': 'BIGINT', 'gc': 'DOUBLE', 'bases': 'VARCHAR'})"
        )
        con.execute("CREATE TABLE kept (sequence_id VARCHAR)")
    if has_contig_reads:
        con.execute(
            "CREATE TABLE cst AS SELECT DISTINCT r.sequence_id AS query_sequence_id, "
            "c.sequence_id AS target_sequence_id, "
            'least(TRY_CAST(t."start" AS BIGINT), TRY_CAST(t.stop AS BIGINT)) AS target_start, '
            'greatest(TRY_CAST(t."start" AS BIGINT), TRY_CAST(t.stop AS BIGINT)) AS target_end, '
            "TRY_CAST(t.len AS BIGINT) AS alignment_length, TRY_CAST(t.qual AS BIGINT) AS mapping_quality "
            f"FROM {tsv(contig_reads)} t JOIN (SELECT name, sequence_id FROM rc) r ON t.seqID = r.name "
            "JOIN ctg c ON t.contigID = c.name"
        )
        con.execute("INSERT INTO kept SELECT DISTINCT target_sequence_id FROM cst")
    if has_contig_counts:
        con.execute(
            "CREATE TABLE cov AS SELECT c.sequence_id, t.sample_id AS library_id, "
            "max(TRY_CAST(t.reads AS BIGINT)) AS mapped_count, max(TRY_CAST(t.avg_fold_cov AS DOUBLE)) AS mean_depth, "
            "max(TRY_CAST(t.cov_perc AS DOUBLE)) / 100 AS covered_fraction, max(TRY_CAST(t.cov_bases AS BIGINT)) AS covered_bases "
            f"FROM {tsv(contig_counts)} t JOIN ctg c ON t.contig_id = c.name "
            f"WHERE t.sample_id IN ({sample_list}) AND TRY_CAST(t.reads AS BIGINT) > 0 GROUP BY ALL"
        )
        con.execute("INSERT INTO kept SELECT DISTINCT sequence_id FROM cov")
    if has_contigs:
        con.execute(
            "CREATE TABLE seq_contigs AS SELECT sequence_id, 'contig' AS kind, any_value(length) AS length, "
            "any_value(bases) AS bases, any_value(gc) AS gc_content, min(name) AS name, "
            f"'{PIPELINE}.assembly' AS run_id FROM ctg WHERE sequence_id IN (SELECT sequence_id FROM kept) "
            "GROUP BY sequence_id"
        )
    if has_contig_taxa:
        con.execute(
            f"CREATE TABLE ctax AS SELECT t.*, c.sequence_id, {lineage_key('t')} AS lin_key, "
            "'nt:' || t.targetID AS entry_id "
            f"FROM {tsv(contig_taxa)} t JOIN ctg c ON t.contigID = c.name "
            "WHERE c.sequence_id IN (SELECT sequence_id FROM kept)"
        )
        lineage_sources.append("ctax")

    # ---- taxonomy
    build_taxonomy(con, lineage_sources)
    allowed = ", ".join(literal(b) for b in BALTIMORE)
    con.execute(
        "CREATE TABLE baltimore AS SELECT family, any_value(baltimoreType) AS genome_composition, "
        f"any_value(CASE WHEN baltimoreGroup IN ({allowed}) THEN baltimoreGroup END) AS baltimore_class "
        "FROM big WHERE family IS NOT NULL AND baltimoreType IS NOT NULL GROUP BY family"
    )

    # ---- write tables
    def write(name: str, query: str) -> None:
        summary.rows[name] = write_table(con, name, query, out_dir)

    extra_columns: dict[str, dict[str, str]] = {}
    sample_query = values(
        [(s, s, "not_recorded") for s in chosen], ("sample_id", "name", "control_type")
    )
    if sample_table:
        source = tsv(sample_table)
        originals = [c for c in relation_columns(con, source) if c != sample_id_column]
        renamed = {c: snake(c) for c in originals}
        extra_columns["sample"] = {
            new: f"Column '{old}' of the study's sample table." for old, new in renamed.items()
        }
        picks = ", ".join(f"m.{quote(old)} AS {quote(new)}" for old, new in renamed.items())
        sample_query = (
            f"SELECT s.*, {picks} FROM ({sample_query}) s LEFT JOIN {source} m "
            f"ON s.sample_id = m.{quote(sample_id_column)}"
        )
    write("sample", sample_query)
    write(
        "library", values([(s, s, "paired") for s in chosen], ("library_id", "sample_id", "layout"))
    )

    runs: list[tuple[str | None, ...]] = [
        (PIPELINE, "pipeline", None, "A Hecatomb pipeline run."),
        (
            f"{PIPELINE}.read_clustering",
            "read_clustering",
            PIPELINE,
            "Clustering of first reads within each sample.",
        ),
        (
            f"{PIPELINE}.search_aa",
            "search",
            PIPELINE,
            "Translated searches of read clusters against protein databases. The output merges the primary and secondary searches.",
        ),
        (
            f"{PIPELINE}.search_nt",
            "search",
            PIPELINE,
            "Untranslated searches of read clusters against nucleotide databases. The output merges the primary and secondary searches.",
        ),
    ]
    parameters: list[tuple[str | None, ...]] = [
        (f"{PIPELINE}.read_clustering", "read_form", "first_of_pair"),
        (f"{PIPELINE}.read_clustering", "count_unit", "fragments"),
        (f"{PIPELINE}.search_aa", "hits_kept", "best hit per read cluster"),
        (f"{PIPELINE}.search_nt", "hits_kept", "best hit per read cluster"),
    ]
    inputs: list[tuple[str | None, ...]] = [
        (f"{PIPELINE}.search_aa", f"{PIPELINE}.read_clustering", "all read clusters"),
        (
            f"{PIPELINE}.search_nt",
            f"{PIPELINE}.search_aa",
            "read clusters the protein searches did not assign to a virus",
        ),
    ]
    references: list[tuple[str | None, ...]] = [
        (f"{PIPELINE}.search_aa", "hecatomb_protein_db"),
        (f"{PIPELINE}.search_aa", "hecatomb_taxonomy"),
        (f"{PIPELINE}.search_nt", "hecatomb_nucleotide_db"),
        (f"{PIPELINE}.search_nt", "hecatomb_taxonomy"),
    ]
    if has_contigs:
        runs.append((f"{PIPELINE}.assembly", "assembly", PIPELINE, "Assembly of contigs."))
    if has_contig_reads:
        runs.append(
            (
                f"{PIPELINE}.map_read_clusters",
                "read_mapping",
                PIPELINE,
                "Mapping of read clusters to contigs.",
            )
        )
        parameters.append((f"{PIPELINE}.map_read_clusters", "mapped_unit", "read clusters"))
        inputs += [
            (f"{PIPELINE}.map_read_clusters", f"{PIPELINE}.read_clustering", "all read clusters"),
            (f"{PIPELINE}.map_read_clusters", f"{PIPELINE}.assembly", "all contigs"),
        ]
    if has_contig_counts:
        runs.append(
            (
                f"{PIPELINE}.map_reads",
                "read_mapping",
                PIPELINE,
                "Mapping of reads to contigs for coverage.",
            )
        )
        parameters.append((f"{PIPELINE}.map_reads", "mapped_unit", "reads"))
        inputs.append((f"{PIPELINE}.map_reads", f"{PIPELINE}.assembly", "all contigs"))
    if has_contig_taxa:
        runs.append(
            (f"{PIPELINE}.contig_search", "search", PIPELINE, "Nucleotide search of contigs.")
        )
        parameters.append((f"{PIPELINE}.contig_search", "hits_kept", "best hit per contig"))
        inputs.append((f"{PIPELINE}.contig_search", f"{PIPELINE}.assembly", "all contigs"))
        references += [
            (f"{PIPELINE}.contig_search", "hecatomb_nucleotide_db"),
            (f"{PIPELINE}.contig_search", "hecatomb_taxonomy"),
        ]
    write(
        "run",
        f"SELECT run_id, '{PIPELINE}' AS tool, step, parent_run_id, description FROM "
        f"({values(runs, ('run_id', 'step', 'parent_run_id', 'description'))})",
    )
    write("run_parameter", values(parameters, ("run_id", "name", "value")))
    write("run_input", values(inputs, ("run_id", "upstream_run_id", "selection")))
    write("run_reference", values(references, ("run_id", "reference_release_id")))
    write(
        "reference_release",
        values(
            [
                ("hecatomb_protein_db", "Hecatomb protein databases", "sequence_database"),
                ("hecatomb_nucleotide_db", "Hecatomb nucleotide databases", "sequence_database"),
                (
                    "hecatomb_taxonomy",
                    "NCBI taxonomy as packaged with the Hecatomb databases",
                    "formal_taxonomy",
                ),
            ],
            ("reference_release_id", "name", "kind"),
        ),
    )

    entries = "SELECT entry_id, alnType, targetID, targetName AS description, TRY_CAST(tlen AS BIGINT) AS length FROM big"
    if has_contig_taxa:
        entries += (
            " UNION ALL SELECT entry_id, 'nt', targetID, NULL, TRY_CAST(tlen AS BIGINT) FROM ctax"
        )
    write(
        "reference_entry",
        "SELECT entry_id AS reference_entry_id, CASE any_value(alnType) WHEN 'aa' THEN 'hecatomb_protein_db' "
        "ELSE 'hecatomb_nucleotide_db' END AS reference_release_id, any_value(targetID) AS accession, "
        "any_value(description) AS description, max(length) AS length, "
        "CASE any_value(alnType) WHEN 'aa' THEN 'protein' ELSE 'nucleotide' END AS molecule "
        f"FROM ({entries}) GROUP BY entry_id",
    )

    sequences = "SELECT * FROM seq_reads"
    if has_contigs:
        sequences += " UNION ALL BY NAME SELECT * FROM seq_contigs"
    write("sequence", sequences)
    write(
        "taxon",
        "SELECT t.* EXCLUDE (idx), b.genome_composition, b.baltimore_class FROM taxon_rows t "
        "LEFT JOIN baltimore b ON t.idx >= 5 AND t.family_name = b.family",
    )
    write(
        "read_count",
        f"SELECT sequence_id, sample AS library_id, '{PIPELINE}.read_clustering' AS run_id, "
        "'first_of_pair' AS read_form, 'not_applicable' AS read_end, sum(count) AS count, "
        "'fragments' AS count_unit, min(name) AS name FROM rc GROUP BY sequence_id, sample",
    )
    if has_contig_counts:
        write(
            "coverage",
            f"SELECT *, '{PIPELINE}.map_reads' AS run_id, 'reads' AS count_unit FROM cov",
        )

    hit_columns = (
        "least(qs, qe) AS query_start, greatest(qs, qe) AS query_end, "
        "CASE WHEN (qs <= qe) = (ts <= te) THEN 'plus' ELSE 'minus' END AS query_strand, "
        "least(ts, te) AS target_start, greatest(ts, te) AS target_end, "
        "TRY_CAST(pident AS DOUBLE) / 100 AS identity, TRY_CAST(alnlen AS BIGINT) AS alignment_length, "
        "TRY_CAST(mismatches AS BIGINT) AS mismatches, TRY_CAST(evalue AS DOUBLE) AS evalue, "
        "TRY_CAST(bits AS DOUBLE) AS bit_score, TRY_CAST(qcov AS DOUBLE) AS query_coverage, "
        "TRY_CAST(tcov AS DOUBLE) AS target_coverage, true AS is_primary"
    )
    coordinates = (
        "TRY_CAST(qstart AS BIGINT) AS qs, TRY_CAST(qend AS BIGINT) AS qe, "
        "TRY_CAST(tstart AS BIGINT) AS ts, TRY_CAST(tend AS BIGINT) AS te, TRY_CAST(qlen AS BIGINT) AS ql"
    )
    alignments = (
        "SELECT DISTINCT sequence_id AS query_sequence_id, entry_id AS target_reference_entry_id, run_id, "
        "search_type, CASE WHEN alnType = 'aa' THEN (CASE WHEN qs <= qe THEN qs - 1 ELSE ql - qs END) % 3 + 1 END AS query_frame, "
        f"{hit_columns} FROM (SELECT *, {coordinates} FROM big)"
    )
    if has_contig_taxa:
        alignments += (
            " UNION ALL BY NAME SELECT DISTINCT sequence_id AS query_sequence_id, "
            f"entry_id AS target_reference_entry_id, '{PIPELINE}.contig_search' AS run_id, "
            f"'nucleotide_to_nucleotide' AS search_type, {hit_columns} "
            f"FROM (SELECT *, {coordinates} FROM ctax)"
        )
    if has_contig_reads:
        alignments += (
            " UNION ALL BY NAME SELECT query_sequence_id, target_sequence_id, "
            f"'{PIPELINE}.map_read_clusters' AS run_id, 'nucleotide_to_nucleotide' AS search_type, "
            "target_start, target_end, alignment_length, mapping_quality, "
            "row_number() OVER (PARTITION BY query_sequence_id ORDER BY mapping_quality DESC, "
            "alignment_length DESC, target_sequence_id) = 1 AS is_primary FROM cst"
        )
    write("alignment", alignments)

    assignments = (
        "SELECT DISTINCT b.sequence_id, b.run_id, CASE WHEN l.taxon_id IS NULL THEN 'none' ELSE 'assigned' END AS outcome, "
        "l.taxon_id, b.search_type, CASE b.taxMethod WHEN 'TopHit' THEN 'top_hit' WHEN 'LCA' THEN "
        "'lowest_common_ancestor' ELSE 'other' END AS rule, l.taxon_id IS NOT NULL AS is_preferred "
        "FROM big b LEFT JOIN leaf l USING (lin_key) "
        "UNION ALL BY NAME SELECT s.sequence_id, r.run_id, 'none' AS outcome, r.search_type, false AS is_preferred "
        f"FROM seq_reads s CROSS JOIN (VALUES ('{PIPELINE}.search_aa', 'translated_nucleotide_to_protein'), "
        f"('{PIPELINE}.search_nt', 'nucleotide_to_nucleotide')) r(run_id, search_type) "
        "WHERE s.sequence_id NOT IN (SELECT sequence_id FROM big)"
    )
    if has_contig_taxa:
        assignments += (
            f" UNION ALL BY NAME SELECT DISTINCT c.sequence_id, '{PIPELINE}.contig_search' AS run_id, "
            "CASE WHEN l.taxon_id IS NULL THEN 'none' ELSE 'assigned' END AS outcome, l.taxon_id, "
            "'nucleotide_to_nucleotide' AS search_type, 'top_hit' AS rule, l.taxon_id IS NOT NULL AS is_preferred "
            "FROM ctax c LEFT JOIN leaf l USING (lin_key) "
            f"UNION ALL BY NAME SELECT s.sequence_id, '{PIPELINE}.contig_search' AS run_id, 'none' AS outcome, "
            "'nucleotide_to_nucleotide' AS search_type, false AS is_preferred FROM seq_contigs s "
            "WHERE s.sequence_id NOT IN (SELECT sequence_id FROM ctax)"
        )
    write("taxon_assignment", assignments)

    con.close()
    build_manifest(
        out_dir,
        bundle_id=bundle_id or f"hecatomb-{'-'.join(chosen)[:60]}",
        title=title or "Hecatomb results",
        extra_columns=extra_columns,
        notes=NOTES,
    )
    return summary
