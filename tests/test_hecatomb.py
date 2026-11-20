"""Import a tiny Hecatomb-style result set and query the bundle."""

from pathlib import Path

import duckdb

from capsomer import ids
from capsomer.importers.hecatomb import import_hecatomb
from capsomer.validate import conforms, validate

VIRAL = "ACGTACGTACGTAAACGGTT"
BACTERIAL = "TTTTGGGGCCCCAAAATTGG"
DARK = "GGGGGGGGAAAACCCCATAT"
CONTIG = VIRAL + "CCCCCCCCCC" + DARK

BIG_HEADER = [
    "seqID",
    "sampleID",
    "count",
    "CPM",
    "alnType",
    "targetID",
    "evalue",
    "pident",
    "fident",
    "nident",
    "mismatches",
    "qcov",
    "tcov",
    "qstart",
    "qend",
    "qlen",
    "tstart",
    "tend",
    "tlen",
    "alnlen",
    "bits",
    "targetName",
    "taxMethod",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "baltimoreType",
    "baltimoreGroup",
]
VIRUS = [
    "Viruses",
    "Phixviricota",
    "Malgrandaviricetes",
    "Petitvirales",
    "Microviridae",
    "NA",
    "NA",
]


def big_row(
    seq: str, sample: str, aln: str, target: str, lineage: list[str], balt: list[str]
) -> str:
    stats = [
        "1e-10",
        "62.5",
        "0.625",
        "4",
        "2",
        "0.9",
        "0.1",
        "18",
        "1",
        "20",
        "5",
        "10",
        "300",
        "18",
        "40",
    ]
    return "\t".join(
        [seq, sample, "1", "1.0", aln, target, *stats, "a protein", "LCA", *lineage, *balt]
    )


def make_results(folder: Path) -> Path:
    folder.mkdir()
    (folder / "seqtable.fasta").write_text(
        f">s1:3:0\n{VIRAL}\n>s1:1:1\n{BACTERIAL}\n>s2:2:0\n{VIRAL}\n>s2:1:1\n{DARK}\n"
    )
    bacteria = [
        "Bacteria",
        "Bacteroidetes",
        "Bacteroidia",
        "Bacteroidales",
        "Prevotellaceae",
        "Prevotella",
        "NA",
    ]
    rows = [
        big_row("s1:3:0", "s1", "aa", "P1", VIRUS, ["ssDNA", "II"]),
        big_row("s2:2:0", "s2", "aa", "P1", VIRUS, ["ssDNA", "II"]),
        big_row("s1:1:1", "s1", "nt", "N1", bacteria, ["NA", "NA"]),
    ]
    (folder / "bigtable.tsv").write_text("\t".join(BIG_HEADER) + "\n" + "\n".join(rows) + "\n")
    (folder / "assembly.fasta").write_text(
        f">contig_1\n{CONTIG}\n>contig_2\nAAAAACCCCCGGGGGTTTTT\n"
    )
    (folder / "contigSeqTable.tsv").write_text(
        "contigID\tseqID\tstart\tstop\tlen\tqual\n"
        "contig_1\ts1:3:0\t1\t20\t20\t60\ncontig_1\ts2:2:0\t1\t20\t20\t60\ncontig_1\ts2:1:1\t31\t50\t20\t60\n"
    )
    (folder / "contig_count_table.tsv").write_text(
        "sample_id\tcontig_id\tlength\treads\tavg_fold_cov\tcov_perc\tcov_bases\n"
        "s1\tcontig_1\t50\t3\t1.2\t40.0\t20\ns2\tcontig_1\t50\t3\t1.2\t80.0\t40\ns1\tcontig_2\t20\t0\t0\t0\t0\n"
    )
    (folder / "samples.tsv").write_text("sampleID\tTreatment Group\ns1\tvaccine\ns2\tcontrol\n")
    return folder


def test_import_conforms_and_connects(tmp_path: Path) -> None:
    results = make_results(tmp_path / "results")
    bundle = tmp_path / "example" / "bundle"
    summary = import_hecatomb(results, bundle, sample_table=results / "samples.tsv")
    assert summary.samples == ["s1", "s2"]
    findings = validate(bundle)
    assert conforms(findings), findings

    con = duckdb.connect()
    for path in bundle.glob("*.parquet"):
        con.execute(f"CREATE VIEW \"{path.stem}\" AS SELECT * FROM read_parquet('{path}')")

    def one(sql: str) -> object:
        row = con.execute(sql).fetchone()
        assert row is not None
        return row[0]

    viral, dark = ids.sequence_id(VIRAL), ids.sequence_id(DARK)
    # The same read sequence in two libraries is one sequence with two counts.
    assert one(f"SELECT count(*) FROM read_count WHERE sequence_id = '{viral}'") == 2
    assert one("SELECT count(*) FROM sequence WHERE kind = 'read_cluster'") == 3
    # Only the contig that the chosen samples touch is kept.
    assert one("SELECT count(*) FROM sequence WHERE kind = 'contig'") == 1
    # A read cluster with no hit is recorded as searched with no result, by both searches.
    assert (
        one(
            f"SELECT count(*) FROM taxon_assignment WHERE sequence_id = '{dark}' AND outcome = 'none'"
        )
        == 2
    )
    # Connectivity: taxa of read clusters along a contig, by position.
    along = con.execute(
        "SELECT a.target_start, t.family_name FROM alignment a "
        "LEFT JOIN taxon_assignment ta ON ta.sequence_id = a.query_sequence_id AND ta.is_preferred "
        "LEFT JOIN taxon t USING (taxon_id) WHERE a.target_sequence_id IS NOT NULL ORDER BY 1"
    ).fetchall()
    assert along == [(1, "Microviridae"), (31, None)]
    assert one("SELECT baltimore_class FROM taxon WHERE name = 'Microviridae'") == "II"
    assert one("SELECT treatment_group FROM sample WHERE sample_id = 's1'") == "vaccine"
    assert (
        one("SELECT query_strand FROM alignment WHERE search_type LIKE 'translated%' LIMIT 1")
        == "minus"
    )


def test_example_questions_run(tmp_path: Path) -> None:
    from capsomer.bundle import run_sql

    results = make_results(tmp_path / "results")
    bundle = tmp_path / "example" / "bundle"
    import_hecatomb(results, bundle)
    questions = Path(__file__).parents[1] / "examples" / "questions.sql"
    report = run_sql(bundle, questions.read_text())
    assert "Question 12" in report
    assert "unmatched_on_a_contig" in report
