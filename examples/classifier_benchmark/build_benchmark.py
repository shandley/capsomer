"""Put the test pieces, what is known about each, and every classifier's calls into one Capsomer dataset.

The dataset holds the pieces as sequences, one viral_call row per piece per
classifier, and a table outside the catalog, x_truth, with the known origin of
each piece. Scoring is then a query (see benchmark.sql).

Usage: build_benchmark.py <work folder> <bundle folder>
The work folder holds bench/ (simulated viruses and real non-viral pieces),
recent_pieces/ and recent/ (genomes deposited after every tool's database
date) and novelty/ (the similarity of the recent genomes to RefSeq viruses).
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

from capsomer.builder import Builder
from capsomer.bundle import literal
from capsomer.importers.classifiers import add_cenotetaker3, add_jaeger, add_phabox2
from capsomer.importers.common import load_fasta
from capsomer.importers.genomad import add_genomad
from capsomer.importers.study import init_bundle
from capsomer.importers.virsorter2 import add_virsorter2

TRUTH_COLUMNS = {
    "sequence_id": "The piece.",
    "run_id": "The run that cut the pieces.",
    "test_set": "simulated_refseq: a RefSeq virus used by ViroForge. recent_cultured and recent_uncultivated: "
    "a viral genome deposited at NCBI after the database date of every classifier compared. non_viral: a real "
    "bacterial, archaeal, fungal or human fragment.",
    "is_viral": "True when the piece comes from a virus.",
    "origin": "viral, or the kind of non-viral sequence.",
    "source": "Accession or name of the sequence the piece was cut from.",
    "has_family": "Simulated viruses only: whether ViroForge's taxonomy has a family for the genome.",
    "refseq_protein_coverage": "Recent genomes only: fraction of the piece covered by translated matches "
    "(e-value 1e-5 or better) to any RefSeq virus. Low values mean a virus unlike those in RefSeq.",
    "refseq_close_coverage": "Recent genomes only: the same fraction, counting only matches with at least "
    "50% amino acid identity.",
    "study": "Recent genomes only: the naming pattern of the deposit, with digits removed.",
}


def union_length(intervals: list[tuple[int, int]]) -> int:
    total, reach = 0, -1
    for start, end in sorted(intervals):
        if end > reach:
            total += end - max(start, reach + 1) + 1
            reach = end
    return total


def piece_coverage(m8: Path, pieces: Path, min_identity: float = 0.0) -> dict[str, float]:
    """Fraction of each piece covered by translated matches of its genome to RefSeq viruses.

    min_identity keeps only matches with at least that amino acid identity.
    """
    hits: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with m8.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if float(row[2]) >= min_identity:
                start, end = sorted((int(row[4]), int(row[5])))
                hits[row[0]].append((start, end))
    found: dict[str, float] = {}
    with pieces.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            first, last = int(row["start"]), int(row["end"])
            inside = [
                (max(start, first), min(end, last))
                for start, end in hits.get(row["source"], ())
                if start <= last and end >= first
            ]
            found[row["name"]] = union_length(inside) / (last - first + 1)
    return found


def add_pieces(builder: Builder, fasta: Path, run_id: str) -> None:
    load_fasta(builder.con, fasta)
    builder.add(
        "sequence", "SELECT sequence_id, 'region' AS kind, length, bases, gc_content FROM _fasta"
    )
    builder.add(
        "sequence_source", f"SELECT sequence_id, {literal(run_id)} AS run_id, name FROM _fasta"
    )


def add_tools(
    builder: Builder, out: Path, fasta: Path, label: str, versions: dict[str, str]
) -> None:
    """Add the classifiers that finished under one output folder, going by its timing table."""
    timing = out / "timing.tsv"
    done = (
        {row[0] for row in csv.reader(timing.open(), delimiter="\t") if row[2] == "0"}
        if timing.exists()
        else set[str]()
    )
    if "genomad" in done:
        add_genomad(builder, out / "genomad", fasta, run_id=f"genomad.{label}", tool_version=versions.get("genomad"), database_version="1.9")  # fmt: skip
    if "virsorter2" in done:
        add_virsorter2(builder, out / "virsorter2", fasta, run_id=f"virsorter2.{label}", tool_version=versions.get("virsorter2"))  # fmt: skip
    if "cenotetaker3" in done:
        add_cenotetaker3(builder, next((out / "cenotetaker3").glob("*/*_virus_summary.tsv")).parent, fasta, run_id=f"cenotetaker3.{label}", tool_version=versions.get("cenotetaker3"), parameters={"prune_prophage": "F", "minimum_hallmark_genes_linear": 1})  # fmt: skip
    if "phabox2" in done:
        add_phabox2(builder, out / "phabox2", fasta, run_id=f"phabox2.{label}", tool_version=versions.get("phabox2"), parameters={"task": "phamer", "database": "phabox_db_v2_2"})  # fmt: skip
    tables = [p for p in (out / "jaeger").glob("*/*.tsv") if not p.name.endswith("_phages.tsv")]
    if tables and "jaeger" in done:
        add_jaeger(builder, tables[0], fasta, run_id=f"jaeger.{label}", tool_version=versions.get("jaeger"), parameters={"model": tables[0].parent.name})  # fmt: skip


def main() -> None:
    work, bundle = Path(sys.argv[1]), Path(sys.argv[2])
    for old in [*bundle.glob("*.parquet"), bundle / "manifest.json"]:
        old.unlink(missing_ok=True)
    bundle.mkdir(parents=True, exist_ok=True)
    table = work / "benchmark_sample.tsv"
    table.write_text("sample_id\tnote\nbenchmark\tPieces cut from sequences of known origin\n")
    init_bundle(bundle, table, bundle_id="classifier-benchmark", title="Viral classifiers on pieces of known origin")  # fmt: skip
    versions = {"genomad": "1.13.0", "virsorter2": "2.2.4", "cenotetaker3": "3.4.4", "phabox2": "2.2", "jaeger": "1.26.5"}  # fmt: skip

    truth: list[tuple[object, ...]] = []
    with (work / "bench" / "labels.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            viral = row["origin"] == "viral"
            truth.append((row["name"], "simulated_refseq" if viral else "non_viral", viral, row["origin"], row["source"], row["has_family"] or None, None, None, None))  # fmt: skip
    if (work / "recent_pieces" / "pieces.tsv").exists():
        m8 = work / "novelty" / "recent_vs_refseq.m8"
        pieces = work / "recent_pieces" / "pieces.tsv"
        coverage = piece_coverage(m8, pieces) if m8.exists() else {}
        close = piece_coverage(m8, pieces, 0.5) if m8.exists() else {}
        kinds: dict[str, tuple[str, str]] = {}
        with (work / "recent" / "genomes.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                study = "".join("N" if c.isdigit() else c for c in row["title"])
                kinds[row["accession"]] = (row["set"], study)
        with (work / "recent_pieces" / "pieces.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                group, study = kinds[row["source"]]
                truth.append((row["name"], "recent_cultured" if group == "cultured_phage" else "recent_uncultivated", True, "viral", row["source"], None, coverage.get(row["name"]), close.get(row["name"]), study))  # fmt: skip

    with Builder(bundle) as b:
        b.add_run("pieces", tool="cut_pieces", step="other", description="Cut 5 kb pieces from sequences of known origin.", parameters={"piece_length": 5000})  # fmt: skip
        add_pieces(b, work / "bench" / "pieces.fasta", "pieces")
        add_tools(b, work / "bench" / "out", work / "bench" / "pieces.fasta", "known", versions)
        if (work / "recent_pieces" / "pieces.fasta").exists():
            add_pieces(b, work / "recent_pieces" / "pieces.fasta", "pieces")
            add_tools(b, work / "recent_pieces" / "out", work / "recent_pieces" / "pieces.fasta", "recent", versions)  # fmt: skip
        b.con.execute(
            "CREATE OR REPLACE TEMP TABLE _truth (name VARCHAR, test_set VARCHAR, is_viral BOOLEAN, origin VARCHAR, "
            "source VARCHAR, has_family BOOLEAN, refseq_protein_coverage DOUBLE, refseq_close_coverage DOUBLE, study VARCHAR)"
        )
        b.con.executemany(
            "INSERT INTO _truth VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", [list(t) for t in truth]
        )
        b.add(
            "x_truth",
            "SELECT s.sequence_id, 'pieces' AS run_id, t.* EXCLUDE (name) FROM _truth t "
            "JOIN sequence_source s ON s.name = t.name AND s.run_id = 'pieces'",
        )
        b.declare(
            "x_truth",
            shape="label",
            one_row_is="what is known about the origin of one test piece",
            description="Known origin of each piece, for scoring classifiers.",
            key=["sequence_id", "run_id"],
        )
        b.describe("x_truth", TRUTH_COLUMNS)


if __name__ == "__main__":
    main()
