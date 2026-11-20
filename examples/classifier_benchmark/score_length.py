"""Score geNomad, Jaeger and PhaBOX2 on pieces of another length.

Usage: score_length.py <work folder> <length>
Reads <work>/len<length>/, as laid out by make_benchmark_set.py and cut_pieces.py, and the
similarity search of the recent genomes against RefSeq viruses in <work>/novelty/.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

from build_benchmark import piece_coverage


def called(folder: Path) -> dict[str, set[str]]:
    """Piece names each tool called viral."""
    found: dict[str, set[str]] = {}
    summary = next((folder / "genomad").glob("*_summary/*_virus_summary.tsv"), None)
    if summary:
        with summary.open() as handle:
            found["genomad"] = {
                r["seq_name"].split("|")[0] for r in csv.DictReader(handle, delimiter="\t")
            }
    tables = [p for p in (folder / "jaeger").glob("*/*.tsv") if not p.name.endswith("_phages.tsv")]
    if tables:
        with tables[0].open() as handle:
            found["jaeger"] = {
                r["contig_id"]
                for r in csv.DictReader(handle, delimiter="\t")
                if r["prediction"] in ("phage", "virus")
            }
    phamer = folder / "phabox2" / "final_prediction" / "phamer_prediction.tsv"
    if phamer.exists():
        with phamer.open() as handle:
            found["phabox2"] = {
                r["Accession"]
                for r in csv.DictReader(handle, delimiter="\t")
                if r["Pred"] == "virus"
            }
    return found


def main() -> None:
    work, length = Path(sys.argv[1]), sys.argv[2]
    base = work / f"len{length}"
    group: dict[str, str] = {}
    with (base / "known" / "labels.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            group[row["name"]] = "known viruses" if row["origin"] == "viral" else "not viral"
    m8 = work / "novelty" / "recent_vs_refseq.m8"
    pieces = base / "recent" / "pieces.tsv"
    any_match, close = piece_coverage(m8, pieces), piece_coverage(m8, pieces, 0.5)
    for name in any_match:
        group[name] = (
            "recent, close to RefSeq" if close[name] >= 0.5
            else "recent, distant" if any_match[name] >= 0.5
            else "recent, very distant"
        )  # fmt: skip
    calls = called(base / "out")
    if {"genomad", "jaeger"} <= set(calls):
        calls["genomad_or_jaeger"] = calls["genomad"] | calls["jaeger"]
        calls["genomad_and_jaeger"] = calls["genomad"] & calls["jaeger"]
    totals: dict[str, int] = defaultdict(int)
    hits: dict[tuple[str, str], int] = defaultdict(int)
    for name, kind in group.items():
        totals[kind] += 1
        for tool, names in calls.items():
            hits[(kind, tool)] += name in names
    tools = sorted(calls)
    print(f"pieces of {length} bases\n" + "\t".join(["group", "pieces", *tools]))
    for kind in sorted(totals):
        print(
            "\t".join(
                [kind, str(totals[kind]), *(f"{hits[(kind, t)] / totals[kind]:.3f}" for t in tools)]
            )
        )


if __name__ == "__main__":
    main()
