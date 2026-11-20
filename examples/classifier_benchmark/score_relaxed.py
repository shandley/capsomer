"""Compare geNomad's default and relaxed filters with Jaeger on pieces of one length.

Usage: score_relaxed.py <work folder> <length>
Needs <work>/len<length>/out_relaxed/genomad, a geNomad run with --relaxed on the same pieces.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

from build_benchmark import piece_coverage
from score_length import called


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
    default = called(base / "out")
    relaxed = called(base / "out_relaxed")["genomad"]
    calls = {
        "genomad_default": default["genomad"],
        "genomad_relaxed": relaxed,
        "jaeger": default["jaeger"],
        "genomad_default_or_jaeger": default["genomad"] | default["jaeger"],
        "genomad_relaxed_or_jaeger": relaxed | default["jaeger"],
    }
    totals: dict[str, int] = defaultdict(int)
    hits: dict[tuple[str, str], int] = defaultdict(int)
    for name, kind in group.items():
        totals[kind] += 1
        for tool, names in calls.items():
            hits[(kind, tool)] += name in names
    print(f"pieces of {length} bases\n" + "\t".join(["group", "pieces", *calls]))
    for kind in sorted(totals):
        print(
            "\t".join(
                [kind, str(totals[kind]), *(f"{hits[(kind, t)] / totals[kind]:.3f}" for t in calls)]
            )
        )


if __name__ == "__main__":
    main()
