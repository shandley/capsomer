"""Cut fixed-length pieces from every sequence of a FASTA file, evenly spaced."""

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("out", type=Path, help="folder for pieces.fasta and pieces.tsv")
    parser.add_argument("--length", type=int, default=5000)
    parser.add_argument("--per-sequence", type=int, default=3)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    records: dict[str, str] = {}
    name = ""
    for line in args.fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        else:
            records[name] += line.strip().upper()
    with (
        (args.out / "pieces.fasta").open("w") as pieces,
        (args.out / "pieces.tsv").open("w") as table,
    ):
        table.write("name\tsource\tstart\tend\n")
        for source, bases in records.items():
            count = min(args.per_sequence, len(bases) // args.length)
            if count == 0:
                continue
            gap = (len(bases) - count * args.length) // (count + 1)
            for i in range(count):
                start = gap * (i + 1) + args.length * i
                piece = bases[start : start + args.length]
                if set(piece) - set("ACGT"):
                    continue
                name = f"{source}__{start + 1}_{start + args.length}"
                pieces.write(f">{name}\n{piece}\n")
                table.write(f"{name}\t{source}\t{start + 1}\t{start + args.length}\n")


if __name__ == "__main__":
    main()
