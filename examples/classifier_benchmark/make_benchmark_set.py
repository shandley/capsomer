"""Cut a labelled test set for viral classifiers from a ViroForge dataset.

Every sequence ViroForge simulated from has a known origin. This script cuts
fixed-length pieces from those sequences, so a classifier can be scored on
pieces the size of typical contigs with no assembly step in between.

Left out: the phiX spike-in (a real phage that is a contaminant here),
ribosomal RNA (too short) and the reagent bacteria, which ViroForge builds
from random sequence when it has no genome database for them.
"""

import argparse
import json
from pathlib import Path

SKIP = {"phix", "rrna", "reagent_bacteria"}


def fasta(path: Path) -> dict[str, str]:
    found: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            found[name] = ""
        else:
            found[name] += line.strip().upper()
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path, help="ViroForge output folder")
    parser.add_argument("out", type=Path, help="folder for pieces.fasta and labels.tsv")
    parser.add_argument("--length", type=int, default=5000, help="length of each piece")
    parser.add_argument("--per-virus", type=int, default=3, help="pieces per viral genome")
    parser.add_argument("--per-other", type=int, default=2, help="pieces per non-viral sequence")
    parser.add_argument(
        "--references",
        type=Path,
        help="ViroForge's folder of reference fragments (bacterial, fungal, archaeal, host). "
        "A dataset uses only some of them; this adds every one as a non-viral sequence.",
    )
    args = parser.parse_args()

    meta = json.loads(next((args.dataset / "metadata").glob("*_metadata.json")).read_text())
    bench = meta["benchmarking"]
    taxonomy = bench.get("taxonomy", {})
    origin = {
        c["genome_id"]: c.get("contaminant_type", "contaminant")
        for c in bench.get("contamination_manifest", {}).get("contaminants", [])
    }
    for s in meta["sequences"]:
        if s.get("sequence_type") == "viral":
            origin[s["genome_id"]] = "viral"
    genomes = fasta(next((args.dataset / "fasta").glob("*.fasta")))
    genomes = {name: bases for name, bases in genomes.items() if origin.get(name) == "viral"}
    if args.references:
        for kind in ("bacterial", "fungal", "archaeal", "host"):
            for name, bases in fasta(args.references / f"{kind}_fragments.fasta").items():
                genomes[name] = bases
                origin[name] = "host_dna" if kind == "host" else f"{kind}_background"
    else:
        genomes = fasta(next((args.dataset / "fasta").glob("*.fasta")))

    args.out.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    with (
        (args.out / "pieces.fasta").open("w") as pieces,
        (args.out / "labels.tsv").open("w") as labels,
    ):
        labels.write("name\tsource\torigin\thas_family\tclass\tstart\tend\n")
        for source, bases in genomes.items():
            kind = origin.get(source)
            if kind is None or kind in SKIP or len(bases) < args.length or bases in seen:
                continue
            seen.add(bases)
            wanted = args.per_virus if kind == "viral" else args.per_other
            count = min(wanted, len(bases) // args.length)
            spare = len(bases) - count * args.length
            gap = spare // (count + 1)
            entry = taxonomy.get(source, {})
            for i in range(count):
                start = gap * (i + 1) + args.length * i
                piece = bases[start : start + args.length]
                if set(piece) - set("ACGT"):
                    continue
                name = f"{source}__{start + 1}_{start + args.length}"
                pieces.write(f">{name}\n{piece}\n")
                has_family = entry.get("is_known", "") if kind == "viral" else ""
                labels.write(
                    f"{name}\t{source}\t{kind}\t{has_family}\t{entry.get('class') or ''}\t"
                    f"{start + 1}\t{start + args.length}\n"
                )


if __name__ == "__main__":
    main()
