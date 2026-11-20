"""Download a random sample of viral genomes deposited at NCBI on or after a date.

A classifier cannot have been trained on a sequence that was not public when
its database was built. Genomes deposited after the newest database date among
the tools being compared are therefore unseen by all of them, although close
relatives of some will have been public earlier.

Two sets are drawn. Cultured phages were found by growing them, so no
classifier chose them. Uncultivated genomes are the more novel, but each was
picked out of a metagenome by some classifier, which favours that classifier
here. One deposit can also hold thousands of them, so the sample takes a
limited number from any one study, recognised by its naming pattern.
"""

import argparse
import json
import random
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UNCULTIVATED = "(MAG[Title] OR uncultured[Title] OR metagenome[Title])"
SETS = {
    "cultured_phage": f"Viruses[Organism] AND (phage[Title] OR bacteriophage[Title]) NOT {UNCULTIVATED}",
    "uncultivated": f"Viruses[Organism] AND {UNCULTIVATED}",
}


def call(endpoint: str, **params: str) -> str:
    query = urllib.parse.urlencode({**params, "tool": "virome_benchmark"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(f"{EUTILS}{endpoint}?{query}", timeout=120) as reply:
                time.sleep(0.4)
                return reply.read().decode()
        except OSError:
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"NCBI did not answer {endpoint}")


def identifiers(term: str) -> list[str]:
    reply = call("esearch.fcgi", db="nuccore", term=term, retmax="100000", idtype="acc")
    return [part.split("</Id>")[0] for part in reply.split("<Id>")[1:]]


def titles(accessions: list[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    for start in range(0, len(accessions), 200):
        batch = accessions[start : start + 200]
        reply = json.loads(call("esummary.fcgi", db="nuccore", id=",".join(batch), retmode="json"))
        for entry in reply.get("result", {}).values():
            if isinstance(entry, dict) and "accessionversion" in entry:
                found[str(entry["accessionversion"])] = str(entry.get("title", ""))
    return found


def study_key(title: str) -> str:
    """A stand-in for the depositing study: the isolate name with digits removed.

    Uncultivated genomes from one deposit share a naming scheme (KCR_vOTU_03142,
    KCR_vOTU_00518), whatever taxon each one was given. Titles with no isolate
    name fall back to the whole title with digits removed.
    """
    named = re.search(r"isolate (\S+)", title)
    return re.sub(r"\d+", "N", named.group(1) if named else title)


def capped(
    accessions: list[str], title: dict[str, str], limit: int, rng: random.Random
) -> list[str]:
    """A random order of the accessions, keeping at most `limit` from one study."""
    kept: list[str] = []
    seen: dict[str, int] = {}
    for accession in rng.sample(accessions, len(accessions)):
        study = study_key(title.get(accession, accession))
        if seen.get(study, 0) < limit:
            seen[study] = seen.get(study, 0) + 1
            kept.append(accession)
    return kept


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--since", required=True, help="first deposit date, as YYYY/MM/DD")
    parser.add_argument("--until", default="3000", help="last deposit date, as YYYY/MM/DD")
    parser.add_argument("--sets", nargs="+", default=sorted(SETS), choices=sorted(SETS))
    parser.add_argument("--per-set", type=int, default=300)
    parser.add_argument("--min-length", type=int, default=5000)
    parser.add_argument("--max-length", type=int, default=500000)
    parser.add_argument(
        "--per-study", type=int, default=40, help="most genomes from one naming pattern"
    )
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    with (
        (args.out / "genomes.fasta").open("w") as fasta,
        (args.out / "genomes.tsv").open("w") as table,
    ):
        table.write("accession\tset\tstudy\ttitle\n")
        for name in args.sets:
            organisms = SETS[name]
            term = (
                f"({organisms}) AND {args.min_length}:{args.max_length}[SLEN] "
                f'AND ("{args.since}"[PDAT] : "{args.until}"[PDAT]) NOT srcdb_refseq[PROP]'
            )
            found = sorted(identifiers(term))
            chosen = capped(found, titles(found), args.per_study, rng)[: args.per_set]
            print(f"{name}: {len(found)} records since {args.since}, {len(chosen)} sampled")
            for start in range(0, len(chosen), 50):
                batch = chosen[start : start + 50]
                text = call(
                    "efetch.fcgi", db="nuccore", id=",".join(batch), rettype="fasta", retmode="text"
                )
                for record in text.split(">")[1:]:
                    header, _, body = record.partition("\n")
                    accession, _, title = header.partition(" ")
                    table.write(f"{accession}\t{name}\t{study_key(title)}\t{title}\n")
                    fasta.write(f">{accession}\n{body.replace(chr(10), '').upper()}\n")


if __name__ == "__main__":
    main()
