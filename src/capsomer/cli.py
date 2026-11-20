"""Command line: build a bundle from pipeline output, or check a bundle."""

import argparse
import sys
from pathlib import Path

from capsomer.bundle import run_sql
from capsomer.importers.hecatomb import import_hecatomb
from capsomer.importers.viroforge import import_viroforge
from capsomer.validate import conforms, format_report, validate


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="capsomer", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    check = commands.add_parser("validate", help="check a bundle against the standard")
    check.add_argument("bundle", type=Path)

    hecatomb = commands.add_parser("import-hecatomb", help="build a bundle from Hecatomb results")
    hecatomb.add_argument(
        "results", type=Path, help="folder holding bigtable.tsv and seqtable.fasta"
    )
    hecatomb.add_argument("bundle", type=Path, help="folder to write the bundle into")
    hecatomb.add_argument("--samples", help="comma-separated sample names to import")
    hecatomb.add_argument("--n-samples", type=int, help="import the first N samples by name")
    hecatomb.add_argument("--sample-table", type=Path, help="tab-separated study sample table")
    hecatomb.add_argument("--sample-id-column", default="sampleID")
    hecatomb.add_argument("--title")

    viroforge = commands.add_parser(
        "import-viroforge", help="build a ground-truth bundle from a ViroForge dataset"
    )
    viroforge.add_argument("dataset", type=Path, help="ViroForge output folder")
    viroforge.add_argument("bundle", type=Path, help="folder to write the bundle into")
    viroforge.add_argument("--sample-id")

    sql = commands.add_parser("sql", help="run SQL against a bundle")
    sql.add_argument("bundle", type=Path)
    sql.add_argument("statement", nargs="?", help="SQL to run")
    sql.add_argument("--file", type=Path, help="file of SQL statements to run instead")

    args = parser.parse_args(argv)
    if args.command == "sql":
        print(run_sql(args.bundle, args.file.read_text() if args.file else args.statement))
        return 0
    if args.command == "import-viroforge":
        truth = import_viroforge(args.dataset, args.bundle, sample_id=args.sample_id)
        print(f"Imported ViroForge sample {truth.sample} into {args.bundle}")
        for name, rows in truth.rows.items():
            print(f"  {name:18} {rows:>12,} rows")
        print(
            f"  reads whose source genome is not in the reference file: {truth.reads_without_source:,}"
        )
    if args.command == "import-hecatomb":
        summary = import_hecatomb(
            args.results,
            args.bundle,
            samples=args.samples.split(",") if args.samples else None,
            n_samples=args.n_samples,
            sample_table=args.sample_table,
            sample_id_column=args.sample_id_column,
            title=args.title,
        )
        print(f"Imported {len(summary.samples)} samples into {args.bundle}")
        for name, rows in summary.rows.items():
            print(f"  {name:18} {rows:>12,} rows")
        print(
            f"  identical read clusters merged within a library: {summary.identical_clusters_merged:,}"
        )
    findings = validate(args.bundle)
    print(format_report(findings))
    return 0 if conforms(findings) else 1


if __name__ == "__main__":
    sys.exit(main())
