"""Command line: build a bundle from tool output, add to it, check it, or query it."""

import argparse
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

from capsomer.builder import Builder
from capsomer.bundle import run_sql
from capsomer.importers.checkm2 import add_checkm2
from capsomer.importers.checkv import add_checkv
from capsomer.importers.clusters import LAYOUTS, add_clusters
from capsomer.importers.contigs import add_contigs
from capsomer.importers.genomad import add_genomad
from capsomer.importers.hecatomb import import_hecatomb
from capsomer.importers.mapping import add_coverm, add_paf, add_read_mapping
from capsomer.importers.study import add_fastp, init_bundle
from capsomer.importers.viroforge import import_viroforge
from capsomer.importers.virsorter2 import add_virsorter2
from capsomer.validate import conforms, format_report, validate


def _pairs(items: list[str] | None) -> dict[str, str]:
    """KEY=VALUE arguments as a dictionary."""
    found: dict[str, str] = {}
    for item in items or []:
        key, _, value = item.partition("=")
        if not value:
            raise SystemExit(f"expected KEY=VALUE, got {item}")
        found[key] = value
    return found


def _add_parser(commands: Any) -> None:
    add = commands.add_parser("add", help="add one tool's results to a bundle")
    tools = add.add_subparsers(dest="tool", required=True)

    def tool(name: str, help_text: str, *, version: bool = True) -> argparse.ArgumentParser:
        parser = tools.add_parser(name, help=help_text)
        parser.add_argument("bundle", type=Path)
        parser.add_argument("--run-id", required=True, help="identifier for this run in the bundle")
        parser.add_argument("--parent", help="the pipeline run this run is a step of")
        parser.add_argument(
            "--input", action="append", help="run whose results this run used (repeat)"
        )
        if version:
            parser.add_argument("--tool-version")
        return parser

    def settings(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--parameter", action="append", metavar="NAME=VALUE", help="a setting of the run"
        )

    fastp = tool("fastp", "read totals from a fastp JSON report", version=False)
    fastp.add_argument("report", type=Path)
    fastp.add_argument("--library", required=True)

    contigs = tool("contigs", "contigs from an assembly FASTA file")
    contigs.add_argument("fasta", type=Path)
    contigs.add_argument("--assembler", required=True, help="name of the assembler")
    contigs.add_argument("--library", action="append", help="library assembled (repeat)")
    settings(contigs)

    genomad = tool("genomad", "geNomad viral calls and taxonomy")
    genomad.add_argument("results", type=Path, help="geNomad output folder")
    genomad.add_argument("fasta", type=Path, help="the FASTA file geNomad was given")
    genomad.add_argument("--database-version")
    genomad.add_argument("--preferred", action="store_true", help="mark as the pipeline's answer")

    checkv = tool("checkv", "CheckV completeness and contamination")
    checkv.add_argument("results", type=Path, help="CheckV output folder")
    checkv.add_argument("fasta", type=Path, help="the FASTA file CheckV was given")
    checkv.add_argument("--database-version")

    virsorter2 = tool("virsorter2", "VirSorter2 viral calls")
    virsorter2.add_argument("results", type=Path, help="VirSorter2 output folder")
    virsorter2.add_argument("fasta", type=Path, help="the FASTA file VirSorter2 was given")
    virsorter2.add_argument("--database-version")
    virsorter2.add_argument(
        "--preferred", action="store_true", help="mark as the pipeline's answer"
    )

    checkm2 = tool("checkm2", "CheckM2 completeness and contamination of genomes")
    checkm2.add_argument("report", type=Path, help="CheckM2 quality_report.tsv")
    checkm2.add_argument(
        "genomes", type=Path, help="folder of genome FASTA files CheckM2 was given"
    )
    checkm2.add_argument("--extension", default="fna", help="extension of the genome files")
    checkm2.add_argument("--binning-run", help="the run that formed the genomes")
    checkm2.add_argument("--database-version")

    clusters = tool("clusters", "clusters of sequences from a two-column table")
    clusters.add_argument("table", type=Path)
    clusters.add_argument("fasta", type=Path, help="the FASTA file that was clustered")
    clusters.add_argument("--clusterer", required=True, help="name of the clustering tool")
    clusters.add_argument("--level", required=True, help="species, genus, or any level you define")
    clusters.add_argument("--layout", choices=sorted(LAYOUTS), default="vclust")
    settings(clusters)

    coverm = tool("coverm", "coverage from a CoverM contig table")
    coverm.add_argument("table", type=Path)
    coverm.add_argument("fasta", type=Path, help="the FASTA file the reads were mapped to")
    coverm.add_argument(
        "--library", action="append", required=True, metavar="SAMPLE=LIBRARY",
        help="sample name in the CoverM header and the library it is (repeat)",
    )  # fmt: skip
    settings(coverm)

    mapping = tool("read-mapping", "a read mapping run, with per-read alignments if reads are held")
    mapping.add_argument("fasta", type=Path, help="the FASTA file the reads were mapped to")
    mapping.add_argument("--mapper", required=True, help="name of the read mapper")
    mapping.add_argument("--library", required=True)
    mapping.add_argument("--alignments", type=Path, help="SAM or BAM file")
    settings(mapping)

    paf = tool("paf", "alignments of sequences to sequences from a PAF file")
    paf.add_argument("paf", type=Path)
    paf.add_argument("query_fasta", type=Path)
    paf.add_argument("target_fasta", type=Path)
    paf.add_argument("--aligner", required=True, help="name of the aligner")
    settings(paf)


def _add(args: argparse.Namespace) -> None:
    run_id: str = args.run_id
    parent: str | None = args.parent
    inputs: list[str] = args.input or []
    with Builder(args.bundle) as builder:
        actions: dict[str, Callable[[], object]] = {
            "fastp": lambda: add_fastp(
                builder, args.report, run_id=run_id, library_id=args.library, parent=parent
            ),
            "contigs": lambda: add_contigs(
                builder, args.fasta, tool=args.assembler, tool_version=args.tool_version,
                libraries=args.library or [], parameters=_pairs(args.parameter), run_id=run_id, parent=parent, inputs=inputs,
            ),
            "genomad": lambda: add_genomad(
                builder, args.results, args.fasta, tool_version=args.tool_version,
                database_version=args.database_version, preferred=args.preferred, run_id=run_id, parent=parent, inputs=inputs,
            ),
            "checkv": lambda: add_checkv(
                builder, args.results, args.fasta, tool_version=args.tool_version,
                database_version=args.database_version, run_id=run_id, parent=parent, inputs=inputs,
            ),
            "virsorter2": lambda: add_virsorter2(
                builder, args.results, args.fasta, tool_version=args.tool_version,
                database_version=args.database_version, preferred=args.preferred,
                run_id=run_id, parent=parent, inputs=inputs,
            ),
            "checkm2": lambda: add_checkm2(
                builder, args.report, args.genomes, extension=args.extension,
                binning_run=args.binning_run, tool_version=args.tool_version,
                database_version=args.database_version, run_id=run_id, parent=parent, inputs=inputs,
            ),
            "clusters": lambda: add_clusters(
                builder, args.table, args.fasta, tool=args.clusterer, level=args.level,
                layout=args.layout, tool_version=args.tool_version,
                parameters=_pairs(args.parameter), run_id=run_id, parent=parent, inputs=inputs,
            ),
            "coverm": lambda: add_coverm(
                builder, args.table, args.fasta, libraries=_pairs(args.library),
                tool_version=args.tool_version, parameters=_pairs(args.parameter), run_id=run_id, parent=parent, inputs=inputs,
            ),
            "read-mapping": lambda: add_read_mapping(
                builder, args.fasta, tool=args.mapper, library_id=args.library,
                alignments=args.alignments, tool_version=args.tool_version,
                parameters=_pairs(args.parameter), run_id=run_id, parent=parent, inputs=inputs,
            ),
            "paf": lambda: add_paf(
                builder, args.paf, args.query_fasta, args.target_fasta, tool=args.aligner,
                tool_version=args.tool_version, parameters=_pairs(args.parameter), run_id=run_id, parent=parent, inputs=inputs,
            ),
        }  # fmt: skip
        result = actions[args.tool]()
    print(
        f"Added {args.tool} run {args.run_id} to {args.bundle}" + (f": {result}" if result else "")
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="capsomer", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    check = commands.add_parser("validate", help="check a bundle against the standard")
    check.add_argument("bundle", type=Path)

    init = commands.add_parser("init", help="start a bundle from a study's sample table")
    init.add_argument("bundle", type=Path, help="folder to write the bundle into")
    init.add_argument("sample_table", type=Path, help="tab-separated, one row per library")
    init.add_argument("--bundle-id")
    init.add_argument("--title")

    _add_parser(commands)

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
    if args.command == "init":
        libraries = init_bundle(
            args.bundle, args.sample_table, bundle_id=args.bundle_id, title=args.title
        )
        print(f"Started {args.bundle} with {libraries} libraries. Add results with: capsomer add")
        return 0
    if args.command == "add":
        _add(args)
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
