# Capsomer

Capsomer is a proposed data standard for the results of a virome analysis. It defines what goes inside the result files: the sequences, the counts, the taxonomy and other labels, and the links among reads, contigs, genomes and clusters.

Status: an early draft. The specification, a schema file and a first reference implementation exist, and all three will change. Comments are welcome as GitHub issues.

## Why

Consortium metadata schemas record that a result file exists and which tool made it. None says what the file contains, so every virome pipeline writes its own tables with its own columns. Results from two pipelines, or two studies, cannot be combined without custom code.

## The design in brief

Every table in a Capsomer dataset is one of four shapes.

| Shape | One row is | Familiar example |
|---|---|---|
| Thing | one item with an identifier | sample data; the sequences behind a count table |
| Count | one thing in one library, with numbers | an OTU table |
| Label | one statement about one thing, made by one tool run | a taxonomy table |
| Link | one statement about two things, made by one tool run | a sequence identity table |

A catalog of named tables sits on top of the shapes, so that a taxonomy table or a read count table means the same thing in every dataset. A dataset is a folder of Parquet tables plus one manifest file that describes them.

The name comes from the capsomer, the repeating protein unit from which a virus capsid is built.

## What the design is meant to do

- Hold read-based and assembly-based results in the same tables, so reads that never assemble stay visible next to the contigs.
- Keep sequences that match nothing in any database, and record when a search found nothing.
- Store calls from several tools side by side when they disagree.
- Reference taxa by identifier within a dated release, so a dataset can be restated under a newer taxonomy.
- Record analytical choices such as clustering thresholds without prescribing them.
- Be readable by a person or an AI agent from the manifest alone.

## What is in this repository

- `spec/data-model-draft.md`: the current draft of the standard, with a log of decisions made.
- `src/capsomer/schema/capsomer.yaml`: the catalog of tables, written in LinkML (the schema language that the MIxS sample vocabulary uses). Each table states its shape and what one row is.
- `src/capsomer/`: a Python library that loads the schema, computes sequence identifiers, writes tables, builds the manifest and checks a dataset against the standard.
- `src/capsomer/importers/`: importers. Some add one tool's results to a dataset: contigs from any assembler, geNomad, CheckV, clusters from vclust or MMseqs2, CoverM coverage, read mappings (SAM or BAM), sequence alignments (PAF) and fastp read totals. Two build a whole dataset: one from Hecatomb version 1 results, and one from the ground truth of a [ViroForge](https://github.com/shandley/viroforge) synthetic dataset.
- `examples/questions.sql`: test questions from the specification, written as SQL.
- `examples/viroforge_skin/`: a worked example. Two assemblers, geNomad, CheckV, two clustering tools and a read mapper were run on simulated reads, and their results were added to the simulation's ground truth, so every result can be checked against what is true.

The draft rests on three background surveys of existing data formats, virome tool outputs and consortium data models. Those are working notes and are not published here.

## Try it

The library needs Python 3.12 and [uv](https://docs.astral.sh/uv/).

```sh
uv sync
uv run pytest

# Start a dataset from a sample table, then add each tool's results.
uv run capsomer init my_dataset samples.tsv
uv run capsomer add contigs my_dataset contigs.fasta --run-id assembly --assembler megahit --library sample1
uv run capsomer add genomad my_dataset genomad_out/ contigs.fasta --run-id genomad --input assembly
uv run capsomer add checkv my_dataset checkv_out/ viruses.fna --run-id checkv --input genomad
uv run capsomer validate my_dataset

# Or build a dataset from Hecatomb results.
uv run capsomer import-hecatomb path/to/hecatomb/results my_dataset --n-samples 4

# Build a ground-truth dataset from a ViroForge simulation.
uv run capsomer import-viroforge path/to/viroforge/output truth_dataset

# Ask questions of it.
uv run capsomer sql my_dataset --file examples/questions.sql
uv run capsomer sql my_dataset "SELECT kind, count(*) FROM sequence GROUP BY 1"
```

A dataset is plain Parquet, so DuckDB, Python, R or any other Parquet reader can open it without this library.

## Background

Capsomer is not tied to any pipeline. A dataset can come from a read classifier that reports only taxa, from an assembly workflow, or from a pipeline that does both. Tools report results by sequence name; an importer reads the FASTA file the tool was given, computes each sequence's identifier from its bases, and joins on that. `capsomer add --help` lists the tools that can be added today.

## Licence

MIT. See `LICENSE`.
