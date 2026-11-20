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
- `src/capsomer/importers/`: importers. One reads Hecatomb version 1 results. One reads the ground truth of a [ViroForge](https://github.com/shandley/viroforge) synthetic dataset, which gives a dataset that no analysis pipeline shaped.
- `examples/questions.sql`: test questions from the specification, written as SQL.

The draft rests on three background surveys of existing data formats, virome tool outputs and consortium data models. Those are working notes and are not published here.

## Try it

The library needs Python 3.12 and [uv](https://docs.astral.sh/uv/).

```sh
uv sync
uv run pytest

# Build a dataset from Hecatomb results, then check it.
uv run capsomer import-hecatomb path/to/hecatomb/results my_dataset --n-samples 4
uv run capsomer validate my_dataset

# Build a ground-truth dataset from a ViroForge simulation.
uv run capsomer import-viroforge path/to/viroforge/output truth_dataset

# Ask questions of it.
uv run capsomer sql my_dataset --file examples/questions.sql
uv run capsomer sql my_dataset "SELECT kind, count(*) FROM sequence GROUP BY 1"
```

A dataset is plain Parquet, so DuckDB, Python, R or any other Parquet reader can open it without this library.

## Background

Capsomer is not tied to any pipeline. A dataset can come from a read classifier that reports only taxa, from an assembly workflow, or from a pipeline that does both. The first importer is for [Hecatomb](https://github.com/shandley/hecatomb), because its results were at hand; importers for other tools will follow.

## Licence

MIT. See `LICENSE`.
