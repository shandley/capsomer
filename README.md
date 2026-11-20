# Capsomer

Capsomer is a proposed data standard for the results of a virome analysis. It defines what goes inside the result files: the sequences, the counts, the taxonomy and other labels, and the links among reads, contigs, genomes and clusters.

Status: a design draft under discussion. There is no schema file and no software yet. Comments are welcome as GitHub issues.

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

- `spec/data-model-draft.md`: the current draft of the standard, with a log of decisions made and decisions still open.

The draft rests on three background surveys of existing data formats, virome tool outputs and consortium data models. Those are working notes and are not published here.

## Background

The design began from the data structure in [Hecatomb](https://github.com/shandley/hecatomb), a virome pipeline that analyses reads and contigs together.

## Licence

Not yet chosen.
