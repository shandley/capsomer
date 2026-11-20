# Capsomer: a data standard for virome results (draft 3)

Status: draft, last revised 2026-11-20. A schema file and a first reference implementation exist in this repository. Where they differ from this document, this document is the intent and the difference is a bug.

A capsomer is the repeating protein unit from which a virus capsid is built. The standard is named for it because every table in a dataset is one of four repeating shapes.

This draft rests on three background surveys: existing data formats (01), what virome tools and reference resources emit (02), and how other consortia modelled their data (03). They are working notes and are not part of this repository. Statements marked with a survey number below come from them.

## What the standard covers

Existing consortium schemas record that a result file exists and what produced it. They do not say what is inside the file. That holds for the Human Virome Program (HVP) metadata dictionary, the Human Microbiome Project schema and the National Microbiome Data Collaborative schema (survey 03).

This standard defines the inside: the sequences, counts, labels and links that a virome analysis produces.

A dataset saved in this format is called a bundle: a folder of tables plus one manifest file that describes them, for one study. A consortium store is built by combining bundles. A bundle carries its own sample and library tables, so it can be read with no access to a consortium database.

The standard is for any virome researcher, inside HVP or not. HVP has published no specification for the content of its result files (survey 03), so nothing here is fitted to one. Identifiers from other systems (an HVP sample identifier, a BioSample accession) are attached through the `identifier_map` table.

## The four shapes

Most biological datasets reduce to a few kinds of table: a table of samples, a count table, a taxonomy table, a sequence identity table. The standard is built from four such shapes, and every table in a bundle is one of them.

| Shape | One row is | Familiar example |
|---|---|---|
| Thing | one item with an identifier | sample data; the sequences behind a count table |
| Count | one thing in one library, with numbers | an OTU table |
| Label | one statement about one thing, made by one run | a taxonomy table |
| Link | one statement about two things, made by one run | a sequence identity table |

A phyloseq object, the R structure widely used for 16S data, is one thing table for samples, one count table, one label table and a tree. A virome analysis needs three additions, and all three fit the same shapes:

- Several levels of thing: read cluster, contig, genome, cluster of genomes, taxon. Links join the levels, so counts can be summed from one level to the next.
- Many labels for one thing. Tools disagree, so a label table holds one row per run.
- Links between things: a read cluster aligned to a contig, a contig placed in a genome.

### Columns each shape requires

**Thing.** An identifier column named `<table>_id`, unique within the table. The other columns hold properties of the item itself. A thing that a run created, such as a contig or a cluster, also carries `run_id`.

**Count.** The identifier of the thing counted, `library_id`, `run_id`, and one or more number columns. There is one row per thing, per library, per run. A thing with no row for a library was not observed there. The manifest declares each number column: its unit, whether it is raw or normalized, and whether it may be summed.

**Label.** The identifier of the thing described, `run_id`, and one or more value columns. A value can be a term, a number, or the identifier of a reference item such as a taxon. Two runs that disagree produce two rows. Five optional columns have the same meaning in every label table: `score`, `score_type`, `confidence`, `is_preferred` (the final answer of the pipeline that produced the row, at most one per thing per pipeline) and `detail` (tool-specific extras as a JSON string).

A label table also records the absence of an answer. A run that examined a thing and reached no answer writes a row with the value empty and the `outcome` column set to `none`. A thing with no row was not examined by that run. Without this, "searched and found nothing" cannot be told apart from "never searched".

**Link.** The identifiers of two things, `run_id`, and columns that describe the pair. When both things come from the same table, the column names carry a role: `query_sequence_id`, `target_sequence_id`.

A label answers "what is this thing?" with a value from a vocabulary or a reference, such as a taxon. A link relates two things the study observed or built, such as two sequences.

### Conventions for every table

- A column named `<table>_id`, or ending in it, refers to that thing table.
- Where a table can describe more than one kind of thing, it has one identifier column per kind, and exactly one is filled in each row.
- No cell holds a list.
- Any table may carry extra columns, provided the manifest describes them. Study metadata, such as a treatment group, goes in extra columns of `sample`.

## Rules

1. Every table is one of the four shapes and states what one row is.
2. An identifier never depends on a prediction. A sequence keeps its identifier whether or not a tool later calls it viral, circular or a provirus.
3. Names are never keys. Taxa are referenced by identifier within a named release. Samples are referenced by identifier, not by file name or column header.
4. Thresholds and database versions are data. Each run records the cutoffs it applied and the reference releases it read.
5. Read-level and assembly-level results use the same tables. A read cluster and a contig are both sequences, and both can carry labels and counts.
6. The standard describes an analysis and does not prescribe one. Clustering levels, identity cutoffs, quality filters and tool choices belong to the researcher. The standard requires that they are recorded, never that they take a particular value.
7. The container does not merge what the analysis kept apart. The one exception is stated under "Read clusters".
8. Every result can be followed to the things it came from. A read, the contig it maps to, the genome that contig belongs to and the taxon each was given are joined by identifiers, so information recorded at one level can be used at another.

## Identifiers

**Sequences.** A sequence identifier is the GA4GH sequence digest of the sequence as written, in upper case. GA4GH is the genomics standards body; its digest is a truncated SHA-512 hash with the prefix `SQ.`. The sequence `ACGT` has the identifier `SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2` (computed locally on 2026-11-19). The standard adds no rule of its own:

- A sequence and its reverse complement have different identifiers.
- Two rotations of one circular genome have different identifiers. Whether a sequence is circular is a prediction, and rule 2 forbids identifiers that depend on one.
- Both relations are recorded as rows in `sequence_relation`.

An optional column on `sequence` holds the digest of whichever strand sorts first. It helps software find reverse-complement pairs, and it is never used as a key.

**Read pairs.** When a pipeline keeps both reads of a pair as one unit, the row holds both reads in two columns, each as sequenced. The identifier is the same hash taken over the first read, a comma, and the second read, with the prefix `RP.` in place of `SQ.` so it cannot be mistaken for the digest of a single sequence. Written in full, GA4GH identifiers carry the namespace `ga4gh:` and read pair identifiers carry `capsomer:`, so the two can never collide. Two optional columns hold each read's plain digest, so a pair can be matched to the same read stored alone in another bundle.

**Genomes.** A genome identifier is derived from the identifiers of its member sequences, sorted and hashed, with the prefix `GN.`. The same set of contigs is then the same genome in any bundle.

**Taxa.** A taxon identifier combines the reference release with the identifier that release uses.

**Everything else.** Sample, library, run and cluster identifiers are local to the bundle. The manifest gives the bundle an identifier of its own, which keeps local identifiers apart when bundles are combined.

## Version 1 catalog

The shapes make any table readable. The catalog makes common tables mean the same thing in every bundle, by fixing each table's name and columns. Version 1 covers what a phyloseq object holds, plus the links between reads, contigs, genomes and clusters.

A valid bundle needs the tables marked required, at least one count table, and the manifest. That is eight tables.

### Things

| Table | One row is | Main columns |
|---|---|---|
| `sample` (required) | a sample | Control type is required. MIxS terms (the Genomic Standards Consortium's sample vocabulary) where they exist. |
| `library` (required) | a sequencing library | Its sample, nucleic acid type, enrichment, amplification, whether it is stranded, single or paired layout, and read and fragment totals before and after quality control. |
| `run` (required) | one execution of one tool | Tool, version, the step it performed, and the pipeline run it belongs to. |
| `reference_entry` | one database entry that something in the bundle matched | Its release, accession, description, length, molecule (nucleotide or protein) and taxon. |
| `reference_release` (required) | one dated version of a database, taxonomy or catalog | Name, version, date, source, checksum, and kind: formal taxonomy, catalog of uncultivated viruses, or sequence database. |
| `sequence` (required) | one nucleotide sequence, or one read pair kept as a unit | Kind (read cluster, read pair cluster, contig or reference), length, the bases, and the name its producer gave it. A read pair has a second column of bases. Optional columns hold GC content and a complexity score, which help separate low-complexity artifacts from unknown biology. |
| `taxon` (required) | one taxon in one release | The release's own identifier, name, rank, parent taxon, genome composition (Baltimore class) and host category. A group in a catalog of uncultivated viruses is a taxon too, under that catalog's release. |
| `genome` | one unit treated as a virus genome | Kind: single sequence, bin, or group of segments. |
| `cluster` | one cluster of genomes or of sequences | Level (strain, species, genus, family, or a level the researcher defines) and the representative. A cluster needs no taxon name, so unclassified groups are first-class. |
| `read` | one original read | Its library, the name the sequencer gave it, the read cluster it belongs to, and which end of the fragment it is. Optional, and large: one row per read. |

`taxon` may also carry one column per major rank, holding the name of the ancestor at that rank. With those columns it reads like the taxonomy table of a phyloseq object. They are for reading; joins use identifiers.

### Counts

| Table | One row is | Number columns |
|---|---|---|
| `read_count` | one read cluster in one library, in one read form | The count, the read form it was counted in, and for independent reads the end of the fragment they came from (see "Read clusters"). |
| `coverage` | one contig in one library | fragments mapped, mean depth, fraction of the contig covered |

### Labels

| Table | One row is | Value columns |
|---|---|---|
| `taxon_assignment` (required) | one taxon given to one sequence, genome or cluster, or the record that a search found none | The taxon, the search type and the rule (top hit, lowest common ancestor). Each step of a pipeline writes its own rows; the final answer is flagged as preferred. |
| `run_parameter` | one setting of one run | Name, value, unit. Every threshold that decides which rows exist is listed here. |
| `identifier_map` | one external identifier for one sample, library, sequence or genome | The authority that issued it and the identifier. |

### Links

| Table | One row is | Pair columns |
|---|---|---|
| `alignment` | one alignment of one sequence to another sequence or to a reference entry | The search type, coordinates on each side (1-based, inclusive, in that side's own units), reading frame for translated searches, identity, length, e-value, bit score, and a flag for the primary alignment. For a read pair, which of the two reads aligned. |
| `genome_member` | one sequence in one genome | Segment label where relevant. |
| `cluster_member` | one member of one cluster | The member is a genome or a sequence in the bundle, or a reference entry. A cluster with no reference member is novel with respect to that release. |
| `sequence_relation` | one stated relation between two sequences | Reverse complement, rotation, or part of a longer sequence with its coordinates. |
| `run_reference` | one reference release read by one run | |
| `run_input` | one run taking its input from another | The selection applied, such as "sequences the earlier run called viral". |
| `taxon_change` | one change to a taxon between two releases | Old taxon, new taxon, type of change. |

`alignment` ties read-based and contig-based analysis together. A read cluster with only a distant translated hit and no alignment to a contig is the low-identity, unassembled signal that assembly-only catalogs drop.

`taxon_change` makes re-annotation a join. The need is real: the 2022 ICTV release removed the three tailed-phage families that held 1,921 species, and the 2024 release renamed 2,884 species (survey 02).

## Levels and summing counts

Counts are stored once, at the finest level measured. Counts at higher levels are computed by following links:

| From | To | Through |
|---|---|---|
| read cluster | contig | `alignment`, primary alignments only |
| contig | genome | `genome_member` |
| genome | cluster | `cluster_member` |
| read cluster or contig | cluster | `cluster_member`, where a researcher clustered sequences directly |
| any of these | taxon | `taxon_assignment`, preferred rows only |
| taxon | higher rank | the parent column of `taxon` |

The primary flag and the preferred flag exist so that each fragment is counted once. Some numbers cannot be summed, mean depth for example. The manifest marks which can.

## Read clusters

A read cluster is a set of reads or read pairs from one library that a pipeline grouped and represented by one sequence, with a count. Hecatomb is one example: it clusters each sample separately and keeps a count for each representative.

- A read cluster is formed within one library. A researcher who groups read clusters or contigs across libraries records that as a `cluster`, at whatever level they choose.
- The clustering level is the researcher's choice. The standard sets no identity threshold and recommends none. Exact matching, Hecatomb version 1's 97% identity, or any other level is valid. The run records the method, the thresholds and how the representative was chosen, in `run_parameter`.
- The unit is whatever the pipeline clustered. The five read forms below cover it, and a pipeline may mix them.
- Reads on opposite strands are different sequences. They stay separate unless the pipeline's clustering merged them.
- The run declares whether it counted fragments or reads. Within a library the counts sum to the number of fragments or reads that entered clustering, and a check compares that sum with the total in `library`.
- A count is not an abundance. It does not separate PCR copies from independent reads that happen to share a sequence. That distinction needs mapping positions, which `alignment` and `coverage` carry.
- Two bundles clustered at different levels have counts with different meanings. Because the level is recorded, software and readers can see this before combining them.

### Read forms

Paired-end sequencing reads each fragment from both ends. Pipelines differ in what they do with the two reads, and the standard accepts every choice. Each `read_count` row states its read form.

| Read form | What the pipeline did | What is stored |
|---|---|---|
| single | The library has one read per fragment. | one sequence |
| first of pair | Used the first read and set the second aside. Hecatomb version 1 does this. | one sequence |
| merged pair | Joined two overlapping reads into one longer sequence. | one sequence |
| one of pair | Treated the two reads as independent. | one sequence per read, with the end it came from (first or second) on the count row |
| pair | Kept the two reads together, so two fragments match only if both reads match. | one row holding both reads |

For independent reads, the same sequence seen as a first read and as a second read gets two count rows, so the reads can be paired again later. The first four forms are single sequences with ordinary identifiers. Only the pair form needs the two-read row described under "Identifiers". A pipeline that merges the pairs that overlap and keeps the rest as pairs writes both forms into the same tables.

The container performs one merge of its own. Two clusters in the same library, the same read form and the same end whose stored sequences are identical, byte for byte, become one row with the counts added. This follows from naming a sequence by its content.

## Searches

Most taxonomy comes from sequence similarity searches, and pipelines chain several of them. A common pattern searches reads against viral proteins, checks the hits against proteins from all domains of life to remove false positives, and then repeats both steps against nucleotide databases for the reads still unplaced. Hecatomb works this way. The standard holds any such chain, in any order, with any tools.

- Each search is a run, with its tool, settings and reference release.
- Every stored hit is an `alignment` row. Each run's taxon call is a `taxon_assignment` row. A read can be called viral by one step and bacterial by the next, and both rows are kept.
- The order of the chain is data. A `run_input` row says that one run took its input from another and what was selected. Every step names the pipeline run it belongs to.
- The pipeline's final answer for each sequence is the row flagged as preferred.

### Search types

Each alignment states what was compared with what, because identity and coordinates mean different things in each.

| Search type | Query | Target | Identity is measured in |
|---|---|---|---|
| nucleotide to nucleotide | nucleotide | nucleotide | nucleotides |
| translated nucleotide to protein | nucleotide, translated in six frames | protein | amino acids |
| protein to protein | protein | protein | amino acids |
| protein to profile | protein | a profile of a protein family | not applicable |

Other types, such as k-mer classification with no alignment, are declared by the run.

### Which hits are kept

The standard applies no filter of its own. The producer decides what to store and declares it in `run_parameter`: the limits the search itself applied (its e-value ceiling, its maximum hits per query) and whether the bundle holds every reported hit or only the best. "Found nothing" then has a precise meaning: nothing under that ceiling. The declaration uses the setting name `hits_kept`. A search run without it raises a warning.

Keeping weak hits matters. A low-identity hit with a poor e-value can be the only trace of an unknown virus, and a filter applied before storage removes it for every later reader.

E-values depend on the size of the database searched, so they compare only within one run. Bit scores compare across runs. The manifest says this for each column.

## Connectivity

The tables are designed to be joined. A result recorded at one level can be carried to another by following identifiers:

| From | To | Through |
|---|---|---|
| a read | its read cluster | `read`, where kept |
| a read cluster | the contig it maps to, and where on it | `alignment`, with coordinates on the contig |
| a contig | its genome, its cluster | `genome_member`, `cluster_member` |
| any of these | its taxon, from each step of each pipeline | `taxon_assignment` |
| a hit | the database entry and that entry's taxon | `reference_entry` |

Three examples of what this allows:

- Taxa along a contig. Join the read clusters aligned to one contig with their own taxon assignments, ordered by position. If the reads on one end carry species A and those on the other end carry species B, the contig may be a recombinant or a chimeric assembly.
- Context for a dark read. A read cluster that matched nothing, but aligns to a contig that has a taxon, inherits a candidate identity.
- Context for a dark contig. A contig that matched nothing as a whole may be covered by reads with weak protein hits, which point to a distant relative.

For this to work, a run that maps reads to contigs declares what it mapped: read clusters or every read. If it mapped every read, the `read` table connects each read's position to its read cluster, and an alignment row may name a read as its query.

## Dark matter

Sequences that match nothing known, or match only weakly, are often called dark matter. They occur at every level: reads, contigs, genomes and clusters. The standard keeps and counts them like any other sequence.

What a bundle records:

- A sequence needs no label to exist. A read cluster or contig that matched nothing still has its `sequence` row and its counts.
- A search that found nothing is recorded. `taxon_assignment` holds a row with no taxon and `outcome` set to `none`, so the unmatched set of each run is explicit.
- Weak matches keep their statistics. Every stored alignment carries its e-value, bit score and identity, so a low-identity hit is available to anyone who looks, whatever cutoff the original analysis applied.
- The bases travel with the bundle. For a dark sequence the bases are its only identity; without them it cannot be searched again or compared with another study. A read cluster or contig without bases raises a warning.
- A cluster needs no taxon name. Within one study, an unmatched genome found in forty samples is one `cluster` row with forty counts.
- Dark reads and contigs can be grouped across samples, as clusters of sequences.
- Omissions are visible. If a pipeline leaves unmatched reads out of the bundle, its counts fall short of the library total, and the shortfall is the number of fragments left out.

### How dark

Four levels can be told apart by query. The standard fixes no threshold between them.

| Level | How it shows in the tables |
|---|---|
| Named | assigned to a taxon in a formal taxonomy such as ICTV |
| Seen before, unnamed | assigned only to a group in a catalog of uncultivated viruses, or clustered with a reference genome that has no formal name |
| Recurrent in this study | no assignment, but found in more than one sample, by identical sequence, by alignment to a shared contig, or by cluster |
| Seen once | no assignment and found in one sample |

Dark status is relative to a date. Each search is a run against a dated release, so a sequence that matched nothing in 2022 and matches in 2026 has two runs with different answers.

### Across studies

Combining bundles joins their tables. It does not recognise the same unknown virus in two studies. An identifier matches only a byte-identical sequence, and the same genome assembled twice will usually differ by a few bases, by strand or by start position. Recognising it takes a comparison step over the combined bundles: a new clustering run over all the genomes, or assignment of each study's genomes to a shared catalog. The result is stored like any other cluster or label. A bundle keeps what that step needs, but it does not perform it.

## Other data types

A 16S dataset is sequences, a count table and a taxonomy label table. A bacterial metagenome is contigs, genomes, coverage and taxonomy. Both fit the catalog as it stands and enter a bundle as more rows in the same tables, tied to the virome through shared sample identifiers.

## Tables outside the catalog

A result the catalog does not cover goes in a table whose name starts with `x_`. It must be one of the four shapes, carry that shape's required columns, and be described in the manifest. Software that knows only the shapes can still read it: any count table converts to a matrix, any label table to annotation columns, any link table to a graph. Tables that prove useful are candidates for the catalog.

## Later catalog entries

None is specified yet. They will be added in this order.

| Order | Table | Shape | Content |
|---|---|---|---|
| 1 | `viral_call` | label | virus, plasmid or cellular, for a sequence |
| 1 | `quality` | label | completeness and contamination of a genome |
| 2 | `taxon_profile` | count | a taxon in a library, for profilers that report taxa without sequences |
| 3 | `host_prediction` | label or link | a host taxon, or a host genome held in the bundle |
| 3 | `lifestyle` | label | temperate or virulent |
| 4 | `origin_class` | label | infects the sampled host, infects its microbiota, passes through from diet or environment, or reagent contaminant |
| 5 | `feature` | thing | a gene or other interval on a sequence |
| 5 | `function` | label | a protein family for a feature |
| on demand | `genome_similarity` | link | identity and aligned fraction for a pair of genomes |

`origin_class` has no counterpart in any surveyed tool. Version 1 already holds its inputs: the host category on `taxon` and the control type on `sample`.

Within-host variation (variant sites and their frequencies) stays deferred. Strain-level profiles can identify individuals, and no consortium policy found in survey 03 addresses them.

## Manifest

One JSON file at the root of the bundle. For every table it gives the shape, the catalog entry and version (or that the table is outside the catalog), what one row is, the key, a description and unit for each column, the row count and the value ranges. For the bundle it gives an identifier, the version of the standard, and the access tier and data use term of each table.

A person or an AI agent who reads only the manifest knows what every table means and how the tables join. The schema is written in LinkML, the schema language MIxS uses, and the manifest is generated from it.

## Conformance checks

A bundle is checked against the standard, and each check has a severity. The three levels follow pointblank, a data validation package for R and Python.

| Severity | Meaning | Examples |
|---|---|---|
| Critical | The bundle cannot be read reliably. | A required table or column is missing. An identifier appears twice in a thing table. The manifest is absent. |
| Error | The bundle can be read, but an analysis of it would be wrong. | An identifier column points at a row that does not exist. Counts exceed the library total. A count, label or link row names no run. A thing has two preferred rows from one pipeline. |
| Warning | The bundle is usable but incomplete. | Counts fall short of the library total. A read cluster or contig has no bases. A run lists no parameters. A search run does not declare which hits were kept. A reference release has no checksum. A column has no description. |

A bundle conforms when it has no critical failures and no errors. Warnings are reported and do not block. For each failed check the report gives the number of failing rows and can list them.

The checks concern structure only. No check judges an analytical choice.

## Questions the model must answer

These drive the schema and become the test set, for people and for AI agents given only the bundle.

Answerable with the version 1 catalog:

1. What share of each library's fragments comes from plant viruses, by family, after removing taxa found in the negative controls?
2. Restate every family-level assignment made under one ICTV release in a later release.
3. Which read clusters have only translated hits below 40% identity and align to no contig?
4. For one contig, list every run's taxon assignment with its score, and the read clusters that align to it.
5. Sum read counts from read clusters up to species-level clusters, and report how many samples each cluster occurs in.
6. Were two bundles clustered at the same level and searched against the same reference releases?
7. Which tables derive from restricted inputs?
8. What fraction of each library's fragments matched nothing in any search, and which of those lie in contigs or clusters found in more than one sample?
9. Which unmatched sequences recur in more than one sample?
10. Which clusters have no reference member and no named taxon, and in how many samples does each occur?
11. Which sequences that matched nothing under an older reference release match under a newer one?
12. Along one contig, which taxa do the aligned read clusters carry, by position?
13. Which unmatched read clusters align to a contig that has a taxon?
14. Which read clusters did one search step call viral and a later step reassign to a cellular taxon?

Waiting on later catalog entries:

15. Which species-level clusters predicted to infect Bacteroides occur in more than 10% of infants at month 3?
16. What is the ratio of temperate to virulent genomes per sample, counting high-quality genomes only?
17. Does the abundance of a phage follow that of its predicted host in the paired bacterial metagenome?

## Open decisions

None at present. The next step is to agree that scope and schema are settled.

## Decisions made

| Date | Decision | Reason |
|---|---|---|
| 2026-11-20 | The standard stands alone and maps to consortium identifiers through `identifier_map`. | No consortium has published a specification for the content of virome result files. A standard that needs no consortium database can be used anywhere. |
| 2026-11-20 | Read clusters follow Hecatomb: per library, counted in fragments. The standard sets no clustering threshold; the run records the method and thresholds used. | A data container should describe an analysis, not dictate it. |
| 2026-11-20 | Tool-specific detail is a JSON string. | Parquet's semi-structured type is not yet readable by two major query engines (survey 01). |
| 2026-11-20 | The standard is four table shapes plus a catalog of named tables. Each catalog table has typed columns. | Biological data reduces to a few kinds of table. A small core is easier to learn and extend, and typed columns let a wrong query fail. |
| 2026-11-20 | The version 1 catalog is the phyloseq equivalent plus the links between reads, contigs, genomes and clusters. | Start small; add host, lifestyle, function and the rest once the core is proven. |
| 2026-11-20 | Sequence identifiers are plain GA4GH digests of the sequence as written, with no strand rule. | A strand rule merged a read with its reverse complement whatever the researcher chose. Reverse complements are recorded as links. |
| 2026-11-20 | A read cluster stores whatever unit the pipeline clustered, including both reads of a pair. | Storing one read per pair merged fragments the pipeline had kept apart. |
| 2026-11-20 | Byte-identical read clusters in one library are merged. | It follows from naming a sequence by its content. |
| 2026-11-20 | The unit of exchange is one bundle per study. A consortium store combines bundles. | A single lab can produce and share a bundle without any consortium infrastructure. |
| 2026-11-20 | The standard serves virome researchers outside HVP as well. | Adoption beyond one consortium is the aim. |
| 2026-11-20 | The name is Capsomer. | The format is built from four repeating table shapes, as a capsid is built from capsomers. |
| 2026-11-20 | All five read forms are accepted, and a read pair kept as a unit is one row holding both reads. | Researchers handle paired reads in different ways, and the container must hold each without merging. |
| 2026-11-20 | The run declares whether counts are fragments or reads. | A pipeline that treats the two reads independently counts reads. |
| 2026-11-20 | Sequences with no match are kept, and a search that found nothing writes a row. | Dark matter must be countable, and "searched and found nothing" must differ from "never searched". |
| 2026-11-20 | Conformance checks carry three severity levels, after pointblank. | A shortfall in counts and a missing table are different kinds of failure. |
| 2026-11-20 | Later tables are added in the order listed under "Later catalog entries". | Viral call and quality are the most widely produced; features and functions are the largest and most varied. |
| 2026-11-20 | Independent reads record which end of the fragment they came from. | The reads can then be paired again later. |
| 2026-11-20 | Read pair identifiers use the prefix `RP.` and the namespace `capsomer:`. | They must not be mistaken for GA4GH sequence digests. |
| 2026-11-20 | The specification and the software are under the MIT licence. | It is short, permissive, and the licence Hecatomb uses. |
| 2026-11-20 | Clusters may group sequences as well as genomes, across libraries. | Dark reads and contigs need grouping across samples, and the level is the researcher's choice. |
| 2026-11-20 | Groups in catalogs of uncultivated viruses are taxa, and reference genomes can be cluster members. | "Seen before but unnamed" must differ from "never seen". |
| 2026-11-20 | Missing bases raise a warning. | The bases are a dark sequence's only identity, but requiring them outright would exclude some datasets. |
| 2026-11-20 | The sequence table names optional GC and complexity columns. | They separate low-complexity artifacts from unknown biology. |
| 2026-11-20 | Hecatomb is an example, not a model the standard adheres to. | The standard must fit any pipeline. |
| 2026-11-20 | Search chains are stored as runs joined by `run_input`, with every step's hits and calls kept and the final call flagged. | Pipelines mix translated and untranslated searches in many orders, and intermediate calls are evidence. |
| 2026-11-20 | Database entries that were matched have their own table, `reference_entry`. | Top hit and lowest common ancestor both need the taxon of each hit. |
| 2026-11-20 | Each alignment states its search type. | Amino acid identity and nucleotide identity are not comparable. |
| 2026-11-20 | The producer declares which hits were kept; the standard applies no filter. | Weak hits can be the only trace of an unknown virus, and the choice belongs to the researcher. |
| 2026-11-20 | Connectivity is a rule: every result can be followed to the things it came from. | Read-level information carried onto contigs supports analyses such as detecting chimeric or recombinant contigs. |
| 2026-11-20 | The table of original reads, `read`, is in version 1 as an optional thing table. | When every read is mapped to contigs, it is the only route from a mapping position back to read-level taxonomy. |
