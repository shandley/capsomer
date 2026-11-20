# Analysis log

A record of every analysis done while building Capsomer and capsomer-workflow, kept so that results can be traced and reused in a publication. Each entry gives the question, the data, the method, the result, its limits, and where the files are. Entries are added as work is done and are not rewritten afterwards; a later correction gets its own line under "Corrections".

Conventions:

- "This repo" is github.com/shandley/capsomer. "Workflow repo" is github.com/shandley/capsomer-workflow.
- Commit hashes name the state of the code that produced a result.
- Sequences and raw tool outputs are not in either repository. They are in git-ignored `work/` folders on the machine that ran the analysis (an 18-core Apple Silicon Mac, macOS). Result tables small enough to keep are in `results/` folders and are tracked.
- Analyses of unpublished research data are recorded in a separate private log, not here.

## Tool and database versions used on 2026-11-20

| Tool | Version | Database or model |
|---|---|---|
| fastp | 1.4.0 | |
| MEGAHIT | 1.2.9 (one thread; it crashes with more on Apple Silicon) | |
| SPAdes (metaSPAdes) | 4.3.0 | |
| geNomad | 1.13.0 | database 1.9 |
| CheckV | 1.1.1 | database 1.5 |
| vclust | 1.3.1 | |
| MMseqs2 | 19.ab57c | |
| minimap2 | 2.31-r1302 | |
| samtools | 1.25 | |
| CoverM | 0.8.0 | |
| VirSorter2 | 2.2.4 (Intel emulation) | database fetched by `virsorter setup` on 2026-11-20 |
| CheckM2 | 1.1.0 | database from Zenodo record 14897628 |
| Cenote-Taker 3 | 3.4.4 (Intel emulation) | database v3.1.1 |
| PhaBOX | 2.2 from GitHub tag v2.2 (Intel emulation) | phabox_db_v2_2 |
| Jaeger | 1.26.5 | model jaeger_448931d3_1.3M_fragment |
| ViroForge | 0.20.0, commit 2363401 | |
| Snakemake | 9.27.2 | |

## 1. Do real tool outputs fit the catalog, and do they agree with known truth?

- **Date:** 2026-11-20. Code: this repo at 7972d8b.
- **Question:** Can the results of a contig-based analysis by widely used tools be imported without bending the catalog, and what do they look like against ground truth?
- **Data:** one ViroForge skin virome (257,505 read pairs) with known source for every read.
- **Method:** fastp, MEGAHIT, metaSPAdes, geNomad, CheckV, vclust, MMseqs2, minimap2 and CoverM, then import of the ViroForge truth and each tool's output into one dataset. Commands in `examples/viroforge_skin/run_tools.sh` and `build_bundle.sh`; queries in `examples/viroforge_skin/questions.sql`.
- **Result:** the dataset conforms. MEGAHIT gave 134 contigs of 1 kb or more and metaSPAdes 100. geNomad called 77 and 53 of them viral. Of the 77 MEGAHIT calls, 73 came from simulated viruses; of the 57 uncalled, 39 were viral. CheckV completeness tracked true completeness (0.98 against 0.97 for high-quality contigs). Two contigs were byte-identical between the assemblers. 38 of 83 species-level clusters held contigs from both assemblers. Several contigs carried reads simulated from eight different genomes.
- **Consequence for the standard:** a sequence can be produced by more than one run, so run and name moved from `sequence` to a new table, `sequence_source`; `run_library` was added (commit 9c6343d).
- **Limits:** one sample; no proviruses or terminal repeats occurred, so those import paths are covered by unit tests only; host reads were not removed.
- **Files:** write-up in `examples/viroforge_skin/README.md`.

## 2. How much of the catalog does each data source fill?

- **Date:** 2026-11-20. Code: this repo at 0707c3d.
- **Question:** Is the catalog shaped by Hecatomb, the first importer?
- **Method:** count catalog tables and columns with at least one non-empty value in each example dataset's manifest.
- **Result:** of 27 tables and 206 columns, the ViroForge truth plus contig tools fill 23 tables and 133 columns; a Hecatomb import fills 15 and 93; together 25 and 155. Only `read_count` and `reference_entry`, both read-level search tables, are filled by Hecatomb alone. `sequence_relation` and `taxon_change` were filled by nothing.
- **Limits:** the Hecatomb figures come from unpublished data and are given here only as counts of columns.

## 3. Second tools for `viral_call` and `quality`

- **Date:** 2026-11-20. Code: this repo at 6ef1ecd.
- **Question:** Does each new catalog table hold the output of a second, independent tool without change?
- **Data:** the skin sample's contigs for VirSorter2; the three genomes that ship with CheckM2 as test data, plus one deliberately mixed bin of two of them.
- **Result:** both imported with no schema change. On the 134 MEGAHIT contigs geNomad and VirSorter2 both called 68 viral, VirSorter2 alone 33 (32 truly viral) and geNomad alone 9. CheckM2 gave completeness 0.99 to 1.00 and contamination 0.002 to 0.008 for the three genomes and contamination 0.36 for the mixed bin.
- **Limits:** CheckM2 scores bacterial genomes; a second quality tool for viruses has not been imported. The CheckM2 dataset has no taxa or counts and so does not conform.

## 4. First end-to-end run of the workflow

- **Date:** 2026-11-20. Code: workflow repo at de04da5.
- **Data:** two ViroForge skin viromes (`skin_erv`, `skin_test`).
- **Result:** 19 of 19 steps completed and the dataset conforms. 134 and 111 contigs; 77 and 72 called viral by geNomad; one contig assembled identically in both samples; 98 species-level clusters, 95 with reads from both samples.

## 5. Maintenance status of VirSorter2 and CheckM2, and alternatives

- **Date:** 2026-11-20. Method: GitHub API queries and web search.
- **Result:** VirSorter2's last release was January 2023 and its last code push November 2024; it pins scikit-learn 0.22.1, which has no Apple Silicon build. CheckM2's last release was February 2025; its macOS failure is issue 97, with an unmerged fix (pull request 164, October 2026). Maintained alternatives exist for viral identification (Cenote-Taker 3, PhaBOX2, Jaeger) and viral quality (ViralQC).
- **Decision:** neither tool is forked, and neither is in the workflow.

## 6. Comparison of five viral classifiers

- **Date:** 2026-11-20. Code: this repo at 196dcbb for the first version; the commit that adds this log for the extended version.
- **Question:** Which classifier, or pair, should the workflow use to decide that a contig is viral, with attention to viruses unlike those in reference databases?
- **Data:** 5 kb pieces of RefSeq viruses from a ViroForge gut simulation (483), of real non-viral DNA (390), of 600 viral genomes deposited at NCBI on or after 2026-10-01 (1,522), and of 381 uncultivated viral genomes from 122 deposits made 2026-04-01 to 2026-09-30 (582). Accessions are in `examples/classifier_benchmark/results/`.
- **Method:** geNomad, Jaeger, VirSorter2, PhaBOX2 and Cenote-Taker 3 with default settings; a novelty measure for each recent piece from a translated MMseqs2 search against 14,568 RefSeq viral genomes; scoring by SQL on a Capsomer dataset. Then the three fast tools on 2 kb and 10 kb pieces, and geNomad with `--relaxed`.
- **Result:** in `examples/classifier_benchmark/README.md`, with every table as a file. In brief: Jaeger found the most very distant viruses at 2 kb and 5 kb (81% and 86%, against 24% and 79% for geNomad); geNomad found the most at 10 kb (90% against 84%); the two together found 83%, 90% and 94%. Jaeger called 7% of non-viral 2 kb pieces viral. geNomad with relaxed filters was not a substitute.
- **Decision:** the workflow takes geNomad's viral calls plus contigs only Jaeger calls viral, with an optional minimum length for the latter (workflow repo, 701f597 and bfcd9bc).
- **Limits:** listed in the README. The main ones: fixed-length pieces are not contigs; the non-viral set is small and contains unlabelled prophages; the classifier that found each uncultivated genome is unknown; Jaeger's training date is unknown.

## 7. Effect of the two-classifier rule on a real workflow run

- **Date:** 2026-11-20. Code: workflow repo at 701f597.
- **Data:** the two skin viromes of entry 4; truth for `skin_erv` from entry 1, joined by sequence identifier across the two datasets.
- **Result:** of 244 distinct contigs, 110 were called viral by both tools, 38 by geNomad only and 9 by Jaeger only. For `skin_erv`, truly viral contigs recovered rose from 73 to 80 of 112, and no non-viral contig was added. CheckV rated all nine Jaeger-only contigs "Not-determined".
- **Limits:** these samples contain only RefSeq viruses.

## 8. Defects found in ViroForge

- **Date:** 2026-11-20. ViroForge 0.20.0 at commit 2363401, `gut-standard` preset, seed 11.
- **Result:** three defects, each filed on github.com/shandley/viroforge:
  - Issue 90: PCR-duplicate read names (`.../1_dup2`) break read pairing; 347,486 of 2,144,615 first reads carry such names.
  - Issue 91: reagent bacteria are random sequence when no database path is given; CheckM2 finds no protein matches in them.
  - Issue 92: background contaminants repeat; 200 bacterial records hold 27 distinct sequences.
- **Consequence here:** reagent bacteria are left out of the classifier test, and the test draws non-viral pieces from ViroForge's full reference files, not from one dataset.

## Corrections

- 2026-11-20: an early statement that the skin sample was "mostly papillomaviruses and polyomaviruses" was wrong. About half of its viral reads come from tailed phages (125,082 of about 251,000 viral fragments).
- 2026-11-20: the first write-up of entry 6 said Jaeger finds the most new viruses. That held only for 5 kb pieces, which is the fragment length Jaeger was trained on, and rested largely on one deposit. The extended analysis replaced it with the length-dependent result above.
- 2026-11-20: the first sampling of recent genomes capped genomes per naming pattern, which let one deposit through under several taxon names. The cap is now per deposit, and a second set of 122 deposits was added.
