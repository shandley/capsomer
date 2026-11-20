# Worked example: contig tools on simulated reads

This example tests whether the catalog holds the results of a contig-based virome analysis made with widely used tools, none of which is Hecatomb.

The reads come from [ViroForge](https://github.com/shandley/viroforge), which simulates a virome sequencing library from known genomes and records the true source of every read. The dataset is a skin virome with host DNA, ribosomal RNA, bacterial and retroviral contaminants: 257,505 read pairs.

## What was run

`run_tools.sh` holds the exact commands and versions.

| Step | Tool | Result |
|---|---|---|
| Read trimming | fastp | 512,810 of 515,010 reads kept |
| Assembly | MEGAHIT | 134 contigs of 1 kb or longer |
| Assembly | metaSPAdes | 100 contigs of 1 kb or longer |
| Viral identification | geNomad | 77 MEGAHIT contigs and 53 metaSPAdes contigs called viral |
| Genome quality | CheckV | 27 of the 130 viral contigs rated high-quality |
| Species-level clustering | vclust, 95% identity over 85% of the shorter sequence | 83 clusters |
| Looser clustering | MMseqs2, 70% identity | 38 clusters |
| Read mapping | minimap2 to the 83 cluster representatives | 97.2% of reads mapped |
| Coverage | CoverM | one row per representative |
| Truth alignment | minimap2, every contig against the simulated genomes | every contig placed |

`build_bundle.sh` imports the ViroForge ground truth and then adds each tool's results with `capsomer add`. The bundle conforms to the standard, with one warning (the reference releases have no checksums). It is 15 MB and is not kept in the repository.

## What the bundle can answer

`questions.sql` joins the tool results to the truth. Run it with `capsomer sql bundle --file questions.sql`.

- **Viral calls against the truth.** Of the 77 MEGAHIT contigs geNomad called viral, 73 came from simulated viruses, 3 from exogenous retroviruses and 1 from the phiX spike-in. Of the 57 it did not call, 39 came from simulated viruses.
- **Completeness against the truth.** CheckV's estimate tracks the true fraction of the genome each contig covers: 0.98 estimated against 0.97 true for high-quality MEGAHIT contigs, and 0.11 against 0.10 for low-quality ones.
- **Two assemblers, one sequence.** Two contigs were assembled byte for byte by both MEGAHIT and metaSPAdes. Each is one `sequence` row with two `sequence_source` rows. This case is why a sequence row carries no run and no name.
- **Clusters across assemblies.** 38 of the 83 species-level clusters hold contigs from both assemblers.
- **Reads on a contig, by true genome.** Each read's position on a contig comes from the mapper and its true genome from the simulator. Several contigs carry reads from eight different simulated genomes, which shows that close relatives collapsed into one contig during assembly.
- **Abundance against the truth.** For the most abundant complete genome, 86,271 reads mapped to its contig and 86,634 were simulated from it.

## What this example does not test

- geNomad found no proviruses and no terminal repeats in these contigs, and CheckV found no host regions. The importer code for those cases is covered by unit tests on made-up files, not by real output.
- Host reads were not removed before assembly.
- `viral_call` and `quality` each have one real tool behind them. A second classifier and a second quality tool still need to be imported.
- The dataset is one sample. About half of its viral reads come from tailed phages and most of the rest from papillomaviruses and polyomaviruses. A multi-sample gut virome would test clustering and abundance across samples.
