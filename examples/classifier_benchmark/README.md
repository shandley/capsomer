# Comparison of five viral classifiers

Which tool should decide that a contig is viral? This folder compares five on sequences of known origin, with the scoring done as queries on a Capsomer dataset. Run on 2026-11-20. The result tables are in `results/`.

| Tool | Version | How it decides |
|---|---|---|
| geNomad | 1.13.0, database 1.9 | marker genes plus a neural network |
| Jaeger | 1.26.5, model `jaeger_448931d3_1.3M_fragment` | neural network on the sequence alone, no reference search |
| VirSorter2 | 2.2.4 | gene features and a classifier per viral group |
| PhaBOX2 (PhaMer) | 2.2, database v2_2 | protein matches fed to a language model |
| Cenote-Taker 3 | 3.4.4, database v3.1.1 | presence of virus hallmark genes |

All were run with their default settings unless a section says otherwise.

## Summary

- No tool is best at every contig length. Jaeger finds the most viruses on pieces of 5 kb and shorter. geNomad finds the most at 10 kb.
- geNomad and Jaeger together find more than either alone at every length tested, at the cost of Jaeger's extra errors, which are highest on short pieces.
- VirSorter2, PhaBOX2 and Cenote-Taker 3 each did worse than geNomad on new viruses, on errors, or on both.

## The test sequences

Every test sequence is a piece of fixed length cut from a sequence of known origin. The main test uses 5 kb pieces.

| Set | 5 kb pieces | Genomes | Where they come from |
|---|---|---|---|
| Known viruses | 483 | 183 | RefSeq viruses used by a ViroForge gut virome simulation |
| Recent, cultured | 892 | 300 | phages deposited at NCBI on or after 2026-10-01 |
| Recent, uncultivated | 630 | 299 | viral genomes from metagenomes, deposited on or after 2026-10-01 |
| Earlier, uncultivated | 582 | 381 | viral genomes from metagenomes, deposited 2026-04-01 to 2026-09-30 |
| Not viral | 390 | 197 | real bacterial, archaeal, fungal and human fragments that ship with ViroForge |

The known viruses cannot show which tool finds new viruses, because every tool was built from RefSeq. The other viral sets address that:

- Recent genomes became public after the database or model date of all five tools, so none of the tools can have been built from them.
- Cultured phages were found by growing them, so no classifier selected them. Most turn out to be close relatives of RefSeq phages.
- Uncultivated genomes are the most novel. Each was picked out of a metagenome by some classifier in the depositing study, which favours that classifier here.
- The recent uncultivated set is dominated by one deposit, a cucumber rhizosphere study (BioProject PRJNA1272261, isolate names starting `KCR_vOTU`). Its records do not say which classifier found the genomes. The first sampling capped genomes per naming pattern, which let this deposit through under several taxon names.
- The earlier uncultivated set was added to correct that. It takes at most 25 genomes from any one deposit and covers 122 deposits. Its genomes postdate the databases of four tools. Jaeger's training date is not documented, so Jaeger may have seen some of them.

Each recent or earlier piece has a novelty measure: how much of it has protein-level similarity to any RefSeq virus (translated search with MMseqs2, e-value 1e-5 or better, against 14,568 RefSeq viral genomes).

| Novelty | Definition |
|---|---|
| Close to RefSeq | half or more of the piece matches at 50% amino acid identity or better |
| Distant | half or more matches, but below 50% identity |
| Very distant | less than half matches anything |

## Results on 5 kb pieces

Share of viral pieces each tool called viral. Recent sets (cultured and uncultivated together):

| Tool | Known viruses (483) | Close (988) | Distant (183) | Very distant (351) |
|---|---|---|---|---|
| Jaeger | 99% | 100% | 96% | 86% |
| geNomad | 98% | 99% | 94% | 79% |
| VirSorter2 | 87% | 95% | 91% | 78% |
| PhaBOX2 | 93% | 100% | 93% | 44% |
| Cenote-Taker 3 | 42% | 52% | 68% | 52% |

Earlier uncultivated set, 122 deposits:

| Tool | Close (189) | Distant (203) | Very distant (190) |
|---|---|---|---|
| Jaeger | 98% | 86% | 85% |
| geNomad | 97% | 85% | 72% |
| PhaBOX2 | 91% | 93% | 41% |
| Cenote-Taker 3 | 70% | 63% | 36% |
| VirSorter2 | 15% | 45% | 52% |

This set includes many RNA and eukaryotic viruses, which VirSorter2's default groups (double-stranded DNA phage and single-stranded DNA) do not target.

Does the ranking depend on one study? Very distant pieces of the recent sets, by deposit:

| Deposit | Pieces | Jaeger | geNomad | VirSorter2 |
|---|---|---|---|---|
| KCR rhizosphere | 247 | 82% | 74% | 74% |
| WSV | 61 | 93% | 92% | 90% |
| Cultured phages | 27 | 100% | 96% | 89% |
| Other uncultivated | 16 | 88% | 81% | 81% |

In the earlier set, eight deposits have at least five very distant pieces. Jaeger found more than geNomad in five of them, the two tied in three, and geNomad led in none. So Jaeger's lead over geNomad on 5 kb pieces is consistent across deposits, though most of its size in the recent set comes from the KCR deposit.

Non-viral pieces called viral, of 390:

| Tool | Called viral | Called by this tool alone |
|---|---|---|
| Cenote-Taker 3 | 9 | 1 |
| geNomad | 18 | 1 |
| Jaeger | 20 | 4 |
| VirSorter2 | 27 | 9 |
| PhaBOX2 | 82 | 67 |

Thirteen bacterial pieces were called viral by four or five tools. Real bacterial chromosomes contain prophages, so those are probably correct calls on mislabelled pieces. The calls a tool made alone are the better estimate of its errors.

## Results at other lengths

Jaeger's model was trained on 5,000-base fragments, the same length as the main test, so the three fast tools were also run on 2 kb and 10 kb pieces of the known, recent and non-viral sets.

Share of very distant pieces called viral:

| Piece length | Pieces | geNomad | Jaeger | Either |
|---|---|---|---|---|
| 2 kb | 313 | 24% | 81% | 83% |
| 5 kb | 351 | 79% | 86% | 90% |
| 10 kb | 185 | 90% | 84% | 94% |

Share of non-viral pieces called viral, probable prophages included:

| Piece length | Pieces | geNomad | Jaeger | Either | PhaBOX2 |
|---|---|---|---|---|---|
| 2 kb | 197 | 1% | 7% | 7% | none called |
| 5 kb | 390 | 5% | 5% | 6% | 21% |
| 10 kb | 193 | 4% | 5% | 6% | 24% |

- geNomad's default rule asks for a hallmark gene on sequences under 2.5 kb, which is why it finds few 2 kb pieces of any kind (30% of known viruses).
- PhaBOX2 skips sequences under 3 kb by default.

Is geNomad with relaxed filters (`--relaxed`) a substitute for adding Jaeger?

| Piece length | Measure | geNomad default | geNomad relaxed | Jaeger | geNomad default or Jaeger |
|---|---|---|---|---|---|
| 2 kb | very distant found | 24% | 71% | 81% | 83% |
| 2 kb | non-viral called | 1% | 19% | 7% | 7% |
| 10 kb | very distant found | 90% | 96% | 84% | 94% |
| 10 kb | non-viral called | 4% | 9% | 5% | 6% |

It is not. On short pieces the relaxed setting finds fewer very distant viruses than Jaeger and makes nearly three times the false calls. At 10 kb it finds slightly more than the pair, with more false calls.

Jaeger's reliability score does not separate its correct calls from its errors: raising the cutoff from 0 to 0.95 removed none of the 11 non-viral 2 kb pieces that only Jaeger called.

## Run time

Seconds on an 18-core Apple Silicon Mac, 10 threads, for the 873 known and non-viral 5 kb pieces:

| Tool | Seconds |
|---|---|
| PhaBOX2 | 16 |
| Jaeger | 18 |
| Cenote-Taker 3 | 114 |
| geNomad | 131 |
| VirSorter2 | 804 |

VirSorter2, Cenote-Taker 3 and PhaBOX2 ran under Intel emulation, which slows them. `results/timing.tsv` has every run.

## What we take from it

- geNomad and Jaeger are the two best tools here, and they are strong at different lengths because they work differently. Jaeger needs no genes, so it handles short contigs. geNomad gains from every extra gene, and also reports taxonomy, plasmids and proviruses.
- Used together with default settings, they found 83% to 94% of very distant pieces across the three lengths, with false calls on 6% to 7% of non-viral pieces (some of which are prophages).
- Jaeger's contribution is concentrated in contigs under about 5 kb. That is also where it errs most. In a sample with much bacterial DNA, a minimum length for contigs accepted on Jaeger's call alone is advisable.
- VirSorter2 found fewer viruses than geNomad at every novelty level, made more errors, and took six to nine times as long.
- PhaBOX2's identification step called a fifth to a quarter of non-viral pieces viral and was weak on very distant viruses.
- Cenote-Taker 3 requires a hallmark gene on every contig, so it misses about half of 5 kb pieces. It was not tested on longer pieces.

## Limits

- Pieces of fixed length are not contigs. Real contigs vary in length and contain assembly errors.
- The non-viral set is small, and its bacterial pieces include prophages that nobody has labelled. Its 2 kb and 10 kb versions hold about 195 pieces each.
- Each uncultivated genome was found by a classifier we cannot identify. The study names are also only naming patterns, not BioProjects.
- A genome deposited after a tool's database date can still have close relatives that were public earlier. The novelty measure, not the deposit date, is what separates new from familiar.
- Jaeger may have been trained on some genomes in the earlier set.
- The viral sets mix phages with eukaryotic and RNA viruses. A phage-only analysis has not been done.
- Default settings only, apart from the geNomad relaxed run. No score cutoff was tuned.
- VirSorter2 and Cenote-Taker 3 were run at 5 kb only.
- Jaeger's documentation names a model that its download list does not offer. The model used was the most recent in the list.

## Files

- `make_benchmark_set.py` cuts the known and non-viral pieces from a ViroForge dataset.
- `fetch_recent_viruses.py` samples genomes deposited at NCBI in a date window, and `cut_pieces.py` cuts them.
- `run_classifiers.sh` runs the tools on a FASTA file and times them.
- `build_benchmark.py` puts the 5 kb pieces, their known origin and every tool's calls into one Capsomer dataset. The known origin goes in a table outside the catalog, `x_truth`.
- `benchmark.sql` holds the queries behind the 5 kb tables. Run it with `capsomer sql bundle --file benchmark.sql`.
- `score_length.py` and `score_relaxed.py` score the 2 kb and 10 kb runs.
- `export_results.py` writes `results/calls.tsv`.

In `results/`:

- `calls.tsv`: one row per 5 kb piece, with its set, origin, source accession, study, novelty measures, and each tool's call and viral score.
- `genomes_recent.tsv` and `genomes_earlier.tsv`: the NCBI accessions sampled, with their titles.
- `benchmark_queries_output.txt`: the output of `benchmark.sql`.
- `length_2000.tsv`, `length_10000.tsv`, `relaxed_2000.tsv`, `relaxed_10000.tsv`: the tables for other lengths.
- `timing.tsv`: run time of every tool on every set.

The similarity search behind the novelty measure:

```sh
mmseqs easy-search genomes.fasta refseq_viral.fna recent_vs_refseq.m8 tmp \
  --search-type 2 -e 1e-5 --max-seqs 50 \
  --format-output "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen"
```

`refseq_viral.fna` was the 14,568 RefSeq viral genomes in ViroForge's local genome store. The sequences, raw tool outputs and the Capsomer dataset are not kept in the repository.
