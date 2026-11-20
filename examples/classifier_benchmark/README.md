# Comparison of five viral classifiers

Which tool should decide that a contig is viral? This folder compares five on sequences of known origin, with the scoring done as queries on a Capsomer dataset. Run on 2026-11-20.

| Tool | Version | How it decides |
|---|---|---|
| geNomad | 1.13.0, database 1.9 | marker genes plus a neural network |
| Jaeger | 1.26.5, model `jaeger_448931d3_1.3M_fragment` | neural network on the sequence alone, no reference search |
| VirSorter2 | 2.2.4 | gene features and a classifier per viral group |
| PhaBOX2 (PhaMer) | 2.2, database v2_2 | protein matches fed to a language model |
| Cenote-Taker 3 | 3.4.4, database v3.1.1 | presence of virus hallmark genes |

All were run with their default settings.

## The test sequences

Every test sequence is a 5 kb piece, about the size of a short contig. There are three kinds.

| Kind | Pieces | Where they come from |
|---|---|---|
| Known viruses | 483 | RefSeq viruses used by a ViroForge gut virome simulation |
| Recent viruses | 1,522 | 600 viral genomes deposited at NCBI on or after 2026-10-01 |
| Not viral | 390 | real bacterial, archaeal, fungal and human fragments that ship with ViroForge |

The known viruses cannot show which tool finds new viruses, because every tool was built from RefSeq. The recent viruses address that. Each was made public after the database or model date of all five tools, so none of them can have been built from it. They are of two sorts:

- 300 cultured phages. Someone grew these, so no classifier selected them.
- 300 uncultivated genomes from metagenomes. These are the most novel, but the depositing study picked each one out with some classifier, which favours that classifier here. At most 40 genomes were taken from any one naming pattern, and one study (names starting `KCR_vOTU`) still supplies most of them.

Each recent piece also has a novelty measure: how much of it has protein-level similarity to any RefSeq virus (a translated search with MMseqs2, e-value 1e-5 or better).

| Novelty | Definition | Pieces |
|---|---|---|
| Close to RefSeq | half or more of the piece matches at 50% amino acid identity or better | 988 |
| Distant | half or more matches, but below 50% identity | 183 |
| Very distant | less than half matches anything | 351 |

## Results

Share of viral pieces each tool called viral:

| Tool | Known viruses | Recent, close | Recent, distant | Recent, very distant |
|---|---|---|---|---|
| Jaeger | 99% | 100% | 96% | 86% |
| geNomad | 98% | 99% | 94% | 79% |
| VirSorter2 | 87% | 95% | 91% | 78% |
| PhaBOX2 | 93% | 100% | 93% | 44% |
| Cenote-Taker 3 | 42% | 52% | 68% | 52% |

Non-viral pieces called viral, of 390:

| Tool | Called viral | Called by this tool alone |
|---|---|---|
| Cenote-Taker 3 | 9 | 1 |
| geNomad | 18 | 1 |
| Jaeger | 20 | 4 |
| VirSorter2 | 27 | 9 |
| PhaBOX2 | 82 | 67 |

Thirteen bacterial pieces were called viral by four or five tools. Real bacterial chromosomes contain prophages, so those are probably correct calls on mislabelled pieces. The calls a tool made alone are the better estimate of its errors.

Time for the 873 known and non-viral pieces on an 18-core Apple Silicon Mac, 10 threads:

| Tool | Seconds |
|---|---|
| PhaBOX2 | 16 |
| Jaeger | 18 |
| Cenote-Taker 3 | 114 |
| geNomad | 131 |
| VirSorter2 | 804 |

VirSorter2 and Cenote-Taker 3 ran under Intel emulation, which slows them.

Two tools together: on the very distant pieces, geNomad alone found 79%. Adding Jaeger's calls raised that to 90%, and all five tools together reached 96%. On the non-viral pieces, geNomad or Jaeger called 25 viral, of which 13 are the probable prophages.

## What we take from it

- geNomad and Jaeger are the two best tools here, and they err in different places because they work differently. Jaeger finds more of the very distant viruses. geNomad makes fewer calls on non-viral DNA and also reports taxonomy, plasmids and proviruses.
- VirSorter2 found fewer viruses than geNomad at every novelty level, made more errors, and took six to nine times as long.
- PhaBOX2's identification step called 67 non-viral pieces viral on its own and was weakest on the very distant viruses.
- Cenote-Taker 3 requires a hallmark gene on every contig. Many 5 kb pieces have none, so it misses about half of them. It may do better on long contigs, which this test does not include.

## Limits

- All pieces are 5 kb. Longer contigs are easier for every tool, and the ranking could change.
- The non-viral set is small, and its bacterial pieces include prophages that nobody has labelled.
- Most of the uncultivated genomes come from one study. If that study used one of these tools to find them, that tool is favoured.
- A genome deposited after a tool's database date can still have close relatives that were public earlier. The novelty measure, not the deposit date, is what separates new from familiar.
- Default settings only. No score cutoff was tuned.
- Jaeger's documentation names a model that its download list does not offer. The model used was the most recent in the list.

## Files

- `make_benchmark_set.py` cuts the known and non-viral pieces from a ViroForge dataset.
- `fetch_recent_viruses.py` samples genomes deposited at NCBI since a date, and `cut_pieces.py` cuts them.
- `run_classifiers.sh` runs the five tools on a FASTA file and times them.
- `build_benchmark.py` puts the pieces, their known origin and every tool's calls into one Capsomer dataset. The known origin goes in a table outside the catalog, `x_truth`.
- `benchmark.sql` holds the queries behind the tables above. Run it with `capsomer sql bundle --file benchmark.sql`.

The similarity search behind the novelty measure:

```sh
mmseqs easy-search recent/genomes.fasta refseq_viral.fna novelty/recent_vs_refseq.m8 tmp \
  --search-type 2 -e 1e-5 --max-seqs 50 \
  --format-output "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen"
```

`refseq_viral.fna` was the 14,568 RefSeq viral genomes in ViroForge's local genome store. The sequences, tool outputs and dataset are not kept in the repository.
