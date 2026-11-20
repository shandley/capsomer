#!/usr/bin/env bash
# Build the example bundle: ViroForge ground truth, then the results of each tool.
# Run run_tools.sh first. Usage: build_bundle.sh <ViroForge dataset folder>
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
dataset="${1:?give the ViroForge dataset folder}"
work="$here/work/contig"
bundle="$here/bundle"
lib=skin_erv
cap() { uv run --project "$here/../.." capsomer "$@"; }

rm -f "${bundle:?}"/*.parquet "${bundle:?}"/manifest.json
cap import-viroforge "$dataset" "$bundle" --sample-id "$lib" > /dev/null

cap add fastp "$bundle" "$work/qc/fastp.json" --run-id fastp --library "$lib" --input viroforge

for assembler in megahit metaspades; do
  version=$( [ "$assembler" = megahit ] && echo 1.2.9 || echo 4.3.0 )
  cap add contigs "$bundle" "$work/$assembler/contigs.fasta" --run-id "$assembler" \
    --assembler "$assembler" --tool-version "$version" --library "$lib" --input fastp \
    --parameter min_contig_length=1000 --parameter "assembler_depth=k-mer depth from the contig header"
  cap add genomad "$bundle" "$work/$assembler/genomad" "$work/$assembler/contigs.fasta" \
    --run-id "genomad.$assembler" --tool-version 1.13.0 --database-version 1.9 --input "$assembler"
  cap add checkv "$bundle" "$work/$assembler/checkv" \
    "$work/$assembler/genomad/contigs_summary/contigs_virus.fna" \
    --run-id "checkv.$assembler" --tool-version 1.1.1 --database-version 1.5 --input "genomad.$assembler"
done

cap add clusters "$bundle" "$work/votu/clusters.tsv" "$work/votu/viral_sequences.fna" \
  --run-id vclust --clusterer vclust --tool-version 1.3.1 --level species --layout vclust \
  --input genomad.megahit --input genomad.metaspades \
  --parameter algorithm=leiden --parameter metric=ani --parameter min_ani=0.95 \
  --parameter min_query_coverage=0.85 --parameter "members=sequences geNomad called virus"
cap add clusters "$bundle" "$work/mmseqs/genus_cluster.tsv" "$work/votu/viral_sequences.fna" \
  --run-id mmseqs --clusterer "MMseqs2 easy-cluster" --tool-version 19.ab57c --level "70% identity" \
  --layout mmseqs --input genomad.megahit --input genomad.metaspades \
  --parameter min_seq_id=0.7 --parameter coverage=0.5 --parameter cov_mode=1 \
  --parameter "members=sequences geNomad called virus"

cap add read-mapping "$bundle" "$work/votu/representatives.fna" --run-id minimap2.reads \
  --mapper minimap2 --tool-version 2.31-r1302 --library "$lib" \
  --alignments "$work/map/reads_to_representatives.bam" --input fastp --input vclust \
  --parameter preset=sr --parameter "targets=representatives of the vclust species clusters"
cap add coverm "$bundle" "$work/map/coverm.tsv" "$work/votu/representatives.fna" \
  --run-id coverm --tool-version 0.8.0 --library "reads_to_representatives=$lib" \
  --input minimap2.reads --parameter mapped_unit=reads \
  --parameter "methods=count mean covered_fraction covered_bases"

cap add paf "$bundle" "$work/map/contigs_to_truth.paf" "$work/map/all_contigs.fasta" \
  "$(ls "$dataset"/fasta/*.fasta)" --run-id minimap2.truth --aligner minimap2 \
  --tool-version 2.31-r1302 --input megahit --input metaspades \
  --parameter preset=asm5 --parameter "hits_kept=every alignment minimap2 reports"

cap validate "$bundle"
