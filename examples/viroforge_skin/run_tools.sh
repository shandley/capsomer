#!/usr/bin/env bash
# Run the contig-side tools on a ViroForge dataset. Usage: run_tools.sh <ViroForge dataset folder>
#
# Needs fastp, megahit, spades.py, genomad, checkv, vclust, mmseqs, minimap2, samtools and coverm
# on the PATH, and the geNomad and CheckV databases under tools/db. One way to get them:
#   micromamba create -p tools/env -c conda-forge -c bioconda megahit spades genomad checkv \
#     minimap2 samtools coverm vclust fastp
#   genomad download-database tools/db && checkv download_database tools/db
# Versions used for the example: fastp 1.4.0, MEGAHIT 1.2.9, SPAdes 4.3.0, geNomad 1.13.0 (database
# 1.9), CheckV 1.1.1 (database 1.5), vclust 1.3.1, MMseqs2 19.ab57c, minimap2 2.31, samtools 1.25,
# CoverM 0.8.0.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
dataset="${1:?give the ViroForge dataset folder}"
work="$here/work/contig"
db="$here/tools/db"
threads="${THREADS:-12}"
mkdir -p "$work"/{qc,megahit,metaspades,votu,mmseqs,map}

r1="$(ls "$dataset"/fastq/*_R1.fastq)"
fastp -i "$r1" -I "${r1%_R1.fastq}_R2.fastq" -o "$work/qc/R1.fastq.gz" -O "$work/qc/R2.fastq.gz" \
  --detect_adapter_for_pe --json "$work/qc/fastp.json" --html "$work/qc/fastp.html" --thread 8

# MEGAHIT 1.2.9 crashes with more than one thread on Apple Silicon, hence -t 1.
megahit -1 "$work/qc/R1.fastq.gz" -2 "$work/qc/R2.fastq.gz" -o "$work/assembly" --min-contig-len 1000 -t 1
cp "$work/assembly/final.contigs.fa" "$work/megahit/contigs.fasta"

spades.py --meta -1 "$work/qc/R1.fastq.gz" -2 "$work/qc/R2.fastq.gz" -o "$work/spades" -t "$threads" -m 64
awk 'BEGIN{RS=">";ORS=""} NR>1{n=index($0,"\n"); h=substr($0,1,n-1); s=substr($0,n+1); gsub("\n","",s);
  if(length(s)>=1000) print ">"h"\n"s"\n"}' "$work/spades/contigs.fasta" > "$work/metaspades/contigs.fasta"

for assembler in megahit metaspades; do
  genomad end-to-end --cleanup --threads "$threads" "$work/$assembler/contigs.fasta" \
    "$work/$assembler/genomad" "$db/genomad_db"
  checkv end_to_end "$work/$assembler/genomad/contigs_summary/contigs_virus.fna" \
    "$work/$assembler/checkv" -t "$threads" -d "$db/checkv-db-v1.5"
done

# Species-level clusters of every sequence geNomad called viral, from both assemblies.
cat "$work"/{megahit,metaspades}/genomad/contigs_summary/contigs_virus.fna > "$work/votu/viral_sequences.fna"
vclust prefilter -i "$work/votu/viral_sequences.fna" -o "$work/votu/prefilter.txt"
vclust align -i "$work/votu/viral_sequences.fna" -o "$work/votu/ani.tsv" --filter "$work/votu/prefilter.txt"
vclust cluster -i "$work/votu/ani.tsv" -o "$work/votu/clusters.tsv" --ids "$work/votu/ani.ids.tsv" \
  --algorithm leiden --metric ani --ani 0.95 --qcov 0.85 -r
# A second, looser clustering with a second tool.
mmseqs easy-cluster "$work/votu/viral_sequences.fna" "$work/mmseqs/genus" "$work/mmseqs/tmp" \
  --min-seq-id 0.7 -c 0.5 --cov-mode 1 --threads 8

# Map reads to the species representatives, and compute coverage.
tail -n +2 "$work/votu/clusters.tsv" | cut -f2 | sort -u > "$work/votu/representatives.txt"
samtools faidx "$work/votu/viral_sequences.fna" -r "$work/votu/representatives.txt" -n 0 \
  > "$work/votu/representatives.fna"
minimap2 -ax sr -t "$threads" "$work/votu/representatives.fna" "$work/qc/R1.fastq.gz" "$work/qc/R2.fastq.gz" \
  | samtools sort -@ 4 -o "$work/map/reads_to_representatives.bam" -
samtools index "$work/map/reads_to_representatives.bam"
coverm contig -b "$work/map/reads_to_representatives.bam" -m count mean covered_fraction covered_bases \
  -o "$work/map/coverm.tsv"

# Align every contig to the simulated genomes, to find each contig's true source.
cat "$work"/{megahit,metaspades}/contigs.fasta > "$work/map/all_contigs.fasta"
minimap2 -x asm5 -c -t "$threads" "$(ls "$dataset"/fasta/*.fasta)" "$work/map/all_contigs.fasta" \
  > "$work/map/contigs_to_truth.paf"
