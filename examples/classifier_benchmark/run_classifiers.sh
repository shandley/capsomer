#!/usr/bin/env bash
# Run five viral classifiers on one FASTA file and record how long each takes.
# Usage: run_classifiers.sh <fasta> <output folder> [tool ...]
# TOOLS is the folder holding one environment per tool (env, env_vs2x, env_ct3, env_phabox,
# env_jaeger) and their databases under db/. See ../viroforge_skin/run_tools.sh for how the
# environments were made.
set -uo pipefail
fasta="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
out="$2"; shift 2
tools="${TOOLS:?set TOOLS to the folder of tool environments}"
threads="${THREADS:-10}"
mkdir -p "$out"; out="$(cd "$out" && pwd)"
wanted="${*:-genomad virsorter2 cenotetaker3 phabox2 jaeger}"

timed() { # name, command...
  local name="$1"; shift
  local start=$SECONDS
  ( "$@" ) > "$out/$name.log" 2>&1
  local status=$?
  printf '%s\t%s\t%s\n' "$name" "$((SECONDS - start))" "$status" >> "$out/timing.tsv"
  echo "$name: exit $status after $((SECONDS - start)) s"
}

for tool in $wanted; do
  case "$tool" in
    genomad)
      PATH="$tools/env/bin:$PATH" timed genomad genomad end-to-end --cleanup --threads "$threads" \
        "$fasta" "$out/genomad" "$tools/db/genomad_db" ;;
    virsorter2)
      PATH="$tools/env_vs2x/bin:$PATH" timed virsorter2 virsorter run -w "$out/virsorter2" -i "$fasta" \
        -d "$tools/db/vs2" -j "$threads" --use-conda-off --min-length 1000 all ;;
    cenotetaker3)
      PATH="$tools/env_ct3/bin:$PATH" timed cenotetaker3 cenotetaker3 -c "$fasta" -r ct3 -p F \
        -t "$threads" --cenote-dbs "$tools/db/ct3" -wd "$out/cenotetaker3" --genbank False ;;
    phabox2)
      PATH="$tools/env_phabox/bin:$PATH" timed phabox2 phabox2 --task phamer \
        --dbdir "$tools/db/phabox_db_v2_2" --outpth "$out/phabox2" --contigs "$fasta" --threads "$threads" ;;
    jaeger)
      PATH="$tools/env_jaeger/bin:$PATH" timed jaeger jaeger predict -i "$fasta" -o "$out/jaeger" \
        -m "${JAEGER_MODEL:-jaeger_448931d3_1.3M_fragment}" --workers "$threads" ;;
    *) echo "unknown tool $tool" ;;
  esac
done
