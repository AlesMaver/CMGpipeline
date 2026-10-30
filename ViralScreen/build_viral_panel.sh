#!/usr/bin/env bash
# Build the viral panel reference used by ViralScreen_CRAM.wdl from viral_panel_accessions.tsv (accession<TAB>name).
# Usage: build_viral_panel.sh <outdir>   (needs curl, samtools and bwa)
set -euo pipefail
OUT=$1; mkdir -p "$OUT"
: > "$OUT/viral_panel.fa"
while read acc name; do
  curl -sf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta" \
    | awk -v n=$name 'NR==1{print ">"n"|"substr($1,2); next}{print}' >> "$OUT/viral_panel.fa"
  sleep 0.4
done < "$(dirname "$0")/viral_panel_accessions.tsv"
samtools faidx "$OUT/viral_panel.fa"
bwa index "$OUT/viral_panel.fa"
