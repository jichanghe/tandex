#!/usr/bin/env bash
# run_hprc_hor.sh - re-run the 20 pilot haplotypes with tandex -H (higher-order repeat mode).
# Run in Ubuntu (WSL):   bash /mnt/c/chm13/run_hprc_hor.sh
# Uses the assemblies already in C:\chm13\hprc; nothing is downloaded. Needs tandex 0.8.0.
# About 25 to 30 minutes per haplotype (about 9 hours in all). Safe to stop and re-run.
set -uo pipefail
TANDEX=${TANDEX:-$HOME/tandex/tandex}
ASM=/mnt/c/chm13/hprc; RES=/mnt/c/chm13/results/hprc; THREADS=${THREADS:-12}
exec > >(tee -a "$RES/run_hprc_hor.log") 2>&1
echo "=== run_hprc_hor.sh started $(date) ==="
V=$($TANDEX -v); echo "tandex $V, $THREADS workers"
case "$V" in 0.[8-9]*|[1-9]*) ;; *) echo "ERROR: need tandex 0.8.0 or later (update ~/tandex first)"; exit 1;; esac
n=0; ok=0
for d in "$RES"/*/; do
  name=$(basename "$d"); f=$ASM/$name.fa.gz
  [ -s "$d/sample.tsv" ] || continue
  n=$((n+1)); echo; echo "--- [$n] $name ---"
  if [ -s "$d/tandexH.tsv.gz" ]; then echo "already done"; ok=$((ok+1)); continue; fi
  [ -s "$f" ] || { echo "assembly not found: $f"; continue; }
  s=$(date +%s)
  "$TANDEX" -H -p 50000 -t "$THREADS" "$f" > "$d/tandexH.tsv.part" || { echo "TANDEX FAILED: $name"; rm -f "$d/tandexH.tsv.part"; continue; }
  echo "TIME $name $(( $(date +%s) - s )) s; calls: $(grep -vc '^#' "$d/tandexH.tsv.part")"
  gzip -f "$d/tandexH.tsv.part" && mv "$d/tandexH.tsv.part.gz" "$d/tandexH.tsv.gz" && ok=$((ok+1))
done
echo; echo "=== finished $(date): $ok of $n haplotypes done ==="
echo "Tell Claude when this line appears (or earlier: a few finished haplotypes are enough to start)"
