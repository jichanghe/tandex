#!/usr/bin/env bash
# run_hprc_all.sh - centromere repeat units for all HPRC release 2 haplotypes (462).
# Run in Ubuntu (WSL):   bash /mnt/c/chm13/run_hprc_all.sh            # everything
#                        bash /mnt/c/chm13/run_hprc_all.sh pilot      # only the 20 pilot haplotypes (a check, about 40 min)
# Per haplotype: download the assembly (about 0.9 GB), cut out the active centromere arrays,
# delete the download, run tandex -H on the arrays only, and work out each array's repeat unit
# by tandex and by an independent self-comparison. About 5 minutes each; about 1.5 days in all.
# Safe to stop (Ctrl+C) and re-run: finished haplotypes are skipped.
# Needs: tandex 0.8.1 in ~/tandex, python3 with numpy, files hprc_all.tsv and hprc_units.py in C:\chm13.
set -uo pipefail
MODE=${1:-all}
TANDEX=${TANDEX:-$HOME/tandex/tandex}; THREADS=${THREADS:-12}
BASE=/mnt/c/chm13; RES=$BASE/results/hprc_all; ASM=$BASE/hprc; TMP=$HOME/tandex_work/hprc_tmp
mkdir -p "$RES" "$TMP"
exec > >(tee -a "$RES/run_hprc_all.log") 2>&1
echo "=== run_hprc_all.sh ($MODE) started $(date) ==="
V=$($TANDEX -v); echo "tandex $V"
case "$V" in 0.8.[1-9]*|0.9*|[1-9]*) ;; *) echo "ERROR: need tandex 0.8.1 or later (update ~/tandex first)"; exit 1;; esac
python3 -c 'import numpy' || { echo "ERROR: numpy missing: pip install numpy"; exit 1; }
[ -s "$BASE/hprc_all.tsv" ] && [ -s "$BASE/hprc_units.py" ] || { echo "ERROR: hprc_all.tsv or hprc_units.py not found in $BASE"; exit 1; }
get() { [ -s "$2" ] || curl -sL --fail --retry 5 --retry-delay 10 -o "$2" "$1" || { echo "  download failed: $1"; rm -f "$2"; return 1; }; }
total=$(wc -l < "$BASE/hprc_all.tsv"); n=0; ok=0
while IFS=$'\t' read -r name sample pop sex phasing fa censat alias gaps; do
  [ -n "$name" ] || continue
  n=$((n+1)); d=$RES/$name
  if [ "$MODE" = pilot ] && [ ! -s "$ASM/$name.fa.gz" ]; then continue; fi
  if [ -s "$d/units.tsv" ]; then ok=$((ok+1)); continue; fi
  echo; echo "--- [$n/$total] $name ($pop) $(date +%H:%M) ---"
  mkdir -p "$d"; printf "%s\t%s\t%s\t%s\n" "$name" "$sample" "$pop" "$sex" > "$d/sample.tsv"
  get "$censat" "$d/cenSat.bed" && get "$alias" "$d/chromAlias.txt" && get "$gaps" "$d/gaps.bed" || { echo "annotation missing, skipped"; continue; }
  if [ ! -s "$d/arrays.fa" ]; then
    f=$ASM/$name.fa.gz; mine=0
    if [ ! -s "$f" ]; then
      f=$TMP/$name.fa.gz; mine=1
      curl -sL --fail --retry 5 --retry-delay 10 -C - -o "$f.part" "$fa" && mv "$f.part" "$f" || { echo "ASSEMBLY DOWNLOAD FAILED"; continue; }
    fi
    python3 "$BASE/hprc_units.py" extract "$d" "$f" || { echo "EXTRACT FAILED"; [ $mine = 1 ] && rm -f "$f"; continue; }
    [ $mine = 1 ] && rm -f "$f"
  fi
  s=$(date +%s)
  "$TANDEX" -H -p 50000 -t "$THREADS" "$d/arrays.fa" > "$d/arrays.tandexH.tsv" || { echo "TANDEX FAILED"; continue; }
  echo "tandex: $(( $(date +%s) - s )) s"
  python3 "$BASE/hprc_units.py" units "$d" "$d/arrays.tandexH.tsv" || { echo "UNITS FAILED"; rm -f "$d/units.tsv"; continue; }
  gzip -f "$d/arrays.tandexH.tsv"; gzip -f "$d/arrays.fa"
  ok=$((ok+1))
done < "$BASE/hprc_all.tsv"
head -1 "$(ls "$RES"/*/units.tsv | head -1)" > "$RES/all_units.tsv"
for f in "$RES"/*/units.tsv; do tail -n +2 "$f"; done >> "$RES/all_units.tsv"
echo; echo "=== finished $(date): $ok haplotypes have results; table: C:\\chm13\\results\\hprc_all\\all_units.tsv ==="
