#!/usr/bin/env bash
# Run ULTRA with a long period limit. On the 7 GB benchmark machine ULTRA could only be run up to
# -p 500 (4.3 GB); longer limits need more memory. NOT TESTED at these settings by the tandex author.
#
# usage: run_ultra_long.sh <genome.fa> <out_prefix> [max_period=2000] [threads=8]
#   ULTRA=/path/to/ultra may be set; otherwise 'ultra' must be on PATH.
set -euo pipefail
FA=$1; OUT=$2; P=${3:-2000}; T=${4:-8}
ULTRA=${ULTRA:-ultra}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

echo "== memory ULTRA says it needs for -p $P -t $T"
"$ULTRA" -p "$P" -t "$T" --mem "$FA" || true
avail=$(awk '/MemAvailable/{printf "%d", $2/1024}' /proc/meminfo 2>/dev/null || echo unknown)
echo "== memory available: $avail MB. Press Ctrl-C within 10 s to stop."
sleep 10

/usr/bin/time -v "$ULTRA" -p "$P" -t "$T" --tsv --hs --disable_summary -o "$OUT.ultra.tsv" "$FA" \
    > /dev/null 2> "$OUT.ultra.time"
grep -E "Elapsed|Maximum resident" "$OUT.ultra.time" || true

# columns: SeqID Start End Period ... ; keep the first four for the evaluation scripts
awk -F'\t' -v OFS='\t' 'NR>1 && $2 ~ /^[0-9]+$/ {print $1, $2, $3, $4}' "$OUT.ultra.tsv" > "$OUT.ultra.calls.tsv"
echo "calls: $OUT.ultra.calls.tsv ($(wc -l < "$OUT.ultra.calls.tsv"))"
cat <<MSG
Evaluate, for example:
  python3 $HERE/evaluate.py sim2.truth.tsv $OUT.ultra.calls.tsv ULTRA-p$P
  python3 $HERE/censat_eval.py chr21 ULTRA-p$P=$OUT.ultra.calls.tsv
MSG
