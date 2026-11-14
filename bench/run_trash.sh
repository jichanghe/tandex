#!/usr/bin/env bash
# Run TRASH (Wlodzimierz et al. 2023) on a FASTA file and convert its output for the evaluation
# scripts. NOT TESTED by the tandex author: TRASH needs R >= 4.4.2 and its R packages, which were
# not available on the benchmark machine. Written from the TRASH README and its example output.
#
# usage: run_trash.sh <genome.fa> <out_dir> [cores]
#   <genome.fa> must end in .fa, .fna or .fasta; <out_dir> is created and must be writable.
#
# One-time setup:
#   git clone https://github.com/vlothec/TRASH && cd TRASH && chmod +x TRASH_install.sh && ./TRASH_install.sh
#   export TRASH_DIR=/path/to/TRASH
set -euo pipefail
FA=$(readlink -f "$1"); OUT=$(mkdir -p "$2" && readlink -f "$2"); CORES=${3:-1}
: "${TRASH_DIR:?set TRASH_DIR to the TRASH checkout}"
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
NAME=$(basename "$FA")

# --win/--m: TRASH cannot report repeats longer than its window (default 1000 bp).
# 5000 lets it see alpha-satellite higher-order units; raise it for longer periods (slower).
/usr/bin/time -v "$TRASH_DIR/TRASH_run.sh" "$FA" --o "$OUT" --par "$CORES" --win 5000 --m 5000 \
    > "$OUT/trash.log" 2> "$OUT/trash.time"
grep -E "Elapsed|Maximum resident" "$OUT/trash.time" || true

python3 "$HERE/trash2tsv.py" "$OUT/all.repeats.from.$NAME.csv" > "$OUT/trash.calls.tsv"
echo "calls: $OUT/trash.calls.tsv ($(wc -l < "$OUT/trash.calls.tsv") arrays)"
cat <<MSG
Evaluate with the same scripts used for the other tools, for example:
  python3 $HERE/evaluate.py sim2.truth.tsv $OUT/trash.calls.tsv TRASH          # simulation
  python3 $HERE/cen_eval.py TRASH=$OUT/trash.calls.tsv                         # Arabidopsis mask (run in the folder with cen180.mask.bed)
  python3 $HERE/censat_eval.py chr21 TRASH=$OUT/trash.calls.tsv                # CHM13 (run in the folder with censat.bed)
MSG
