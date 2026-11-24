#!/usr/bin/env bash
# run_recheck.sh - re-check centromere units with tandex 0.8.2. Run in Ubuntu (WSL):
#     bash /mnt/c/chm13/run_recheck.sh
# Part 1: chr12 and chr15 arrays of all 461 haplotypes (about 2 hours).
# Part 2: every array of 40 randomly chosen haplotypes, to check nothing else got worse (about 2.5 hours).
# Uses the arrays already cut out; nothing is downloaded. Safe to stop and re-run.
set -uo pipefail
TANDEX=${TANDEX:-$HOME/tandex/tandex}; THREADS=${THREADS:-12}
RES=/mnt/c/chm13/results/hprc_all; B=/mnt/c/chm13
exec > >(tee -a "$RES/run_recheck.log") 2>&1
echo "=== run_recheck.sh started $(date) ==="
V=$($TANDEX -v); echo "tandex $V"
case "$V" in 0.8.[2-9]*|0.9*|[1-9]*) ;; *) echo "ERROR: need tandex 0.8.2 (update ~/tandex first)"; exit 1;; esac
echo "##### part 1: chr12 and chr15, all haplotypes #####"
python3 "$B/hprc_recheck.py" "$TANDEX" "$THREADS" "$RES" chr12,chr15 all "$RES/recheck_chr12_15.tsv"
echo "##### part 2: all arrays, 40 haplotypes #####"
python3 "$B/hprc_recheck.py" "$TANDEX" "$THREADS" "$RES" all "$B/recheck_haps40.txt" "$RES/recheck_sample40.tsv"
echo "=== finished $(date) ===  (tell Claude)"
