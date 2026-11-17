#!/usr/bin/env bash
# run_fix.sh - two leftovers.  Run in Ubuntu (WSL):  bash /mnt/c/chm13/run_fix.sh
#  1  download the published ape satellite annotation (12 small files) to C:\chm13\ape_annot
#  2  retry TRASH on CHM13 chr21 with its R packages taken from conda, then score all tools
set -uo pipefail
RES=/mnt/c/chm13/results/next; WORK=$HOME/tandex_work; BENCH=$HOME/tandex/bench
A=/mnt/c/chm13/ape_annot; mkdir -p "$A" "$RES"
exec > >(tee -a "$RES/run_fix.log") 2>&1
echo "=== run_fix.sh started $(date) ==="

echo "##### 1: annotation files #####"
base=https://genomeark.s3.amazonaws.com/species
while read -r dir id date ver; do
  curl -sL --fail --retry 3 -o "$A/$id.bb"  "$base/$dir/$id/assembly_curated/repeats/${id}_${ver}_CenSat_v2.0.bb" || echo "FAILED: $id annotation"
  curl -sL --fail --retry 3 -o "$A/$id.fai" "$base/$dir/$id/assembly_curated/$id.pri.cur.$date.fasta.gz.fai"   || echo "FAILED: $id fai"
done <<'LIST'
Pan_troglodytes mPanTro3 20231122 v2.0
Pan_paniscus mPanPan1 20231122 v2.0
Gorilla_gorilla mGorGor1 20231122 v2.0
Pongo_abelii mPonAbe1 20231205 v2.0
Pongo_pygmaeus mPonPyg2 20231122 v2.0
Symphalangus_syndactylus mSymSyn1 20240514 v2.1
LIST
ls -la "$A"
echo ">>> part 1 done: you can tell Claude now; part 2 may take hours"

echo "##### 2: TRASH #####"
cd "$WORK/rivals" || { echo "no $WORK/rivals: run run_next.sh first"; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate tandex_rivals || { echo "conda environment tandex_rivals missing"; exit 1; }
PK="Biostrings stringr stringdist seqinr doParallel circlize"
conda install -y -c conda-forge -c bioconda bioconductor-biostrings r-stringr r-stringdist r-seqinr r-doparallel r-circlize > conda_r.log 2>&1 \
  || { echo "conda could not install the R packages:"; tail -15 conda_r.log; }
ok=1; for p in $PK; do Rscript -e "suppressPackageStartupMessages(library($p))" > /dev/null 2>&1 || { echo "R package missing: $p"; ok=0; }; done
if [ "$ok" = 1 ]; then
  rm -rf trash_out; mkdir -p trash_out
  s=$(date +%s)
  timeout 10h TRASH/TRASH_run.sh "$PWD/chr21.fa" --o "$PWD/trash_out" --par 12 --win 5000 --m 5000 --def > trash_run.log 2>&1
  echo "TIME TRASH --par 12 --win 5000: $(( $(date +%s) - s )) s (exit code $?)"
  f=trash_out/all.repeats.from.chr21.fa.csv
  if [ -s "$f" ]; then python3 "$BENCH/trash2tsv.py" "$f" > trash.calls.tsv && echo "TRASH arrays: $(wc -l < trash.calls.tsv)"
  else echo "TRASH produced no output:"; tail -12 trash_run.log; ls trash_out | head; fi
else
  echo "TRASH not run"
fi
args="tandex=tandex.tsv ULTRA-p2000=ultra2000.calls.tsv"; [ -s trash.calls.tsv ] && args="$args TRASH=trash.calls.tsv"
python3 "$BENCH/censat_eval.py" chr21 $args | tee "$RES/chr21_rivals_censat.txt"
[ -s trash.calls.tsv ] && gzip -c trash.calls.tsv > "$RES/chr21.trash.calls.tsv.gz"
echo "=== finished $(date) ==="
