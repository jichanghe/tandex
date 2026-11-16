#!/usr/bin/env bash
# run_next.sh - three follow-up jobs for the tandex paper. Run in Ubuntu (WSL):
#     bash /mnt/c/chm13/run_next.sh            # parts D, E and C (A and B are already done)
#     bash /mnt/c/chm13/run_next.sh A          # only one part: A, B, C, D or E
#
#  A  list the annotation files of the six ape genomes on GenomeArk (1 minute)
#  B  measure how tandex scales with 1, 2, 4, 8 and 12 workers (about 40 minutes)
#  C  run ULTRA with long periods and TRASH on T2T-CHM13 chr21 and score them (hours; best effort)
#  D  score the tandex ape calls against the published CenSat v2.0 satellite annotation (10 minutes)
#  E  record the processor layout and the CPU time used with 1, 6 and 12 workers (15 minutes)
#
# Needs tandex 0.7.2 in ~/tandex. Each part is independent: if one fails the others still run.
# Everything is written to C:\chm13\results\next. Large temporary files go to ~/tandex_work.
set -uo pipefail
el() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.1f", b - a }'; }

PARTS=${1:-DEC}
TANDEX=${TANDEX:-$HOME/tandex/tandex}
BENCH=$HOME/tandex/bench
RES=/mnt/c/chm13/results/next
WORK=$HOME/tandex_work
mkdir -p "$RES" "$WORK"
exec > >(tee -a "$RES/run_next.log") 2>&1

echo "=== run_next.sh ($PARTS) started $(date) ==="
echo "cores: $(nproc)"; free -g | head -2
[ -x "$TANDEX" ] || { echo "ERROR: $TANDEX not found"; exit 1; }
V=$($TANDEX -v); echo "tandex version: $V"
case "$V" in 0.7.[2-9]*|0.[8-9]*|[1-9]*) ;; *) echo "ERROR: need tandex 0.7.2 or later (update ~/tandex first)"; exit 1;; esac

# ----------------------------------------------------------------------------------------------
part_A() {
  echo; echo "##### A: GenomeArk file listing #####"
  python3 - "$RES/genomeark_listing.tsv" <<'PY'
import sys, re, urllib.request, urllib.parse
species = ['Pan_troglodytes/mPanTro3', 'Pan_paniscus/mPanPan1', 'Gorilla_gorilla/mGorGor1',
           'Pongo_abelii/mPonAbe1', 'Pongo_pygmaeus/mPonPyg2', 'Symphalangus_syndactylus/mSymSyn1']
out = open(sys.argv[1], 'w'); total = 0
for sp in species:
    token = None; n = 0
    for page in range(40):
        q = {'list-type': '2', 'prefix': f'species/{sp}/assembly_curated/', 'max-keys': '1000'}
        if token: q['continuation-token'] = token
        url = 'https://genomeark.s3.amazonaws.com/?' + urllib.parse.urlencode(q)
        try:
            xml = urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception as e:
            print(f'  {sp}: listing failed: {e}'); break
        for key, size in re.findall(r'<Key>([^<]+)</Key>.*?<Size>(\d+)</Size>', xml, re.S):
            out.write(f'{key}\t{size}\n'); n += 1
        m = re.search(r'<NextContinuationToken>([^<]+)</NextContinuationToken>', xml)
        if not m: break
        token = m.group(1)
    print(f'  {sp}: {n} files'); total += n
print(f'listing written: {total} files')
PY
  grep -i -E "censat|satellite|\.fai|chrom.?sizes" "$RES/genomeark_listing.tsv" | head -80
}

# ----------------------------------------------------------------------------------------------
part_B() {
  echo; echo "##### B: scaling with workers #####"
  local src=/mnt/c/chm13/apes/chimpanzee.fa sub=$WORK/chimp_sub.fa
  [ -s "$src" ] || { echo "B skipped: $src not found"; return; }
  if [ ! -s "$sub" ]; then
    echo "making a subset of about 600 Mb on the Linux file system ..."
    awk 'BEGIN{n=0} /^>/{ if (n >= 600000000) exit } {print; if ($0 !~ /^>/) n += length($0)}' "$src" > "$sub.part" && mv "$sub.part" "$sub"
  fi
  ls -la "$sub"; grep -c '>' "$sub"
  : > "$RES/scaling.tsv"
  for t in 12 8 4 2 1; do
    [ "$t" -le "$(nproc)" ] || continue
    echo "tandex -p 50000 -t $t ..."
    s=$(date +%s.%N)
    TANDEX_PROFILE=1 "$TANDEX" -p 50000 -t "$t" "$sub" 2> "$RES/scaling.t$t.profile" | grep -vc '^#' > "$WORK/calls.t$t"
    e=$(date +%s.%N)
    printf "workers\t%s\tseconds\t%.1f\tcalls\t%s\t%s\n" "$t" "$(el "$s" "$e")" "$(cat "$WORK/calls.t$t")" "$(grep 'profile total' "$RES/scaling.t$t.profile")" | tee -a "$RES/scaling.tsv"
  done
  echo "reading the same file from the Windows drive with 12 workers ..."
  cp "$sub" /mnt/c/chm13/apes/chimp_sub.fa
  s=$(date +%s.%N); "$TANDEX" -p 50000 -t 12 /mnt/c/chm13/apes/chimp_sub.fa > /dev/null; e=$(date +%s.%N)
  printf "workers\t12 (file on C:)\tseconds\t%.1f\n" "$(el "$s" "$e")" | tee -a "$RES/scaling.tsv"
  rm -f /mnt/c/chm13/apes/chimp_sub.fa
}

# ----------------------------------------------------------------------------------------------
part_C() {
  echo; echo "##### C: ULTRA (long periods) and TRASH on CHM13 chr21 #####"
  local gz=/mnt/c/chm13/chm13v2.0.fa.gz bed=/mnt/c/chm13/chm13v2.0_censat_v2.1.bed
  [ -s "$gz" ] && [ -s "$bed" ] || { echo "C skipped: CHM13 genome or censat bed not found in C:\\chm13"; return; }
  mkdir -p "$WORK/rivals" && cd "$WORK/rivals" || return
  cp "$bed" censat.bed
  if [ ! -s chr21.fa ]; then
    echo "extracting chr21 ..."
    # awk stops reading once chr21 is complete; zcat then ends with a broken pipe, which is expected
    zcat "$gz" 2> /dev/null | awk '/^>/{ if (f) exit; f = ($1 == ">chr21") } f' > chr21.part || true
    [ -s chr21.part ] && mv chr21.part chr21.fa
  fi
  ls -la chr21.fa || { echo "C stopped: chr21 could not be extracted"; return; }

  echo "--- tandex on this machine"
  s=$(date +%s.%N); "$TANDEX" -p 50000 -t 12 chr21.fa > tandex.tsv; e=$(date +%s.%N)
  printf "TIME tandex -p 50000 -t 12: %.1f s\n" "$(el "$s" "$e")"

  # tools are built inside a separate conda environment so that nothing else on this machine changes
  if command -v conda > /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    if ! conda env list | grep -q '^tandex_rivals '; then
      echo "creating conda environment tandex_rivals (cmake, compilers, R) ..."
      conda create -y -n tandex_rivals -c conda-forge cmake make cxx-compiler c-compiler "r-base>=4.4.2" > conda_create.log 2>&1 \
        || echo "conda environment could not be created (see $WORK/rivals/conda_create.log)"
    fi
    conda activate tandex_rivals 2> /dev/null || echo "could not activate tandex_rivals; using system tools"
  else
    echo "conda not found; using system cmake / R if present"
  fi

  echo "--- ULTRA -p 2000"
  (
    if [ ! -x ULTRA/ultra ]; then
      rm -rf ULTRA; git clone --depth 1 https://github.com/TravisWheelerLab/ULTRA > ultra_build.log 2>&1 \
        && (cd ULTRA && cmake . >> ../ultra_build.log 2>&1 && make -j 8 >> ../ultra_build.log 2>&1)
    fi
    [ -x ULTRA/ultra ] || { echo "ULTRA build FAILED (see $WORK/rivals/ultra_build.log)"; tail -5 ultra_build.log; exit 1; }
    echo "memory ULTRA says it needs:"; ULTRA/ultra -p 2000 -t 4 --mem chr21.fa 2>&1 | tail -3
    s=$(date +%s.%N)
    timeout 6h ULTRA/ultra -p 2000 -t 4 --tsv --hs --disable_summary -o ultra2000.tsv chr21.fa > /dev/null 2> ultra_run.log
    rc=$?; e=$(date +%s.%N)
    printf "TIME ULTRA -p 2000 -t 4: %.1f s (exit code %s)\n" "$(el "$s" "$e")" "$rc"
    [ "$rc" = 0 ] || { tail -5 ultra_run.log; exit 1; }
    awk -F'\t' -v OFS='\t' 'NR>1 && $2 ~ /^[0-9]+$/ {print $1, $2, $3, $4}' ultra2000.tsv > ultra2000.calls.tsv
    echo "ULTRA calls: $(wc -l < ultra2000.calls.tsv)"
  ) || echo "ULTRA part did not complete"

  echo "--- TRASH"
  (
    command -v Rscript > /dev/null || { echo "TRASH skipped: R (Rscript) not available"; exit 1; }
    Rscript --version 2>&1 | head -1
    if [ ! -d TRASH ]; then git clone --depth 1 https://github.com/vlothec/TRASH > trash_build.log 2>&1; fi
    DEF=""
    if [ ! -e TRASH/.installed ]; then
      # TRASH's own installer could not compile Biostrings here, so take the R packages from conda
      # and tell TRASH to use the default R library (--def)
      conda install -y -n tandex_rivals -c conda-forge -c bioconda bioconductor-biostrings r-stringr r-stringdist \
          r-seqinr r-doparallel r-circlize >> trash_build.log 2>&1 || echo "conda could not install the R packages"
      (cd TRASH && chmod +x TRASH_install.sh && ./TRASH_install.sh --def) >> trash_build.log 2>&1 && touch TRASH/.installed \
        || { echo "TRASH install FAILED (see $WORK/rivals/trash_build.log)"; tail -8 trash_build.log; exit 1; }
    fi
    DEF="--def"
    mkdir -p trash_out
    s=$(date +%s.%N)
    timeout 8h TRASH/TRASH_run.sh "$PWD/chr21.fa" --o "$PWD/trash_out" --par 12 --win 5000 --m 5000 $DEF > trash_run.log 2>&1
    rc=$?; e=$(date +%s.%N)
    printf "TIME TRASH --par 12 --win 5000: %.1f s (exit code %s)\n" "$(el "$s" "$e")" "$rc"
    f=trash_out/all.repeats.from.chr21.fa.csv
    [ -s "$f" ] || { echo "TRASH produced no output"; tail -8 trash_run.log; ls trash_out | head; exit 1; }
    python3 "$BENCH/trash2tsv.py" "$f" > trash.calls.tsv && echo "TRASH arrays: $(wc -l < trash.calls.tsv)"
  ) || echo "TRASH part did not complete"

  echo "--- satellite coverage on chr21 (censat v2.1)"
  args="tandex=tandex.tsv"
  [ -s ultra2000.calls.tsv ] && args="$args ULTRA-p2000=ultra2000.calls.tsv"
  [ -s trash.calls.tsv ] && args="$args TRASH=trash.calls.tsv"
  python3 "$BENCH/censat_eval.py" chr21 $args | tee "$RES/chr21_rivals_censat.txt"
  for f in ultra2000.calls.tsv trash.calls.tsv; do [ -s "$f" ] && gzip -c "$f" > "$RES/chr21.$f.gz"; done
  for f in ultra_build.log ultra_run.log trash_build.log trash_run.log conda_create.log; do [ -s "$f" ] && tail -40 "$f" > "$RES/$f.tail.txt"; done
}

# ----------------------------------------------------------------------------------------------
part_D() {
  echo; echo "##### D: tandex calls against the published ape satellite annotation (CenSat v2.0) #####"
  local A=$WORK/ape_annot; mkdir -p "$A"
  if [ ! -x "$WORK/venv/bin/python" ]; then python3 -m venv "$WORK/venv" > /dev/null 2>&1; fi
  "$WORK/venv/bin/python" -c 'import pyBigWig' 2> /dev/null || "$WORK/venv/bin/pip" install -q pyBigWig > "$A/pip.log" 2>&1
  "$WORK/venv/bin/python" -c 'import pyBigWig' 2> /dev/null || { echo "D stopped: pyBigWig could not be installed (see $A/pip.log)"; tail -3 "$A/pip.log"; return; }
  local base=https://genomeark.s3.amazonaws.com/species
  : > "$RES/apes_censat.txt"
  while read -r name dir id date ver; do
    [ -n "$name" ] || continue
    local calls=/mnt/c/chm13/results/apes/$name.tandex.tsv.gz sizes=/mnt/c/chm13/results/apes/$name.chromsizes.tsv
    [ -s "$calls" ] && [ -s "$sizes" ] || { echo "$name: calls not found, skipped"; continue; }
    [ -s "$A/$id.bb" ]  || curl -sL --fail --retry 3 -o "$A/$id.bb"  "$base/$dir/$id/assembly_curated/repeats/${id}_${ver}_CenSat_v2.0.bb" || { echo "$name: annotation download failed"; continue; }
    [ -s "$A/$id.fai" ] || curl -sL --fail --retry 3 -o "$A/$id.fai" "$base/$dir/$id/assembly_curated/$id.pri.cur.$date.fasta.gz.fai" || { echo "$name: .fai download failed"; continue; }
    "$WORK/venv/bin/python" "$BENCH/ape_censat_eval.py" "$name" "$A/$id.bb" "$A/$id.fai" "$sizes" "$calls" | tee -a "$RES/apes_censat.txt"
    echo | tee -a "$RES/apes_censat.txt"
  done <<'LIST'
chimpanzee Pan_troglodytes mPanTro3 20231122 v2.0
bonobo Pan_paniscus mPanPan1 20231122 v2.0
gorilla Gorilla_gorilla mGorGor1 20231122 v2.0
sumatran_orangutan Pongo_abelii mPonAbe1 20231205 v2.0
bornean_orangutan Pongo_pygmaeus mPonPyg2 20231122 v2.0
siamang Symphalangus_syndactylus mSymSyn1 20240514 v2.1
LIST
}

# ----------------------------------------------------------------------------------------------
part_E() {
  echo; echo "##### E: processor layout and CPU time with 12 workers #####"
  lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|MHz" | tee "$RES/cpu.txt"
  local sub=$WORK/chimp_sub.fa
  [ -s "$sub" ] || { echo "E skipped: run part B first"; return; }
  for t in 1 6 12; do
    echo "tandex -p 50000 -t $t (wall, user and system CPU seconds):"
    ( TIMEFORMAT="workers $t  wall %R  user %U  sys %S"; time "$TANDEX" -p 50000 -t "$t" "$sub" > /dev/null ) 2>&1 | tee -a "$RES/cpu.txt"
  done
}

case "$PARTS" in *A*) part_A;; esac
case "$PARTS" in *B*) part_B;; esac
case "$PARTS" in *D*) part_D;; esac
case "$PARTS" in *E*) part_E;; esac
case "$PARTS" in *C*) part_C;; esac
echo; echo "=== finished $(date) ==="
echo "Results are in C:\\chm13\\results\\next  (tell Claude when this line appears)"
