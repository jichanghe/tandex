#!/usr/bin/env bash
# run_apes.sh - download the six T2T ape genomes, run tandex on each, and write small summaries.
# Run in Ubuntu (WSL):   bash /mnt/c/chm13/run_apes.sh
# Safe to re-run: finished downloads and finished tandex runs are skipped.
# Needs: ~/tandex/tandex (version 0.7.0), curl, python3. About 25 GB of free disk on C:.
#
# Assemblies (Yoo et al. 2025, Nature, "Complete sequencing of ape genomes"), GenBank primary assemblies.
set -uo pipefail

TANDEX=${TANDEX:-$HOME/tandex/tandex}
OUT=/mnt/c/chm13/apes
RES=/mnt/c/chm13/results/apes
THREADS=${THREADS:-$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))}
mkdir -p "$OUT" "$RES"
LOG=$RES/run_apes.log
exec > >(tee -a "$LOG") 2>&1

echo "=== run_apes.sh started $(date) ==="
echo "cores: $(nproc)   workers used: $THREADS"
free -g | head -2
df -h /mnt/c | tail -1
[ -x "$TANDEX" ] || { echo "ERROR: $TANDEX not found. Build it first: cd ~/tandex && make"; exit 1; }
echo "tandex version: $($TANDEX -v)"

GENOMES="
chimpanzee:GCA_028858775.2
bonobo:GCA_029289425.2
gorilla:GCA_029281585.2
sumatran_orangutan:GCA_028885655.2
bornean_orangutan:GCA_028885625.2
siamang:GCA_028878055.2
"

for g in $GENOMES; do
  name=${g%%:*}; acc=${g##*:}
  fa=$OUT/$name.fa
  tsv=$RES/$name.tandex.tsv.gz
  echo; echo "--- $name ($acc) ---"

  if [ ! -s "$fa" ]; then
    zip=$OUT/$name.zip
    echo "downloading $acc ..."
    curl -L --fail --retry 5 --retry-delay 10 -C - -o "$zip" \
      "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" \
      || { echo "DOWNLOAD FAILED for $name ($acc); skipping"; continue; }
    python3 - "$zip" "$fa" <<'PY' || { echo "UNPACK FAILED for $name; skipping"; rm -f "$fa"; continue; }
import sys, zipfile, shutil
z = zipfile.ZipFile(sys.argv[1])
fna = [n for n in z.namelist() if n.endswith('.fna')]
assert fna, 'no .fna file in the download: ' + ', '.join(z.namelist()[:10])
with z.open(fna[0]) as src, open(sys.argv[2] + '.part', 'wb') as dst:
    shutil.copyfileobj(src, dst, 1 << 24)
PY
    mv "$fa.part" "$fa" && rm -f "$zip"
  fi
  ls -la "$fa"

  if [ ! -s "$tsv" ]; then
    echo "running tandex -p 50000 -t $THREADS ..."
    start=$(date +%s)
    "$TANDEX" -p 50000 -t "$THREADS" "$fa" > "$RES/$name.tandex.tsv.part" \
      || { echo "TANDEX FAILED for $name"; rm -f "$RES/$name.tandex.tsv.part"; continue; }
    end=$(date +%s)
    echo "TIME $name $((end - start)) s with $THREADS workers"
    gzip -f "$RES/$name.tandex.tsv.part" && mv "$RES/$name.tandex.tsv.part.gz" "$tsv"
  fi

  if [ ! -s "$RES/$name.summary.txt" ]; then
    python3 - "$fa" "$tsv" "$name" "$RES" <<'PY'
import sys, gzip, collections
fa, tsv, name, res = sys.argv[1:5]
size = {}; cur = None
with open(fa, 'rb') as f:
    for l in f:
        if l[:1] == b'>': cur = l[1:].split()[0].decode(); size[cur] = 0
        else: size[cur] += len(l) - 1
G = sum(size.values())
classes = [('1-6', 1, 6), ('7-100', 7, 100), ('101-1000', 101, 1000), ('1001-50000', 1001, 10**9)]
iv = {c[0]: collections.defaultdict(list) for c in classes}; allv = collections.defaultdict(list)
byp = collections.Counter(); ncalls = 0; arrays = []
for l in gzip.open(tsv, 'rt'):
    if l[0] == '#': continue
    x = l.rstrip('\n').split('\t'); s, e, p = int(x[1]), int(x[2]), int(x[3]); ncalls += 1
    allv[x[0]].append((s, e)); byp[p] += e - s
    for c, a, b in classes:
        if a <= p <= b: iv[c][x[0]].append((s, e))
    if e - s >= 10000:
        arrays.append((e - s, x[0], s, e, p, x[4], x[6], x[12] if len(x) > 12 else '', x[10] if p <= 400 else ''))
def union(d):
    t = 0
    for v in d.values():
        v.sort(); cs, ce = v[0]
        for s, e in v[1:]:
            if s > ce: t += ce - cs; cs, ce = s, e
            elif e > ce: ce = e
        t += ce - cs
    return t
with open(f'{res}/{name}.summary.txt', 'w') as o:
    o.write(f'genome\t{name}\nsequences\t{len(size)}\ngenome_bp\t{G}\ncalls\t{ncalls}\n')
    o.write(f'bp_in_repeats\t{union(allv)}\n')
    for c, a, b in classes: o.write(f'bp_period_{c}\t{union(iv[c]) if iv[c] else 0}\n')
    o.write('# top 60 periods by bases called (overlapping calls counted separately)\n')
    for p, n in byp.most_common(60): o.write(f'period\t{p}\t{n}\n')
arrays.sort(reverse=True)
with open(f'{res}/{name}.arrays.tsv', 'w') as o:
    o.write('#length\tseq\tstart\tend\tperiod\tcopies\tidentity\ttype\tconsensus_if_period_le_400\n')
    for a in arrays[:5000]: o.write('\t'.join(map(str, a)) + '\n')
with open(f'{res}/{name}.chromsizes.tsv', 'w') as o:
    for k, v in size.items(): o.write(f'{k}\t{v}\n')
print(f'summary written: {name}: genome {G} bp, {ncalls} calls')
PY
  fi
done

echo; echo "=== finished $(date) ==="
echo "Results are in C:\\chm13\\results\\apes  (tell Claude when this line appears)"
