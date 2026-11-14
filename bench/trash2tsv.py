#!/usr/bin/env python3
"""trash2tsv.py all.repeats.from.<fasta>.csv > calls.tsv
TRASH lists one row per repeat monomer. Merge adjacent monomers of one region into array calls
(seq, start, end, period) so the evaluation scripts in this folder can read them.
start is converted from 1-based inclusive to 0-based. Period = median monomer width."""
import csv, statistics, sys

rows = []
with open(sys.argv[1], newline='') as f:
    for r in csv.DictReader(f):
        rows.append((r['seq.name'], r.get('region.name', ''), int(r['start']), int(r['end']), int(r['width'])))
rows.sort()
cur = None
def flush(c):
    if c: print(f"{c[0]}\t{c[2] - 1}\t{c[3]}\t{int(statistics.median(c[4]))}")
for seq, region, s, e, w in rows:
    if cur and cur[0] == seq and cur[1] == region and s - cur[3] <= 10:
        cur[3] = max(cur[3], e); cur[4].append(w)
    else:
        flush(cur); cur = [seq, region, s, e, [w]]
flush(cur)
