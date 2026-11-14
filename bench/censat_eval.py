#!/usr/bin/env python3
"""censat_eval.py chroms name=calls.tsv ...  : coverage of the T2T-CHM13 censat v2.1 satellite annotation.
chroms: comma list or 'all'. Per class: % of annotated bases covered by any call; for alpha satellite
(hor, dhor, mon) also % covered by calls whose period is a multiple of 171 (+-3%)."""
import sys
from collections import defaultdict
chroms = None if sys.argv[1] == 'all' else set(sys.argv[1].split(','))
CL = ['hor', 'dhor', 'mon', 'hsat1A', 'hsat1B', 'hsat2', 'hsat3', 'bsat', 'gsat', 'censat', 'rDNA']
ann = defaultdict(lambda: defaultdict(list))
for l in open('censat.bed'):
    x = l.split('\t')
    if x[0] == 'chrom' or (chroms and x[0] not in chroms): continue
    n = x[3].split('_')[0].split('(')[0]
    if n in CL: ann[n][x[0]].append((int(x[1]), int(x[2])))
def union(iv):
    o = []
    for s, e in sorted(iv):
        if o and s <= o[-1][1]: o[-1][1] = max(o[-1][1], e)
        else: o.append([s, e])
    return o
def inter(a, b):
    i = j = t = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0]); e = min(a[i][1], b[j][1])
        if e > s: t += e - s
        if a[i][1] < b[j][1]: i += 1
        else: j += 1
    return t
tot = {c: sum(e - s for v in ann[c].values() for s, e in v) for c in CL}
print("class\tMb\t" + "\t".join(a.split('=')[0] for a in sys.argv[2:]))
rows = {c: [] for c in CL}; alpha = []; allrow = []
for a in sys.argv[2:]:
    name, f = a.split('=')
    calls = defaultdict(list); a171 = defaultdict(list)
    for l in open(f):
        if l[0] == '#' or l.startswith('SeqID'): continue
        x = l.split('\t')
        try: c, s, e, p = x[0], int(x[1]), int(x[2]), int(float(x[3]))
        except ValueError: continue
        if chroms and c not in chroms: continue
        calls[c].append((s, e))
        k = round(p / 171)
        if k >= 1 and abs(p - k * 171) <= 0.03 * p + 1: a171[c].append((s, e))
    U = {c: union(v) for c, v in calls.items()}; A = {c: union(v) for c, v in a171.items()}
    cov_all = 0
    for cl in CL:
        cov = sum(inter(union(v), U.get(c, [])) for c, v in ann[cl].items()); cov_all += cov
        rows[cl].append(100 * cov / tot[cl] if tot[cl] else float('nan'))
    allrow.append(100 * cov_all / sum(tot.values()))
    ta = sum(tot[c] for c in ('hor', 'dhor', 'mon'))
    alpha.append(100 * sum(inter(union(v), A.get(c, [])) for cl in ('hor', 'dhor', 'mon') for c, v in ann[cl].items()) / ta if ta else float('nan'))
for cl in CL:
    if tot[cl]: print(f"{cl}\t{tot[cl]/1e6:.2f}\t" + "\t".join(f"{v:.1f}%" for v in rows[cl]))
print(f"all satellite classes\t{sum(tot.values())/1e6:.2f}\t" + "\t".join(f"{v:.1f}%" for v in allrow))
print(f"alpha satellite, period n x 171\t{sum(tot[c] for c in ('hor','dhor','mon'))/1e6:.2f}\t" + "\t".join(f"{v:.1f}%" for v in alpha))
