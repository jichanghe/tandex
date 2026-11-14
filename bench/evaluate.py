#!/usr/bin/env python3
"""evaluate.py truth.tsv calls.bed4 [label] [max_true_period]  (calls: seq start end period, 0-based half-open)
Prints base-level recall/precision, locus recall, period accuracy, and recall by divergence/period class."""
import sys, bisect
from collections import defaultdict

truth = [l.split('\t') for l in open(sys.argv[1])]
truth = [(int(t[1]), int(t[2]), int(t[3]), float(t[4])) for t in truth]
if len(sys.argv) > 4: truth = [t for t in truth if t[2] <= int(sys.argv[4])]
calls = []
for l in open(sys.argv[2]):
    if l.startswith('#') or not l.strip(): continue
    f = l.split('\t'); calls.append((int(f[1]), int(f[2]), int(float(f[3]))))
calls.sort()
label = sys.argv[3] if len(sys.argv) > 3 else sys.argv[2]

def union(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out
U = union([(s, e) for s, e, p in calls]); Us = [u[0] for u in U]
def cov(s, e):
    c = 0; i = max(0, bisect.bisect_right(Us, s) - 1)
    while i < len(U) and U[i][0] < e:
        c += max(0, min(e, U[i][1]) - max(s, U[i][0])); i += 1
    return c
starts = [c[0] for c in calls]; maxlen = max([c[1] - c[0] for c in calls] + [1])
def best_period_ok(s, e, p):
    i = bisect.bisect_left(starts, s - maxlen)
    while i < len(calls) and calls[i][0] < e:
        cs, ce, cp = calls[i]; ov = min(e, ce) - max(s, cs)
        if ov >= 0.5 * (e - s) and abs(cp - p) <= max(0, 0.05 * p): return True
        i += 1
    return False

tb = sum(e - s for s, e, p, d in truth); tp = sum(cov(s, e) for s, e, p, d in truth)
pb = sum(e - s for s, e in U)
found = per = 0
bydiv = defaultdict(lambda: [0, 0]); byper = defaultdict(lambda: [0, 0])
def pclass(p): return '1-6' if p <= 6 else '7-100' if p <= 100 else '101-2000' if p <= 2000 else '>2000'
for s, e, p, d in truth:
    ok = cov(s, e) >= 0.5 * (e - s); pok = ok and best_period_ok(s, e, p)
    found += ok; per += pok
    bydiv[d][0] += 1; bydiv[d][1] += pok; byper[pclass(p)][0] += 1; byper[pclass(p)][1] += pok
n = len(truth)
print(f"{label}\tcalls={len(calls)}\tbase_recall={tp/tb:.3f}\tbase_precision={(tp/pb if pb else 0):.3f}\t"
      f"locus_recall={found/n:.3f}\tlocus+period={per/n:.3f}\tcalls_per_found_locus={(len(calls)/found if found else 0):.2f}")
print("  locus+period by divergence: " + "  ".join(f"{int(d*100)}%:{v[1]/v[0]:.2f}" for d, v in sorted(bydiv.items())))
print("  locus+period by period:     " + "  ".join(f"{k}:{byper[k][1]/byper[k][0]:.2f}" for k in ['1-6', '7-100', '101-2000', '>2000'] if k in byper))
