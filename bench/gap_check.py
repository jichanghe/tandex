#!/usr/bin/env python3
"""gap_check.py calls.tsv chrom.fa chrom class lagfrom lagto
Are the annotated satellite bases that tandex leaves uncovered satellite-like? Tool-independent test:
share of bases equal to the base one monomer downstream (best lag), in 500-bp windows."""
import sys, numpy as np
f, fa, c, cl, l0, l1 = sys.argv[1:7]; lags = range(int(l0), int(l1) + 1)
seq = np.frombuffer(''.join(l.strip() for l in open(fa) if l[0] != '>').upper().encode(), dtype=np.uint8)
calls = sorted((int(x[1]), int(x[2])) for x in (l.split('\t', 4) for l in open(f) if l[0] != '#') if x[0] == c)
def auto(s, e): return max((seq[s:e] == seq[s + g:e + g]).mean() for g in lags if e + g <= len(seq))
cov = np.zeros(len(seq), dtype=bool)
for a, b in calls: cov[a:b] = True
tot = miss = sat = 0; cal = []
for l in open('censat.bed'):
    x = l.split('\t')
    if x[0] != c or x[3].split('_')[0].split('(')[0] != cl: continue
    s, e = int(x[1]), int(x[2]); tot += e - s
    for w in range(s, e - 499, 500):
        m = int((~cov[w:w + 500]).sum())
        if m >= 400:
            miss += 500
            if auto(w, w + 500) >= 0.40: sat += 500
        elif m == 0: cal.append(auto(w, w + 500))
print(f"{c} {cl}: annotated {tot/1e6:.2f} Mb; fully uncovered 500-bp windows {miss/1e3:.0f} kb, of which satellite-like {sat/1e3:.0f} kb "
      f"({100*sat/max(1,miss):.0f}%); covered windows passing the same test {100*np.mean(np.array(cal)>=0.40):.0f}% (n={len(cal)})")
