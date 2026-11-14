#!/usr/bin/env python3
"""Recall of each tool against the Broad tandem repeat catalog v1.0 on GRCh38 chr21.
usage: catalog_eval.py name=calls.tsv:conscol ...   (conscol = 0-based column of the consensus)
A catalog locus is 'found' if >=50% of it is covered.
  exact      : a call covering >=50% has period == catalog motif length
  compatible : exact, or one motif is the other repeated (lengths differ by an integer factor)
               with >=80% identity in the best rotation, e.g. CCTTCCTTCCTTCGTT vs CCTT"""
import bisect, sys
def compatible(m, c):
    if not c or not m: return False
    a, b = (m, c) if len(m) >= len(c) else (c, m)
    if len(a) % len(b): return False
    if len(a) > 400: return len(a) == len(b)
    best = 0
    for r in range(len(b)):
        u = (b[r:] + b[:r]) * (len(a) // len(b))
        best = max(best, sum(x == y for x, y in zip(a, u)))
    return best >= 0.8 * len(a)
def load(f, cc):
    out = []
    for l in open(f):
        if l[0] == '#' or l.startswith('SeqID'): continue
        x = l.rstrip('\n').split('\t')
        try: out.append((int(x[1]), int(x[2]), int(float(x[3])), x[cc].upper() if cc < len(x) else ''))
        except ValueError: pass
    return sorted(out)
cat = [(int(x[1]), int(x[2]), x[3].strip().upper()) for x in (l.split('\t') for l in open('../truth/cat21.bed'))]
bins = [(25, 49), (50, 99), (100, 10**9)]
print("tool\tcalls\t" + "\t".join(f"{a}-{b if b < 10**9 else 'max'} bp: found / exact / compatible" for a, b in bins) + "\t9-24 bp found")
for arg in sys.argv[1:]:
    name, rest = arg.split('='); f, cc = rest.split(':')
    C = load(f, int(cc)); st = [c[0] for c in C]; mx = max(c[1] - c[0] for c in C)
    res = {b: [0, 0, 0, 0] for b in bins}; small = [0, 0]
    for s, e, m in cat:
        i = bisect.bisect_left(st, s - mx); cov = []; ex = co = False
        while i < len(C) and C[i][0] < e:
            cs, ce, cp, cons = C[i]
            if ce > s:
                cov.append((max(s, cs), min(e, ce)))
                if min(e, ce) - max(s, cs) >= 0.5 * (e - s):
                    if cp == len(m): ex = True
                    if cp == len(m) or compatible(m, cons): co = True
            i += 1
        c = 0; last = s
        for a, b in sorted(cov):
            a = max(a, last)
            if b > a: c += b - a; last = b
        found = c >= 0.5 * (e - s)
        if e - s < 25: small[0] += 1; small[1] += found
        for b in bins:
            if b[0] <= e - s <= b[1]:
                res[b][0] += 1; res[b][1] += found; res[b][2] += found and ex; res[b][3] += found and co
    print(f"{name}\t{len(C)}\t" + "\t".join(f"{100*res[b][1]/res[b][0]:.1f}% / {100*res[b][2]/res[b][0]:.1f}% / {100*res[b][3]/res[b][0]:.1f}% (n={res[b][0]})" for b in bins) + f"\t{100*small[1]/small[0]:.1f}%")
