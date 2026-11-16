#!/usr/bin/env python3
"""Bases in tandem repeats by period family, from a tandex TSV (plain or .gz).
usage: classify_periods.py NAME calls.tsv[.gz] GENOME_BP
Families are assigned from the period alone. Nested calls of different families cover the same
bases, so the family columns can add up to more than repeat_bp."""
import sys, gzip, collections
name, tsv, G = sys.argv[1], sys.argv[2], int(sys.argv[3])
def nm(p, q, tol):
    m = round(p / q)
    return m >= 1 and abs(p - m * q) <= min(tol * p, 8)
def cls(p):
    if p >= 165 and nm(p, 171, 0.012): return 'alpha_171n'
    if p in (32, 64, 96, 128, 31, 33): return 'p32'
    if p <= 6: return 'micro_1-6' if p != 5 else 'p5'
    if p % 5 == 0 and p <= 200: return 'p5'
    if p <= 100: return 'mini_7-100'
    if p <= 1000: return 'p101-1000'
    return 'p1001+'
iv = collections.defaultdict(lambda: collections.defaultdict(list)); allv = collections.defaultdict(list)
n = 0
op = gzip.open if tsv.endswith('.gz') else open
for l in op(tsv, 'rt'):
    if l[0] == '#': continue
    x = l.split('\t', 5); s, e, p = int(x[1]), int(x[2]), int(x[3]); n += 1
    iv[cls(p)][x[0]].append((s, e)); allv[x[0]].append((s, e))
def union(d, minlen=0):
    t = 0; k = 0
    for v in d.values():
        v.sort(); cs, ce = v[0]
        for s, e in v[1:]:
            if s > ce:
                if ce - cs >= minlen: t += ce - cs; k += 1
                cs, ce = s, e
            elif e > ce: ce = e
        if ce - cs >= minlen: t += ce - cs; k += 1
    return t, k
tot, _ = union(allv); m1, k1 = union(allv, 1000000); m100, k100 = union(allv, 100000)
out = [name, G, n, tot, f'{100*tot/G:.2f}', k100, m100, k1, m1]
for c in ['micro_1-6', 'p5', 'mini_7-100', 'p32', 'alpha_171n', 'p101-1000', 'p1001+']:
    out.append(union(iv[c])[0] if iv[c] else 0)
print('\t'.join(map(str, out)))
