#!/usr/bin/env python3
"""Simulate a random genome with implanted tandem repeats. Writes FASTA + truth TSV.
usage: simulate.py out_prefix seed n_repeats [min_spacer max_spacer]"""
import random, sys

prefix, seed, nrep = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
smin, smax = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (1000, 4000)
rnd = random.Random(seed)
PERIODS = [1, 2, 3, 4, 5, 6, 8, 12, 17, 25, 40, 64, 100, 171, 300, 500, 1000, 1800, 3000, 6000]
DIVS = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25]

def rseq(n): return ''.join(rnd.choice('ACGT') for _ in range(n))

def primitive_unit(p):
    while True:
        u = rseq(p)
        if p == 1 or all(u != u[q:] + u[:q] for q in range(1, p)) and (p > 3 or len(set(u)) > 1):
            return u

def mutate(s, div):
    """per-base divergence: 80% substitutions, 10% insertions, 10% deletions"""
    out = []
    for ch in s:
        if rnd.random() < div:
            r = rnd.random()
            if r < 0.8: out.append(rnd.choice([c for c in 'ACGT' if c != ch]))
            elif r < 0.9: out.append(ch); out.append(rnd.choice('ACGT'))
        else: out.append(ch)
    return ''.join(out)

parts, truth, pos = [], [], 0
for i in range(nrep):
    sp = rseq(rnd.randint(smin, smax)); parts.append(sp); pos += len(sp)
    p = rnd.choice(PERIODS); div = rnd.choice(DIVS)
    if p == 1 and div > 0.1: div = 0.1
    minlen = max(40, int(2.5 * p)) if p > 6 else 30
    copies = max(minlen / p, rnd.choice([3, 5, 10, 20, 50]) if p < 500 else rnd.choice([2.5, 4, 8, 15]))
    unit = primitive_unit(p)
    arr = (unit * (int(copies) + 2))[:int(copies * p)]
    arr = mutate(arr, div)
    parts.append(arr); truth.append((pos, pos + len(arr), p, div, copies)); pos += len(arr)
parts.append(rseq(rnd.randint(smin, smax)))
g = ''.join(parts)
with open(prefix + '.fa', 'w') as f:
    f.write('>sim\n')
    for i in range(0, len(g), 80): f.write(g[i:i + 80] + '\n')
with open(prefix + '.truth.tsv', 'w') as f:
    for t in truth: f.write('sim\t%d\t%d\t%d\t%.2f\t%.1f\n' % t)
print(len(g), 'bp,', len(truth), 'repeats', file=sys.stderr)
