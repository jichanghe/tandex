#!/usr/bin/env python3
"""hor_check.py arrays.fa calls.tsv
Compares the higher-order repeat (HOR) unit tandex reports for each array with an independent
estimate from autocorrelation. arrays.fa holds one record per array, named chr|name|start|end,
with up to 200 kb of flank on each side (made from censat 'hor' intervals of 500 kb or more).
tandex unit: among call periods >= 300 bp that cover at least half the array (periods within 1%
merged), the shortest whose identity is within one point of the best.
Autocorrelation unit: in six 20-kb windows, the shortest shift (300 to 6500 bp) whose ungapped
self-identity is within 0.02 of the best; the majority over the windows. Units are in monomers
(period / 171)."""
import sys, collections, numpy as np
seqs = {}; name = None
for l in open(sys.argv[1]):
    if l[0] == '>': name = l[1:].strip(); seqs[name] = []
    else: seqs[name].append(l.strip())
calls = collections.defaultdict(list)
for l in open(sys.argv[2]):
    if l[0] == '#': continue
    x = l.split('\t', 13); calls[x[0]].append((int(x[1]), int(x[2]), int(x[3]), float(x[6])))
key = lambda c: (int(c.split('|')[0][3:]) if c.split('|')[0][3:].isdigit() else 99)
agree = 0
print('chr\tarray\tMb\tautocorrelation_n-mer\twindows\ttandex_n-mer\ttandex_identity\tagree')
for nm in sorted(seqs, key=key):
    c, fam, s, e = nm.split('|'); s = int(s); e = int(e); L = e - s; off = min(200000, s)
    a = np.frombuffer(''.join(seqs[nm]).upper().encode(), dtype=np.uint8)[off:off + L]
    votes = collections.Counter()
    for w in range(6):
        st = int((w + 0.5) * L / 6) - 10000; seg = a[st:st + 20000]; best = (0, 0)
        for d in range(300, 6500):
            v = (seg[:-d] == seg[d:]).mean()
            if v > best[0] + 0.02: best = (v, d)
        votes[round(best[1] / 171)] += 1
    auto, cnt = votes.most_common(1)[0]
    per = collections.Counter(); idn = collections.Counter()
    for x, y, p, i in calls[nm]:
        if p < 300: continue
        ov = min(y, off + L) - max(x, off)
        if ov > 0:
            k = next((q for q in per if abs(p - q) <= 0.01 * q), p); per[k] += ov; idn[k] += i * ov
    big = {k: idn[k] / per[k] for k in per if per[k] >= 0.5 * L}
    t = tid = 0
    if big:
        top = max(big.values()); p = min(k for k in big if big[k] >= top - 1.0); t = round(p / 171); tid = big[p]
    agree += t == auto
    print(c, fam, '%.2f' % (L / 1e6), auto, f'{cnt}/6', t, '%.1f' % tid, 'yes' if t == auto else 'no', sep='\t')
print(f'# agree {agree} of {len(seqs)}')
