#!/usr/bin/env python3
"""hprc_cen_hor.py CALLS_NAME DIR...   (CALLS_NAME: tandex.tsv.gz or tandexH.tsv.gz)
Per annotated active HOR array (HPRC cenSat 'active_hor', >= 500 kb, not at a contig end or gap)
of each haplotype directory: its length and the HOR unit tandex reports, taken as the shortest
call period >= 300 bp that covers at least half the array and whose identity is within one point
of the best such period (periods within 1% merged). Writes cen_hor.tsv; prints a summary."""
import sys, gzip, collections, re, statistics as st
cname = sys.argv[1]; rows = []
for d in sys.argv[2:]:
    d = d.rstrip('/')
    try: hap, sample, pop, sex = open(d + '/sample.tsv').read().split()
    except Exception: continue
    alias = {}
    for l in open(d + '/chromAlias.txt'):
        x = l.split()
        if l[0] != '#' and len(x) > 1: alias[x[0]] = re.sub(r'_.*', '', x[1])
    ext = {}; arr = collections.defaultdict(list)
    for l in open(d + '/cenSat.bed'):
        x = l.rstrip('\n').split('\t')
        if len(x) < 4 or not x[1].isdigit(): continue
        s, e = int(x[1]), int(x[2])
        lo, hi = ext.get(x[0], (s, e)); ext[x[0]] = (min(lo, s), max(hi, e))
        if x[3].startswith('active_hor') and e - s >= 500000: arr[x[0]].append((s, e, x[3][11:-1]))
    gaps = collections.defaultdict(list)
    for l in open(d + '/gaps.bed'):
        x = l.split()
        if len(x) > 2: gaps[x[0]].append((int(x[1]), int(x[2])))
    calls = collections.defaultdict(list)
    for l in gzip.open(d + '/' + cname, 'rt'):
        if l[0] == '#': continue
        x = l.split('\t', 7)
        if x[0] in arr and int(x[3]) >= 300: calls[x[0]].append((int(x[1]), int(x[2]), int(x[3]), float(x[6])))
    for seq, ivs in arr.items():
        for s, e, fam in ivs:
            internal = s - ext[seq][0] > 50000 and ext[seq][1] - e > 50000 and not any(a < e + 1000 and b > s - 1000 for a, b in gaps[seq])
            per = collections.Counter(); idn = collections.Counter()
            for a, b, p, i in calls[seq]:
                ov = min(b, e) - max(a, s)
                if ov <= 0: continue
                k = next((q for q in per if abs(p - q) <= 0.01 * q), p); per[k] += ov; idn[k] += i * ov
            L = e - s
            big = {k: idn[k] / per[k] for k in per if per[k] >= 0.5 * L}
            p = 0; ident = 0.0
            if big:
                top = max(big.values()); p = min(k for k in big if big[k] >= top - 1.0); ident = big[p]
            rows.append([hap, sample, pop, alias.get(seq, '?'), fam, seq, s, e, L, int(internal), p, round(p / 171), '%.1f' % ident])
with open('cen_hor.tsv', 'w') as o:
    o.write('hap\tsample\tpop\tchr\tannotated_family\tseq\tstart\tend\tlength\tinternal\thor_period\tn_mer\thor_identity\n')
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
by = collections.defaultdict(list)
for r in rows:
    if r[9]: by[r[3]].append(r)
print('chr\tarrays\tMb_min\tMb_median\tMb_max\twith_unit\tn-mers (n: arrays)\tidentity_median')
key = lambda c: (int(c[3:]) if c[3:].isdigit() else 99, c)
for c in sorted(by, key=key):
    v = by[c]; L = sorted(r[8] / 1e6 for r in v); ok = [r for r in v if r[10]]
    nm = collections.Counter(r[11] for r in ok)
    print(c, len(v), '%.2f' % L[0], '%.2f' % st.median(L), '%.2f' % L[-1], f'{len(ok)}/{len(v)}',
          ' '.join(f'{n}:{k}' for n, k in sorted(nm.items())), '%.1f' % (st.median(float(r[12]) for r in ok) if ok else 0), sep='\t')
