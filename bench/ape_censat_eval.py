#!/usr/bin/env python3
"""ape_censat_eval.py NAME annotation.(bb|bed) annotation_names.fai our.chromsizes.tsv calls.tsv[.gz]
Coverage of a T2T CenSat satellite annotation by tandem repeat calls, per annotation class.
The annotation and the calls may name the sequences differently (GenomeArk names against GenBank
accessions); sequences are paired by identical length. Classes 'ct' (centromeric transition, not
satellite) and gaps are left out. For each class the commonest call period (by bases) is shown."""
import sys, re, gzip, bisect, collections
name, annf, fai, sizes, callsf = sys.argv[1:6]
theirs = {}
for l in open(fai):
    x = l.split('\t'); theirs[x[0]] = int(x[1])
ours = {}
for l in open(sizes):
    x = l.split('\t'); ours[x[0]] = int(x[1])
bylen = collections.defaultdict(list)
for k, v in ours.items(): bylen[v].append(k)
m = {}; amb = 0
for k, v in theirs.items():
    if len(bylen.get(v, [])) == 1: m[k] = bylen[v][0]
    elif len(bylen.get(v, [])) > 1: amb += 1
mapped_bp = sum(ours[v] for v in m.values())
print(f'# {name}: {len(m)} of {len(ours)} of our sequences paired by length ({100 * mapped_bp / sum(ours.values()):.1f}% of bases); ambiguous lengths: {amb}')

def cls(n):
    if n.startswith('active_hor'): return 'active_hor'
    return re.split(r'[\(_]', n)[0]
ann = collections.defaultdict(lambda: collections.defaultdict(list))   # class -> seq -> intervals
flat = collections.defaultdict(list)                                    # seq -> (start, end, class)
def add(c, s, e, n):
    if c not in m: return
    k = cls(n)
    if k.lower() in ('ct', 'gap', ''): return
    ann[k][m[c]].append((s, e)); flat[m[c]].append((s, e, k))
if annf.endswith('.bed'):
    for l in open(annf):
        x = l.rstrip('\n').split('\t')
        if len(x) > 3 and x[1].isdigit(): add(x[0], int(x[1]), int(x[2]), x[3])
else:
    import pyBigWig
    bb = pyBigWig.open(annf)
    for c, ln in bb.chroms().items():
        for s, e, rest in (bb.entries(c, 0, ln) or []): add(c, s, e, rest.split('\t')[0])

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
for v in flat.values(): v.sort()
starts = {k: [a[0] for a in v] for k, v in flat.items()}
maxlen = {k: max(a[1] - a[0] for a in v) for k, v in flat.items()}
calls = collections.defaultdict(list); per = collections.defaultdict(collections.Counter)
op = gzip.open if callsf.endswith('.gz') else open
for l in op(callsf, 'rt'):
    if l[0] == '#': continue
    x = l.split('\t', 5); c = x[0]
    if c not in flat: continue
    s, e, p = int(x[1]), int(x[2]), int(x[3]); calls[c].append((s, e))
    v = flat[c]; i = bisect.bisect_left(starts[c], s - maxlen[c])
    while i < len(v) and v[i][0] < e:
        ov = min(e, v[i][1]) - max(s, v[i][0])
        if ov > 0: per[v[i][2]][p] += ov
        i += 1
cu = {c: union(v) for c, v in calls.items()}
print('class\tannotated_Mb\tcovered_%\tcommonest_call_periods(bp:share of call bases)')
allann = collections.defaultdict(list); rows = []
for k in ann:
    tot = cov = 0
    for c, iv in ann[k].items():
        u = union(iv); tot += sum(e - s for s, e in u); cov += inter(u, cu.get(c, [])); allann[c] += iv
    sp = sum(per[k].values()) or 1
    top = ' '.join(f'{p}:{100 * n / sp:.0f}%' for p, n in per[k].most_common(3))
    rows.append((tot, k, cov, top))
for tot, k, cov, top in sorted(rows, reverse=True):
    if tot >= 50000: print(f'{k}\t{tot / 1e6:.2f}\t{100 * cov / tot:.1f}\t{top}')
tot = cov = 0
for c, iv in allann.items():
    u = union(iv); tot += sum(e - s for s, e in u); cov += inter(u, cu.get(c, []))
print(f'ALL satellite classes\t{tot / 1e6:.2f}\t{100 * cov / tot:.1f}\t')
