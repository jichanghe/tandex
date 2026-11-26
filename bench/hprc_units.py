#!/usr/bin/env python3
"""hprc_units.py extract DIR ASSEMBLY.fa.gz   ->  DIR/arrays.fa
   hprc_units.py units   DIR CALLS.tsv[.gz]  ->  DIR/units.tsv
DIR holds one haplotype's HPRC annotation: cenSat.bed, gaps.bed, chromAlias.txt, sample.tsv.
extract: writes every annotated active HOR array ('active_hor', 500 kb or more, not at a contig
  end or next to a gap) with up to 200 kb of flank, one record each, named
  seq|start|end|chr|family|left_flank.
units: for each array, the higher-order repeat unit (in monomers of about 171 bp) by two methods:
  tandex  - among call periods >= 300 bp covering at least half the array (periods within 1%
            merged), the shortest whose coverage is at least 0.9 of the best coverage and whose
            identity is within two points of the best (before 0.8.2: within one point, no coverage term);
  self-comparison (independent of tandex) - in 12 windows of 20 kb, for every unit of k monomers
            (k = 2..38) the identity between the sequence and itself k monomers further on, measured
            in 200-bp pieces each at its own best shift; the smallest k within 0.01 of the best;
            majority over the windows. 'confident' means 8 or more of 12 windows."""
import sys, gzip, collections, re
FLANK = 200000
def opn(f): return gzip.open(f, 'rt') if f.endswith('.gz') else open(f)
def extract(d, asm):
    alias = {}
    for l in open(d + '/chromAlias.txt'):
        x = l.split()
        if l[0] != '#' and len(x) > 1: alias[x[0]] = re.sub(r'_.*', '', x[1])
    ext = {}; arr = collections.defaultdict(list)
    for l in open(d + '/cenSat.bed'):
        x = l.rstrip('\n').split('\t')
        if len(x) < 4 or not x[1].isdigit(): continue
        s, e = int(x[1]), int(x[2]); lo, hi = ext.get(x[0], (s, e)); ext[x[0]] = (min(lo, s), max(hi, e))
        if x[3].startswith('active_hor') and e - s >= 500000: arr[x[0]].append((s, e, x[3][11:-1]))
    gaps = collections.defaultdict(list)
    for l in open(d + '/gaps.bed'):
        x = l.split()
        if len(x) > 2: gaps[x[0]].append((int(x[1]), int(x[2])))
    keep = {}
    for seq, ivs in arr.items():
        v = [(s, e, f) for s, e, f in ivs if s - ext[seq][0] > 50000 and ext[seq][1] - e > 50000 and not any(a < e + 1000 and b > s - 1000 for a, b in gaps[seq])]
        if v: keep[seq] = v
    n = 0
    with open(d + '/arrays.fa.part', 'w') as o:
        cur = None; name = None
        def flush():
            nonlocal n
            if cur is None: return
            s_all = b''.join(cur)
            for s, e, f in keep[name]:
                off = min(FLANK, s)
                o.write(f'>{name}|{s}|{e}|{alias.get(name, "?")}|{f}|{off}\n'); o.write(s_all[s - off:e + FLANK].decode() + '\n'); n += 1
        for l in gzip.open(asm, 'rb'):
            if l[:1] == b'>':
                flush(); name = l[1:].split()[0].decode(); cur = [] if name in keep else None
            elif cur is not None: cur.append(l.strip())
        flush()
    import os; os.replace(d + '/arrays.fa.part', d + '/arrays.fa')
    print(f'{n} arrays written')
def auto_unit(a):
    import numpy as np
    L = len(a); votes = collections.Counter()
    for w in range(12):
        st = max(0, int((w + 0.5) * L / 12) - 10000); seg = a[st:st + 20000]
        if len(seg) < 14000: continue
        n = (len(seg) - 6700) // 200; idk = {}
        for k in range(2, 39):
            best = np.zeros(n)
            for dd in range(int(169.5 * k) - 6, int(171.5 * k) + 7):
                best = np.maximum(best, (seg[:n * 200] == seg[dd:dd + n * 200]).reshape(n, 200).sum(1))
            idk[k] = best.mean() / 200
        top = max(idk.values()); votes[min(k for k in idk if idk[k] >= top - 0.01)] += 1
    return votes
def units(d, callsf):
    import numpy as np
    hap, sample, pop, sex = (open(d + '/sample.tsv').read().split() + ['-'] * 4)[:4]
    calls = collections.defaultdict(list)
    for l in opn(callsf):
        if l.startswith('#seq'): continue
        x = l.split('\t', 7)
        if len(x) > 6 and int(x[3]) >= 300: calls[x[0]].append((int(x[1]), int(x[2]), int(x[3]), float(x[6])))
    name = None; out = []
    def one(name, seq):
        sq, s, e, c, fam, off = name.split('|'); s = int(s); e = int(e); off = int(off); L = e - s
        per = collections.Counter(); idn = collections.Counter()
        for a, b, p, i in calls[name]:
            ov = min(b, off + L) - max(a, off)
            if ov > 0:
                k = next((q for q in per if abs(p - q) <= 0.01 * q), p); per[k] += ov; idn[k] += i * ov
        big = {k: (min(1, per[k] / L), idn[k] / per[k]) for k in per if per[k] >= 0.5 * L}
        p = 0; ident = 0.0
        if big:
            top = max(v[1] for v in big.values()); cmax = max(v[0] for v in big.values())
            ok = [k for k in big if big[k][1] >= top - 2.0 and big[k][0] >= 0.9 * cmax]
            p = min(ok) if ok else min(k for k in big if big[k][1] >= top - 2.0); ident = big[p][1]
        v = auto_unit(np.frombuffer(seq.upper().encode(), dtype=np.uint8)[off:off + L]).most_common(2)
        a, cnt = v[0] if v else (0, 0); second = f'{v[1][0]}:{v[1][1]}' if len(v) > 1 else '-'
        t = round(p / 171)
        out.append([hap, sample, pop, sex, c, fam, sq, s, e, L, p, t, '%.1f' % ident, a, cnt, second, 'yes' if t == a else 'no', 'yes' if cnt >= 8 else 'no'])
    for l in open(d + '/arrays.fa'):
        if l[0] == '>': name = l[1:].strip()
        else: one(name, l.strip())
    with open(d + '/units.tsv', 'w') as o:
        o.write('hap\tsample\tpop\tsex\tchr\tannotated_family\tseq\tstart\tend\tlength\ttandex_period\ttandex_n_mer\ttandex_identity\tself_n_mer\tself_windows\tself_second\tagree\tconfident\n')
        for r in out: o.write('\t'.join(map(str, r)) + '\n')
    print(f'{len(out)} arrays; agree {sum(r[16] == "yes" for r in out)}; confident and agree {sum(r[16] == "yes" and r[17] == "yes" for r in out)}; confident and differ {sum(r[16] == "no" and r[17] == "yes" for r in out)}')
if __name__ == '__main__':
    if len(sys.argv) > 3 and sys.argv[1] == 'extract': extract(sys.argv[2].rstrip('/'), sys.argv[3])
    elif len(sys.argv) > 3 and sys.argv[1] == 'units': units(sys.argv[2].rstrip('/'), sys.argv[3])
    else: sys.exit(__doc__)
