#!/usr/bin/env python3
"""hprc_verify.py cen_hor_H.tsv ASSEMBLY_DIR out.tsv
Independent check of the higher-order repeat unit of each centromere array listed in the table
(made by hprc_cen_hor.py), by self-comparison of the assembly sequence (see unit()); the unit is
the majority over 12 windows of 20 kb, in monomers of about 171 bp. Needs numpy. Reads each <haplotype>.fa.gz in ASSEMBLY_DIR once."""
import sys, gzip, csv, collections, time
import numpy as np
table, asm, out = sys.argv[1:4]
rows = [r for r in csv.DictReader(open(table), delimiter='\t') if r['internal'] == '1']
byhap = collections.defaultdict(list)
for r in rows: byhap[r['hap']].append(r)
def unit(a):
    """votes over 12 windows of 20 kb. In each window, for every unit of k monomers (k = 2..38) the
    identity between the sequence and itself k monomers further on is measured in 200-bp pieces,
    each piece at its own best shift near k x 171 bp (so that small insertions and deletions do not
    count against long units); the unit is the smallest k within 0.01 of the best identity."""
    L = len(a); votes = collections.Counter()
    for w in range(12):
        st = max(0, int((w + 0.5) * L / 12) - 10000); seg = a[st:st + 20000]
        if len(seg) < 14000: continue
        n = (len(seg) - 6700) // 200
        idk = {}
        for k in range(2, 39):
            best = np.zeros(n)
            for d in range(int(169.5 * k) - 6, int(171.5 * k) + 7):
                eq = (seg[:n * 200] == seg[d:d + n * 200]).reshape(n, 200).sum(1)
                best = np.maximum(best, eq)
            idk[k] = best.mean() / 200
        top = max(idk.values())
        votes[min(k for k in idk if idk[k] >= top - 0.01)] += 1
    return votes
with open(out, 'w') as o:
    o.write('hap\tsample\tpop\tchr\tannotated_family\tlength\ttandex_n_mer\ttandex_identity\tauto_n_mer\tauto_windows\tauto_second\tagree\n')
    for hap in sorted(byhap):
        t0 = time.time(); want = collections.defaultdict(list)
        for r in byhap[hap]: want[r['seq']].append(r)
        seqs = {}; cur = None
        try: f = gzip.open(f'{asm}/{hap}.fa.gz', 'rb')
        except Exception as e: print('cannot open', hap, e, flush=True); continue
        for l in f:
            if l[:1] == b'>':
                name = l[1:].split()[0].decode(); cur = [] if name in want else None
                if cur is not None: seqs[name] = cur
            elif cur is not None: cur.append(l.strip().upper())
        f.close(); n = ok = 0
        for name, rs in want.items():
            if name not in seqs: print('sequence not found:', hap, name, flush=True); continue
            s = np.frombuffer(b''.join(seqs[name]), dtype=np.uint8)
            for r in rs:
                v = unit(s[int(r['start']):int(r['end'])]).most_common(2)
                if not v: continue
                a, c = v[0]; second = f'{v[1][0]}:{v[1][1]}' if len(v) > 1 else ''
                agree = str(a) == r['n_mer']; n += 1; ok += agree
                o.write('\t'.join([hap, r['sample'], r['pop'], r['chr'], r['annotated_family'], r['length'], r['n_mer'], r['hor_identity'], str(a), f'{c}/12', second, 'yes' if agree else 'no']) + '\n')
            o.flush()
        print(f'{hap}: {ok} of {n} arrays agree ({time.time() - t0:.0f} s)', flush=True)
print('written', out)
