#!/usr/bin/env python3
"""hprc_recheck.py TANDEX THREADS RESULTS_DIR CHRS|all HAPS|all OUT.tsv
Re-runs tandex -H on arrays already cut out by hprc_units.py (DIR/arrays.fa.gz) and works out
each array's unit again, reusing the independent self-comparison already in DIR/units.tsv (it
does not depend on tandex). CHRS: comma-separated chromosomes (for example chr12,chr15) or all.
HAPS: a file with one haplotype name per line, or all.
Unit rule (from version 0.8.2 of this analysis): among call periods >= 300 bp covering at least
half the array (periods within 1% merged), the shortest whose coverage is at least 0.9 of the best
coverage and whose identity is within two points of the best identity."""
import sys, os, gzip, csv, collections, subprocess, tempfile, time
tandex, threads, res, chrs, haps, out = sys.argv[1:7]
chrs = None if chrs == 'all' else set(chrs.split(','))
dirs = sorted(d for d in os.listdir(res) if os.path.isfile(f'{res}/{d}/units.tsv') and os.path.isfile(f'{res}/{d}/arrays.fa.gz'))
if haps != 'all':
    want = set(l.strip() for l in open(haps) if l.strip()); dirs = [d for d in dirs if d in want]
def unit(calls, off, L):
    per = collections.Counter(); idn = collections.Counter()
    for a, b, p, i in calls:
        if p < 300: continue
        ov = min(b, off + L) - max(a, off)
        if ov > 0:
            k = next((q for q in per if abs(p - q) <= 0.01 * q), p); per[k] += ov; idn[k] += i * ov
    big = {k: (min(1, per[k] / L), idn[k] / per[k]) for k in per if per[k] >= 0.5 * L}
    if not big: return 0, 0.0
    top = max(v[1] for v in big.values()); cmax = max(v[0] for v in big.values())
    ok = [k for k in big if big[k][1] >= top - 2.0 and big[k][0] >= 0.9 * cmax]
    p = min(ok) if ok else min(k for k in big if big[k][1] >= top - 2.0)  # no period is both near-best: drop the coverage term
    return p, big[p][1]
first = not os.path.exists(out)
with open(out, 'a') as o:
    if first: o.write('hap\tsample\tpop\tchr\tlength\told_tandex_n_mer\tnew_tandex_period\tnew_tandex_n_mer\tnew_tandex_identity\tself_n_mer\tself_windows\tconfident\told_agree\tnew_agree\n')
    done = set()
    if not first:
        for l in open(out):
            done.add(l.split('\t')[0])
    for d in dirs:
        if d in done: continue
        t0 = time.time()
        old = {}
        for r in csv.DictReader(open(f'{res}/{d}/units.tsv'), delimiter='\t'):
            old[(r['seq'], r['start'], r['end'])] = r
        tmp = tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False); n = 0; name = None; keep = False
        for l in gzip.open(f'{res}/{d}/arrays.fa.gz', 'rt'):
            if l[0] == '>':
                name = l[1:].strip(); f = name.split('|'); keep = (chrs is None or f[3] in chrs) and (f[0], f[1], f[2]) in old
                if keep: tmp.write(l); n += 1
            elif keep: tmp.write(l)
        tmp.close()
        if n == 0: os.unlink(tmp.name); continue
        p = subprocess.run([tandex, '-H', '-p', '50000', '-t', threads, tmp.name], capture_output=True, text=True)
        os.unlink(tmp.name)
        if p.returncode: print(d, 'TANDEX FAILED', p.stderr[-300:], flush=True); continue
        calls = collections.defaultdict(list)
        for l in p.stdout.splitlines():
            if l.startswith('#'): continue
            x = l.split('\t', 13); calls[x[0]].append((int(x[1]), int(x[2]), int(x[3]), float(x[6])))
        ok_old = ok_new = conf = 0
        for nm in calls.keys() | set():
            pass
        # every array, including ones with no call
        for l in gzip.open(f'{res}/{d}/arrays.fa.gz', 'rt'):
            if l[0] != '>': continue
            nm = l[1:].strip(); f = nm.split('|')
            if not ((chrs is None or f[3] in chrs) and (f[0], f[1], f[2]) in old): continue
            r = old[(f[0], f[1], f[2])]; L = int(f[2]) - int(f[1])
            pp, idn = unit(calls.get(nm, []), int(f[5]), L); t = round(pp / 171)
            na = 'yes' if str(t) == r['self_n_mer'] else 'no'
            o.write('\t'.join(map(str, [d, r['sample'], r['pop'], r['chr'], L, r['tandex_n_mer'], pp, t, '%.1f' % idn, r['self_n_mer'], r['self_windows'], r['confident'], r['agree'], na])) + '\n')
            if r['confident'] == 'yes':
                conf += 1; ok_old += r['agree'] == 'yes'; ok_new += na == 'yes'
        o.flush()
        print(f'{d}: {n} arrays; confident {conf}: agree before {ok_old}, now {ok_new} ({time.time() - t0:.0f} s)', flush=True)
print('written', out)
