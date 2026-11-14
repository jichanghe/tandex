#!/usr/bin/env python3
"""Step 2: are the tandex calls that TRF (default settings) does not report real tandem repeats?
Independent checks that use no tandex code:
  1. unit test: edit-distance identity between consecutive repeat units (edlib),
     compared with 20 length- and period-matched control regions drawn from unannotated sequence;
  2. confirmation by other tools (TRF with lenient settings, ULTRA);
  3. overlap with the Broad tandem repeat catalog v1.0 (GRCh38)."""
import bisect, random, sys, edlib
rnd = random.Random(1)
seq = ''.join(l.strip() for l in open('chr21.fa') if l[0] != '>').upper()
N = len(seq)

def load(f, sc=None):
    out = []
    for l in open(f):
        if l[0] == '#' or l.startswith('SeqID'): continue
        x = l.rstrip('\n').split('\t')
        try: out.append((int(x[1]), int(x[2]), int(float(x[3])), float(x[sc]) if sc is not None else 0.0))
        except ValueError: pass
    return sorted(out)
def union(iv):
    o = []
    for s, e in sorted(iv):
        if o and s <= o[-1][1]: o[-1][1] = max(o[-1][1], e)
        else: o.append([s, e])
    return o
class Cov:
    def __init__(self, iv): self.u = union(iv); self.s = [x[0] for x in self.u]
    def bp(self, s, e):
        c = 0; i = max(0, bisect.bisect_right(self.s, s) - 1)
        while i < len(self.u) and self.u[i][0] < e:
            c += max(0, min(e, self.u[i][1]) - max(s, self.u[i][0])); i += 1
        return c
    def frac(self, s, e): return self.bp(s, e) / (e - s)

T = load('chr21.tandex.tsv', 5)
trf = Cov([(s, e) for s, e, p, _ in load('chr21.trf.bed', 4)])
trf30 = Cov([(s, e) for s, e, p, _ in load('chr21.trf30.bed', 4)])
ultra = Cov([(s, e) for s, e, p, _ in load('chr21.ultra.tsv')])
cat = Cov([(int(x[1]), int(x[2])) for x in (l.split('\t') for l in open('../truth/cat21.bed'))])
anycall = Cov([(s, e) for s, e, p, _ in T] + [tuple(x) for x in trf.u])  # controls avoid tandex and TRF-default calls

def shift_ident(s, e, p):
    """identity between consecutive repeat units, each aligned globally to the next one.
    (Units are separate segments, so unlike a plain shifted self-alignment an unrelated
    sequence gets no credit: random DNA scores about 0.45-0.5.)"""
    B = p if p >= 12 else p * -(-12 // p)          # compare blocks of >= 12 bp
    n = (e - s) // B
    pairs = max(1, min(n - 1, max(1, 1000 // B)))
    ed = tot = 0
    for j in range(pairs):
        a = seq[s + j * B:s + (j + 1) * B]; b = seq[s + (j + 1) * B:min(e, s + (j + 2) * B)]
        if len(b) < 6: break
        mode = 'NW' if len(b) == B else 'SHW'       # last unit may be partial
        ed += edlib.align(b, a, mode=mode, task='distance')['editDistance']; tot += len(b)
    return 1 - ed / tot if tot else None
def control(length, p):
    while True:
        s = rnd.randrange(5_000_000, N - length - 1)
        need = min(length, max(2 * p, 1000 + 2 * p))  # only the part the test actually reads
        if 'N' in seq[s:s + need] or anycall.bp(s - 50, s + need + 50): continue
        return shift_ident(s, s + length, p)
def passes(s, e, p, n=20):
    v = shift_ident(s, e, p)
    if v is None: return None
    return all(v > control(e - s, p) for _ in range(n))

only = [t for t in T if trf.frac(t[0], t[1]) < 0.5]
shared = [t for t in T if trf.frac(t[0], t[1]) >= 0.5]
print(f"tandex calls: {len(T)}; also in TRF default: {len(shared)}; tandex-only: {len(only)} ({sum(t[1]-t[0] for t in only)} bp)")

def report(name, calls):
    n = len(calls)
    if not n: return
    res = [passes(s, e, p) for s, e, p, sc in calls]
    tested = [r for r in res if r is not None]
    ok_shift = sum(tested)
    c30 = sum(trf30.frac(s, e) >= 0.5 for s, e, p, sc in calls)
    cu = sum(ultra.frac(s, e) >= 0.5 for s, e, p, sc in calls)
    cc = sum(cat.frac(s, e) >= 0.3 for s, e, p, sc in calls)
    anyc = sum((r is True) or trf30.frac(s, e) >= 0.5 or ultra.frac(s, e) >= 0.5 or cat.frac(s, e) >= 0.3
               for r, (s, e, p, sc) in zip(res, calls))
    nothing = [(c, r) for r, c in zip(res, calls) if not ((r is True) or trf30.frac(c[0], c[1]) >= 0.5 or ultra.frac(c[0], c[1]) >= 0.5 or cat.frac(c[0], c[1]) >= 0.3)]
    print(f"{name}\tn={n}\tbp={sum(c[1]-c[0] for c in calls)}\tunit_test={ok_shift}/{len(tested)} ({100*ok_shift/max(1,len(tested)):.1f}%)\t"
          f"TRF-lenient={100*c30/n:.1f}%\tULTRA={100*cu/n:.1f}%\tcatalog={100*cc/n:.1f}%\tany_support={100*anyc/n:.1f}%\tno_support={len(nothing)} ({sum(c[0][1]-c[0][0] for c in nothing)} bp)")
    return nothing

print("\nCalibration: tandex calls that TRF default also reports (random sample of 2000)")
report("shared", rnd.sample(shared, min(2000, len(shared))))
print("\ntandex-only calls by tandex score")
un = []
for lo, hi in [(40, 50), (50, 100), (100, 10**9)]:
    un += report(f"score {lo}-{hi if hi < 10**9 else 'max'}", [t for t in only if lo <= t[3] < hi]) or []
print("\ntandex-only calls by period")
for lo, hi in [(1, 6), (7, 100), (101, 2000), (2001, 10**9)]:
    report(f"period {lo}-{hi if hi < 10**9 else 'max'}", [t for t in only if lo <= t[2] <= hi])
print("\nLargest tandex-only calls with no independent support:")
for (s, e, p, sc), r in sorted(un, key=lambda x: -x[0][3])[:10]:
    print(f"  chr21:{s}-{e} period={p} score={sc:.0f} shift_identity={shift_ident(s,e,p)}")
