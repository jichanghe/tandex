#!/usr/bin/env python3
"""Regression tests for tandex. usage: run_tests.py ./tandex"""
import gzip, os, random, subprocess, sys, tempfile

exe = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else './tandex'
rnd = random.Random(7)
tmp = tempfile.mkdtemp()
rs = lambda n: ''.join(rnd.choice('ACGT') for _ in range(n))

def mutate(s, d):
    out = []
    for c in s:
        if rnd.random() < d:
            r = rnd.random()
            if r < 0.8: out.append(rnd.choice([b for b in 'ACGT' if b != c]))
            elif r < 0.9: out.append(c + rnd.choice('ACGT'))
        else: out.append(c)
    return ''.join(out)

def write(name, records, gz=False):
    path = os.path.join(tmp, name)
    text = ''.join('>%s\n%s\n' % (n, '\n'.join(s[i:i + 70] for i in range(0, len(s), 70))) for n, s in records)
    (gzip.open(path, 'wt') if gz else open(path, 'w')).write(text)
    return path

def run(args, stdin=None):
    r = subprocess.run([exe] + args, input=stdin, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return [l.split('\t') for l in r.stdout.splitlines() if l and not l.startswith('#')]

def find(calls, seq, s, e, period, tol=0.05, cover=0.9):
    """a call on seq with the given period covering most of [s, e)"""
    for c in calls:
        cs, ce, p = int(c[1]), int(c[2]), int(c[3])
        if c[0] == seq and abs(p - period) <= tol * period and min(e, ce) - max(s, cs) >= cover * (e - s):
            return c
    return None

tests = []
def test(f): tests.append(f); return f

@test
def perfect_microsatellites():
    seq = rs(500) + 'CAG' * 30 + rs(500) + 'AC' * 25 + rs(500) + 'A' * 40 + rs(500)
    calls = run([write('a.fa', [('x', seq)])])
    c = find(calls, 'x', 500, 590, 3); assert c and c[10] in ('CAG', 'AGC', 'GCA') and float(c[6]) == 100.0, calls
    assert find(calls, 'x', 1090, 1140, 2) and find(calls, 'x', 1640, 1680, 1), calls

@test
def diverged_minisatellite():
    unit = rs(33); arr = mutate(unit * 20, 0.10)
    calls = run([write('b.fa', [('x', rs(800) + arr + rs(800))])])
    assert find(calls, 'x', 800, 800 + len(arr), 33, cover=0.8), calls

@test
def long_period_two_and_a_half_copies():
    unit = rs(3000); arr = (unit * 3)[:7500]
    calls = run([write('c.fa', [('x', rs(2000) + arr + rs(2000))])])
    c = find(calls, 'x', 2000, 9500, 3000); assert c and abs(float(c[4]) - 2.5) < 0.05, calls

@test
def higher_order_repeat_reports_both_levels():
    mono = rs(171); hor = ''.join(mutate(mono, 0.2) for _ in range(6))
    arr = ''.join(mutate(hor, 0.01) for _ in range(30))
    calls = run([write('d.fa', [('x', rs(3000) + arr + rs(3000))])])
    assert find(calls, 'x', 3000, 3000 + len(arr), len(hor), tol=0.02), calls
    assert find(calls, 'x', 3000, 3000 + len(arr), 171, tol=0.05, cover=0.8), calls

@test
def random_sequence_gives_almost_nothing():
    calls = run([write('e.fa', [('x', rs(2_000_000))])])
    assert len(calls) <= 3, len(calls)

@test
def lowercase_n_runs_crlf_and_empty_records():
    path = os.path.join(tmp, 'f.fa')
    open(path, 'w').write('>empty\n>tiny\nACGT\n>n desc\r\n' + 'N' * 50 + 'acgtt' * 12 + 'N' * 50 + '\r\n' + 'cag' * 20 + '\r\n')
    calls = run([path])
    assert find(calls, 'n', 50, 110, 5) and find(calls, 'n', 160, 220, 3), calls
    assert all(c[0] == 'n' for c in calls)

@test
def gzip_and_stdin_match_plain():
    seq = rs(1000) + 'GATA' * 40 + rs(1000)
    plain = run([write('g.fa', [('x', seq)])])
    assert plain == run([write('g.fa.gz', [('x', seq)], gz=True)])
    assert plain == run(['-'], stdin='>x\n' + seq + '\n')

@test
def bed_and_units_output():
    seq = rs(1000) + 'GATTACA' * 20 + rs(1000)
    u = os.path.join(tmp, 'u.bed')
    bed = run(['-f', 'bed', '-u', u, write('h.fa', [('x', seq)])])
    hit = [b for b in bed if int(b[6]) == 7]
    assert hit and hit[0][5] == '.' and len(hit[0]) == 9, bed
    s, e = int(hit[0][1]), int(hit[0][2])
    units = [l.split('\t') for l in open(u)]
    mine = [x for x in units if x[3] == 'x:%d-%d:7' % (s, e)]
    # units tile the call exactly; a chance match in the flank may add a short partial unit
    assert len(mine) in (20, 21) and int(mine[0][1]) == s and int(mine[-1][2]) == e, mine[:3]
    assert all(int(a[2]) == int(b[1]) for a, b in zip(mine, mine[1:]))
    assert sum(int(x[2]) - int(x[1]) == 7 for x in mine) >= 19

@test
def workers_match_single_process():
    parts = []
    for i in range(40):
        p = rnd.choice([2, 5, 17, 60, 400]); parts += [rs(120_000), mutate(rs(p) * rnd.randint(8, 40), 0.05)]
    seq = ''.join(parts) + rs(50_000)
    path = write('i.fa', [('x', seq)])
    a = {tuple(c[:4]) for c in run([path])}; b = {tuple(c[:4]) for c in run(['-t', '3', path])}
    assert len(a) >= 38 and len(a ^ b) <= max(1, len(a) // 50), (len(a), len(a ^ b))

@test
def evidence_and_score_filters():
    seq = rs(1000) + 'ACGTTGCATGCCAGT' * 2 + rs(1000)   # two copies of a 15-mer: score 60, evidence 30
    path = write('j.fa', [('x', seq)])
    assert find(run(['-e', '0', path]), 'x', 1000, 1030, 15), run(['-e', '0', path])
    assert find(run([path]), 'x', 1000, 1030, 15)
    assert not find(run(['-e', '45', path]), 'x', 1000, 1030, 15)
    assert not find(run(['-s', '80', path]), 'x', 1000, 1030, 15)

@test
def evalue_column_and_filter():
    seq = rs(200_000) + 'CAG' * 30 + rs(200_000) + 'ACGTTGCATGCCAGT' * 2 + rs(200_000)
    path = write('k.fa', [('x', seq)])
    calls = run([path])
    strong = find(calls, 'x', 200_000, 200_090, 3); weak = find(calls, 'x', 400_090, 400_120, 15)
    assert strong and weak and len(strong) == 14, calls
    assert float(strong[13]) < 1e-10 and float(weak[13]) > 1e-3, (strong[13], weak[13])
    kept = run(['-E', '1e-6', path])
    assert find(kept, 'x', 200_000, 200_090, 3) and not find(kept, 'x', 400_090, 400_120, 15)

@test
def pure_motif_next_to_related_repeat_is_kept():
    seq = rs(1000) + 'TTTC' * 14 + 'TC' * 110 + rs(1000)
    calls = run([write('l.fa', [('x', seq)])])
    assert find(calls, 'x', 1000, 1056, 4), calls
    assert find(calls, 'x', 1060, 1276, 2), calls

@test
def workers_with_small_sequence_after_large_array():
    # a long array (published to the shared board) followed by a short sequence handled in-process
    unit = rs(178); arr = mutate(unit * 30000, 0.03)
    path = write('m.fa', [('big', rs(100_000) + arr + rs(100_000)), ('small', rs(500) + 'CAG' * 30 + rs(500))])
    calls = run(['-t', '2', path])
    assert find(calls, 'big', 100_000, 100_000 + len(arr), 178, cover=0.95), [c[:4] for c in calls][:5]
    assert find(calls, 'small', 500, 590, 3), [c[:4] for c in calls if c[0] == 'small']

@test
def short_option_finds_small_microsatellites():
    seq = rs(5000) + 'CA' * 7 + rs(5000) + 'AAT' * 5 + rs(5000) + 'ACGTTGCATG' * 2 + rs(5000)
    path = write('n.fa', [('x', seq)])
    plain = run([path]); short = run(['-S', path])
    assert not find(plain, 'x', 5000, 5014, 2) and not find(plain, 'x', 10014, 10029, 3), plain
    assert find(short, 'x', 5000, 5014, 2) and find(short, 'x', 10014, 10029, 3), short
    # the lower bar applies to periods of 4 bp or less only
    assert not find(short, 'x', 15029, 15049, 10), short

@test
def hor_option_reports_the_full_unit():
    # a 12-monomer higher-order repeat; some copies lack two monomers, as in real centromeres
    rnd = random.Random(7)
    base = rs(171)
    def mutate(x, rate):
        return ''.join(rnd.choice('ACGT') if rnd.random() < rate else c for c in x)
    mono = [mutate(base, 0.25) for _ in range(12)]
    arr = []
    for i in range(60):
        m = mono if i % 7 else mono[:5] + mono[7:]
        arr.append(mutate(''.join(m), 0.01))
    seq = rs(3000) + ''.join(arr) + rs(3000)
    path = write('h.fa', [('x', seq)])
    end = 3000 + sum(len(a) for a in arr)
    def full_unit(rows):
        return [r for r in rows if r[0] == 'x' and abs(int(r[3]) - 2052) <= 20 and int(r[2]) - int(r[1]) > 0.9 * (end - 3000) and float(r[6]) >= 96]
    assert not full_unit(run(['-p', '5000', path])), 'default mode is not expected to follow the unit through the short copies'
    assert full_unit(run(['-p', '5000', '-H', path])), run(['-p', '5000', '-H', path])

failed = 0
for t in tests:
    try: t(); print('ok   ', t.__name__)
    except AssertionError as ex:
        failed += 1; print('FAIL ', t.__name__, str(ex)[:300])
print('%d/%d passed' % (len(tests) - failed, len(tests)))
sys.exit(1 if failed else 0)
