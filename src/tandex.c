/* tandex - fast de novo tandem repeat detector for DNA (FASTA in, TSV/BED out)
 *
 * Method
 *  1. Seeding: for four k-mer sizes (4, 6, 8, 11) every position votes for the distance d to the
 *     previous occurrence of its k-mer. Runs of votes for one distance (or a small window of
 *     distances) are chained into a candidate (anchor, period d).
 *  2. Refinement: a seed unit is cut from the cleanest part of the chain and both flanks are
 *     aligned against the endlessly repeated unit with an adaptive-band X-drop extension. The
 *     traceback votes a new consensus (substitutions, insertions, deletions) until it is stable.
 *  3. Selection: periods are reduced to their primitive divisor; competing periods over the same
 *     region are compared per base; a discarded call may survive as a trimmed high-purity core.
 *  4. Filtering: score, copy number and "evidence" = score minus what the first copy contributes
 *     (a unit always matches itself, so only further copies count as evidence of repetition).
 *
 * Build: gcc -O3 -o tandex tandex.c -lm -lz
 * Licence: MIT
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <zlib.h>
#include <time.h>

#define VERSION "0.8.0"

static int MA = 2, MI = 3, GP = 5, XD = 40, MINSCORE = 40, EVMIN = 30, BAND = 10, MAXP = 20000, KEEPALL = 0;
static int THREADS = 1, BEDOUT = 0;
static FILE *UNITS = NULL;
static double MINCOPIES = 1.8, MAXEVALUE = -1;
/* -S: periods of SHORTP or less need only SHORTSC score and SHORTEV evidence. Few periods share the
 * lower bar, so short microsatellites are found at a modest cost in chance calls. */
static int SHORT = 0;
#define SHORTP 4
#define SHORTSC 22
#define SHORTEV 14
#define MSC(p) ((SHORT && (p) <= SHORTP && MINSCORE > SHORTSC) ? SHORTSC : MINSCORE)
#define EVM(p) ((SHORT && (p) <= SHORTP && EVMIN > SHORTEV) ? SHORTEV : EVMIN)
/* per-sequence null model for the E-value: ln(chance calls per bp with evidence >= x) = EA + EB * x */
static double EA, EB, ENULL;
/* selection margins (tuned on simulated data and on the Broad TR catalog, human chr21) */
#define RED 0.97     /* a shorter primitive period must reach this share of the per-base score */
#define LONGER 1.05  /* a longer period must beat a shorter one by this factor */
#define COREF 1.2    /* a rescued core must be this much denser than the call that displaced it */
#define CHUNK ((int64_t)2000000)
#define NEG (INT_MIN / 4)
#define CAP ((int64_t)1 << 20)

static inline int64_t mod(int64_t a, int64_t p) { a %= p; return a < 0 ? a + p : a; }

/* ---------------- alignment against a periodic pattern ---------------- */
typedef struct { int (*v)[4]; int (*ins)[4]; int *del; } votes_t;
typedef struct { int64_t tlen, qlen, score, mat, mis, ins, del; } ext_t;

static int64_t ROWMAX = INT64_MAX;
/* extension may not run past these positions (set when a big array next door is already called) */
static int64_t LIMLO = 0, LIMHI = INT64_MAX;
/* per 32-kb region: windowed refinements inside a large call that only re-found a period already
 * called there. After REDMAX of them, further multiples of that call's period are skipped. */
static int *redun;
#define REDBIN (16 * 2048)
#define REDMAX 8
#define REDCOVER 100000
static uint8_t *TB; static size_t TBcap; static int64_t *LO; static size_t LOcap;
/* unit boundaries (per-copy coordinates), collected during traceback when BON is set */
static int BON = 0; static int64_t BPH, *BND; static size_t nB, mB;
static inline void bnd_push(int64_t x) { if (nB == mB) { mB = mB ? mB * 2 : 1024; BND = realloc(BND, mB * sizeof(int64_t)); } BND[nB++] = x; }

/* one DP segment; returns 1 if the extension terminated here */
static int dp_segment(const uint8_t *s, int64_t a, int dir, int64_t avail,
                      const uint8_t *P, int p, int64_t ph, votes_t *V, ext_t *o)
{
    int W = 2 * BAND + 1;
    int64_t R = avail < CAP ? avail : CAP;
    if ((size_t)(R + 1) * W > TBcap) { TBcap = (size_t)(R + 1) * W; TB = realloc(TB, TBcap); }
    if ((size_t)(R + 1) > LOcap) { LOcap = R + 1; LO = realloc(LO, LOcap * sizeof(int64_t)); }
    int bufA[64], bufB[64], *prev = bufA, *cur = bufB;
    LO[0] = 0;
    for (int c = 0; c < W; c++) { prev[c] = -GP * c; TB[c] = c ? 2 : 3; }
    int64_t bs = 0, bi = 0, bj = 0, pam = 0, i;
    int term = 0;
#define T(i) (dir > 0 ? s[a + (i)] : s[a - 1 - (i)])
#define PH(j) (dir > 0 ? mod(ph + (j), p) : mod(ph - 1 - (j), p))
    for (i = 1; i <= R; i++) {
        int64_t lo = pam + 1 - BAND, plo = LO[i - 1];
        if (lo < plo) lo = plo;
        LO[i] = lo;
        uint8_t tb = T(i - 1);
        int rmax = NEG; int64_t ram = lo;
        uint8_t *tbr = TB + (size_t)i * W;
        int64_t qph = lo >= 1 ? PH(lo - 1) : 0;
        for (int c = 0; c < W; c++) {
            int64_t j = lo + c, pc = j - 1 - plo;
            int best = NEG; uint8_t op = 0;
            if (j >= 1 && pc >= 0 && pc < W && prev[pc] > NEG) {
                best = prev[pc] + ((tb == P[qph] && tb < 4) ? MA : -MI);
            }
            pc++;
            if (pc >= 0 && pc < W && prev[pc] > NEG && prev[pc] - GP > best) { best = prev[pc] - GP; op = 1; }
            if (c > 0 && cur[c - 1] > NEG && cur[c - 1] - GP > best) { best = cur[c - 1] - GP; op = 2; }
            cur[c] = best; tbr[c] = op;
            if (best > rmax) { rmax = best; ram = j; }
            if (j >= 1) { qph += dir; if (qph >= p) qph = 0; else if (qph < 0) qph = p - 1; }
            else qph = PH(0);
        }
        if (rmax > bs) { bs = rmax; bi = i; bj = ram; }
        if (rmax < bs - XD) { term = 1; break; }
        pam = ram;
        int *t = prev; prev = cur; cur = t;
    }
    if (!term && R == avail) term = 1;
    /* traceback */
    int64_t ti = bi, j = bj; int lastop = -1;
    while (ti > 0 || j > 0) {
        uint8_t op = TB[(size_t)ti * W + (j - LO[ti])];
        if (op == 0) {
            uint8_t tb = T(ti - 1); int64_t q = PH(j - 1);
            if (tb == P[q]) o->mat++; else o->mis++;
            if (V && tb < 4) V->v[q][tb]++;
            if (BON && q == BPH) bnd_push(dir > 0 ? a + ti - 1 : a - ti);
            ti--; j--;
        } else if (op == 1) {
            uint8_t tb = T(ti - 1);
            o->ins++;
            if (V && tb < 4 && lastop != 1) V->ins[dir > 0 ? mod(ph + j, p) : mod(ph - j, p)][tb]++;
            ti--;
        } else if (op == 2) {
            o->del++;
            if (V) V->del[PH(j - 1)]++;
            if (BON && PH(j - 1) == BPH) bnd_push(dir > 0 ? a + ti : a - ti);
            j--;
        } else break;
        lastop = op;
    }
    o->tlen += bi; o->qlen += bj; o->score += bs;
    if (bi == 0) term = 1;
    return term;
}

/* Re-registration for long periods. Arrays of higher-order repeats hold copies that lack or repeat
 * whole sub-units, which shifts the phase by far more than the band. Before each segment the next
 * JWIN bases are compared with the pattern at the current phase; if another phase fits them clearly
 * better (at least 3 fewer mismatches, and at most 5), the alignment jumps to that phase at a cost of JPEN. */
#define JMIN 300
#define JWIN 64
#define JSEG 1024
#define JPEN 30
static int JUMPS = 0; /* set by -H */
static int win_mis(const uint8_t *s, int64_t x, const uint8_t *P, int p, int64_t q, int maxmis)
{
    int m = 0;
    q = mod(q, p);
    for (int i = 0; i < JWIN; i++) {
        if (s[x + i] != P[q] || s[x + i] > 3) { if (++m > maxmis) return m; }
        if (++q == p) q = 0;
    }
    return m;
}
static ext_t extend(const uint8_t *s, int64_t n, int64_t a, int dir, const uint8_t *P, int p, votes_t *V)
{
    ext_t o; memset(&o, 0, sizeof o);
    int64_t shift = 0; int jm = JUMPS && p >= JMIN;
    for (;;) {
        int64_t pos = dir > 0 ? a + o.tlen : a - o.tlen;
        int64_t avail = dir > 0 ? (n < LIMHI ? n : LIMHI) - pos : pos - LIMLO;
        int64_t ph = dir > 0 ? mod(o.qlen + shift, p) : mod(-o.qlen + shift, p);
        if (avail > ROWMAX - o.tlen) avail = ROWMAX - o.tlen;
        if (avail <= 0) break;
        if (jm) {
            if (o.tlen > 0 && avail >= JWIN) {
                int64_t x = dir > 0 ? pos : pos - JWIN, q0 = dir > 0 ? ph : ph - JWIN;
                int cur = JWIN;
                for (int dl = -3; dl <= 3; dl++) { int m = win_mis(s, x, P, p, q0 + dl, cur); if (m < cur) cur = m; }
                if (cur >= 3) { /* not a near-perfect fit here: is there a phase that fits clearly better? */
                    int bm = cur - 2 < 6 ? cur - 2 : 6; int64_t bq = -1;
                    for (int64_t q = 0; q < p; q++) { int m = win_mis(s, x, P, p, q, bm - 1); if (m < bm) { bm = m; bq = q; } }
                    if (bq >= 0) {
                        int64_t nph = dir > 0 ? bq : mod(bq + JWIN, p);
                        shift += nph - ph; ph = nph; o.score -= JPEN;
                    }
                }
            }
            int64_t before = o.tlen, seg = avail > JSEG ? JSEG : avail;
            int t = dp_segment(s, pos, dir, seg, P, p, ph, V, &o);
            if (o.tlen == before || (t && seg == avail)) break; /* no progress, or the end of what may be aligned */
            continue;
        }
        if (dp_segment(s, pos, dir, avail, P, p, ph, V, &o)) break;
    }
    return o;
}

/* ---------------- candidates / results ---------------- */
typedef struct { int64_t start, end, anchor; int d, hits, run; } cand_t;
typedef struct { int64_t start, end, score, mat, mis, indel, anchor; int p, aph, core, xd, d0, nested, imported, nrej, mrej, *rej; double copies, ident; uint8_t *cons; int dead; } res_t;

static cand_t *C; static size_t nC, mC;
static res_t *Rs; static size_t nR, mR;
#define DENSITY(x) ((double)(x)->score / (double)((x)->end - (x)->start))

typedef struct { int64_t first, last, bestpos; int cnt, run, bestrun, lastd, bestd; int64_t stamp; } chain_t;
/* ch[0]: exact-distance chains; ch[1], ch[2]: coarse (width w) chains, staggered by w/2 */
typedef struct { int k, dmin, dmax, minhits, w; int64_t *tab; chain_t *ch[3]; int *act[3]; int nact[3]; } tier_t;

static void emit(tier_t *t, chain_t *c)
{
    int d = c->bestd;
    if (c->cnt < t->minhits) return;
    if (c->last + t->k - c->first < d / 4) return; /* e.g. two copies of a dispersed repeat: support far shorter than the distance */
    if (nC == mC) { mC = mC ? mC * 2 : 1024; C = realloc(C, mC * sizeof(cand_t)); }
    cand_t *x = &C[nC++];
    x->start = c->first - d; if (x->start < 0) x->start = 0;
    x->end = c->last + t->k; x->anchor = c->bestpos - d; x->d = d; x->hits = c->cnt; x->run = c->bestrun;
}

/* hits are collected for k-mer starts in [cs, ce); earlier sequence only fills the k-mer table */
static void seed(tier_t *t, const uint8_t *s, int64_t n, int64_t base, int64_t stamp, int64_t cs, int64_t ce)
{
    int k = t->k; uint64_t mask = ((uint64_t)1 << (2 * k)) - 1, h = 0; int valid = 0;
    t->nact[0] = t->nact[1] = t->nact[2] = 0;
    int64_t e0 = cs - t->dmax - 64, e1 = ce + k - 1;
    if (e0 < 0) e0 = 0;
    if (e1 > n) e1 = n;
    for (int64_t e = e0; e < e1; e++) {
        if (s[e] > 3) { valid = 0; continue; }
        h = ((h << 2) | s[e]) & mask;
        if (++valid < k) continue;
        int64_t i = e - k + 1, g = base + i, prev = t->tab[h] - 1;
        t->tab[h] = g + 1;
        if (prev < base || i < cs) continue;
        int64_t d64 = g - prev;
        if (d64 > t->dmax || d64 < t->dmin) continue;
        int d = (int)d64;
        for (int f = 0; f < (t->w ? 3 : 1); f++) {
            int key = f == 0 ? d : f == 1 ? d / t->w : (d + t->w / 2) / t->w;
            chain_t *c = &t->ch[f][key];
            if (c->stamp != stamp) {
                c->stamp = stamp; t->act[f][t->nact[f]++] = key;
                c->first = i; c->cnt = 0; c->run = 0; c->bestrun = 0; c->last = -10;
            } else if (i - c->last > 2 * (int64_t)d + 30) {
                emit(t, c);
                c->first = i; c->cnt = 0; c->run = 0; c->bestrun = 0; c->last = -10;
            }
            c->cnt++;
            c->run = (i == c->last + 1 && d == c->lastd) ? c->run + 1 : 1;
            if (c->run > c->bestrun) { c->bestrun = c->run; c->bestpos = i - c->run + 1; c->bestd = d; }
            c->last = i; c->lastd = d;
        }
    }
    for (int f = 0; f < (t->w ? 3 : 1); f++)
        for (int x = 0; x < t->nact[f]; x++) emit(t, &t->ch[f][t->act[f][x]]);
}

/* grid index of results for overlap queries */
#define BIN 2048
typedef struct { int r, next; } node_t;
static int *head; static node_t *nodes; static size_t nN, mN;
static void grid_add(int r)
{
    for (int64_t b = Rs[r].start / BIN; b <= (Rs[r].end - 1) / BIN; b++) {
        if (nN == mN) { mN = mN ? mN * 2 : 4096; nodes = realloc(nodes, mN * sizeof(node_t)); }
        nodes[nN].r = r; nodes[nN].next = head[b]; head[b] = (int)nN++;
    }
}

static int near_multiple(int d, int q)
{
    int m = (int)((d + q / 2) / q);
    if (m < 1) return 0;
    double tol = 0.03 * d; if (tol < 1) tol = d < 20 ? 0 : 1; /* 4 and 5 are different motifs */
    return fabs((double)d - (double)m * q) <= tol ? m : 0;
}

/* failed refinements, remembered so that look-alike candidates are not retried forever */
typedef struct { int64_t a; int d, next; } fail_t;
static int *fhead; static fail_t *fails; static size_t nF, mF;
static void fail_add(int64_t a, int d)
{
    if (nF == mF) { mF = mF ? mF * 2 : 4096; fails = realloc(fails, mF * sizeof(fail_t)); }
    fails[nF].a = a; fails[nF].d = d; fails[nF].next = fhead[a / BIN]; fhead[a / BIN] = (int)nF++;
}

/* 0: refine, 1: skip, 2: refine only if a quick probe beats result *probe,
 * 3: refine only if a quick probe shows a real repeat at all (shorter period nested in *probe) */
static int covered(const cand_t *c, int64_t nb, int *probe)
{
    int64_t a = c->anchor < 0 ? 0 : c->anchor; int need = 0;
    for (int x = head[a / BIN]; x >= 0; x = nodes[x].next) {
        res_t *r = &Rs[nodes[x].r];
        if (r->start > a || r->end < a + c->d) continue;
        int m = near_multiple(c->d, r->p);
        if (m == 1 || near_multiple(c->d, r->d0) == 1) return 1;   /* same period, or a candidate like the one that gave r */
        for (int k = 0; k < r->nrej; k++) if (near_multiple(c->d, r->rej[k]) == 1) return 1; /* period already tried here and rejected */
        /* a multiple of r's period (a higher-order repeat?) with only short exact matches is not tried;
         * one with long exact matches goes on to the probe below, whatever r's identity */
        if (m > 1 && c->d > 100 && ((!JUMPS && r->ident >= 0.85) || c->run < 24)) return 1;
        /* longer period: must beat r in a quick probe. Similar or shorter period inside a clean array:
         * probe as well, so that off-period chains in a big satellite are not aligned end to end.
         * Inside a noisy call nothing is assumed and the candidate is refined. */
        if (c->d > r->p) { *probe = nodes[x].r; need = 2; }
        else if (r->ident >= 0.85) {
            if ((int64_t)c->d * 3 > (int64_t)r->p * 2) { *probe = nodes[x].r; need = 2; }
            else if (need == 0) { *probe = nodes[x].r; need = 3; }
        }
    }
    if (need) return need;
    int nf = 0;
    int64_t b0 = (a - 4 * (int64_t)c->d) / BIN, b1 = (a + 4 * (int64_t)c->d) / BIN;
    if (b0 < 0) b0 = 0;
    if (b1 >= nb) b1 = nb - 1;
    if (b1 - b0 > 64) { b0 = a / BIN - 32 < 0 ? 0 : a / BIN - 32; b1 = a / BIN + 32 >= nb ? nb - 1 : a / BIN + 32; }
    for (int64_t b = b0; b <= b1; b++)
        for (int x = fhead[b]; x >= 0; x = fails[x].next)
            if (llabs(fails[x].a - a) <= 4 * (int64_t)c->d && near_multiple(c->d, fails[x].d) == 1 && ++nf >= 2) return 1;
    return 0;
}

static int build_cons(const uint8_t *P, int p, votes_t *V, uint8_t *out)
{
    int n = 0;
    for (int ph = 0; ph < p; ph++) {
        int tot = 0, it = 0, bi = 0, bv = 0;
        for (int b = 0; b < 4; b++) { tot += V->v[ph][b]; it += V->ins[ph][b]; if (V->ins[ph][b] > V->ins[ph][bi]) bi = b; }
        if (it >= 2 && it * 2 > tot + V->del[ph]) out[n++] = bi;
        if (V->del[ph] > tot) continue;
        bv = P[ph] < 4 ? P[ph] : 0;
        for (int b = 0; b < 4; b++) if (V->v[ph][b] > V->v[ph][bv]) bv = b;
        out[n++] = bv;
    }
    return n;
}

/* cheap test for a longer-period candidate inside an existing call: align a few units against the
 * raw seed unit and require a clearly higher score density than the existing call has */
static double local_density(const res_t *q, const uint8_t *s, int64_t st, int64_t en);
static int probe_ok(const uint8_t *s, int64_t n, const cand_t *c, const res_t *r, int absolute)
{
    int p = c->d; int64_t a = c->anchor < 0 ? 0 : c->anchor;
    if (a + p > n) return 0;
    uint8_t *P = malloc(p);
    for (int i = 0; i < p; i++) P[i] = s[a + i] > 3 ? 0 : s[a + i];
    int xd0 = XD; XD += p / 8 < 400 ? p / 8 : 400;
    ROWMAX = 3 * (int64_t)p < 30000 ? 3 * (int64_t)p : (p > 30000 ? p : 30000);
    ext_t e = extend(s, n, a + p, 1, P, p, NULL);
    if (e.tlen < p) { ext_t l = extend(s, n, a, -1, P, p, NULL); e.tlen += l.tlen; e.score += l.score; }
    ROWMAX = INT64_MAX; XD = xd0; free(P);
    if (e.tlen < 0.9 * p) return 0;
    if (absolute) return (double)e.score / (double)e.tlen >= 0.2; /* about 65% identity to the raw seed unit */
    /* compare with how r's own motif does on the probed stretch, not with r's average */
    int64_t pe = a + p + e.tlen < n ? a + p + e.tlen : n;
    return (double)e.score / (double)e.tlen > 1.15 * local_density(r, s, a, pe);
}

static void refine(const uint8_t *s, int64_t n, const cand_t *c)
{
    int p = c->d, cap = 2 * p + 16;
    int64_t a = c->anchor; if (a < 0) a = 0; if (a + p > n) return;
    uint8_t *P = malloc(cap), *P2 = malloc(cap);
    memcpy(P, s + a, p);
    for (int i = 0; i < p; i++) if (P[i] > 3) P[i] = 0;
    int d0 = p; uint8_t *S = malloc(p); memcpy(S, P, p);
    ext_t L, Rr;
    votes_t V;
    int mi = MI, gp = GP, xd0 = XD, xd = XD + (c->d / 8 < 400 ? c->d / 8 : 400); /* longer units tolerate longer bad stretches */
    for (int iter = 0;; iter++) {
        if (iter == 0) { MI = mi > 1 ? mi - 1 : mi; GP = gp > 1 ? gp - 1 : gp; XD = xd; } /* lenient: seed unit is one noisy copy */
        else { MI = mi; GP = gp; XD = xd; }
        V.v = calloc(p, sizeof *V.v); V.ins = calloc(p, sizeof *V.ins); V.del = calloc(p, sizeof(int));
        if (iter == 0) { /* the seed unit matches itself: skip it, and bail out early on noise */
            Rr = extend(s, n, a + p, 1, P, p, &V); L = extend(s, n, a, -1, P, p, &V);
            if (Rr.qlen + L.qlen < (MINCOPIES - 1.1) * p || Rr.score + L.score < MSC(p) / 2) { free(V.v); free(V.ins); free(V.del); free(P); free(P2); free(S); MI = mi; GP = gp; XD = xd0; return; }
            for (int i = 0; i < p; i++) V.v[i][P[i]]++;
            Rr.tlen += p; Rr.qlen += p; Rr.mat += p; Rr.score += (int64_t)MA * p;
        } else { Rr = extend(s, n, a, 1, P, p, &V); L = extend(s, n, a, -1, P, p, &V); }
        int p2 = build_cons(P, p, &V, P2);
        free(V.v); free(V.ins); free(V.del);
        if (iter == 3 || p2 < 1 || p2 > cap || (iter > 0 && p2 == p && !memcmp(P, P2, p))) break;
        if (L.score + Rr.score < MSC(p) / 2) break;
        memcpy(P, P2, p2); p = p2;
    }
    MI = mi; GP = gp; XD = xd;
    /* primitive period */
    for (int q = 1; q <= p / 2; q++) {
        if (p % q) continue;
        int mism = 0;
        for (int i = 0; i < p; i++) mism += P[i] != P[(i + q) % p];
        if (mism > p / 3) continue;
        for (int i = 0; i < q; i++) {
            int cnt[4] = {0}, b = 0;
            for (int j = i; j < p; j += q) cnt[P[j]]++;
            for (int x = 1; x < 4; x++) if (cnt[x] > cnt[b]) b = x;
            P2[i] = b;
        }
        ext_t R2 = extend(s, n, a, 1, P2, q, NULL), L2 = extend(s, n, a, -1, P2, q, NULL);
        /* compare evidence (score beyond the first copy) per base */
        double sq = (double)(R2.score + L2.score - (int64_t)MA * q) / (double)(R2.tlen + L2.tlen + 1);
        double sp = (double)(Rr.score + L.score - (int64_t)MA * p) / (double)(Rr.tlen + L.tlen + 1);
        if (sq >= RED * sp && R2.score + L2.score >= MSC(q)) {
            memcpy(P, P2, q); p = q; Rr = R2; L = L2; break;
        }
    }
    int64_t score = L.score + Rr.score;
    double copies = (double)(L.qlen + Rr.qlen) / p;
    if (score >= MSC(p) && score - (int64_t)MA * p >= EVM(p) && copies >= MINCOPIES) {
        if (nR == mR) { mR = mR ? mR * 2 : 1024; Rs = realloc(Rs, mR * sizeof(res_t)); }
        res_t *r = &Rs[nR];
        r->start = a - L.tlen; r->end = a + Rr.tlen; r->score = score; r->p = p; r->copies = copies;
        r->mat = L.mat + Rr.mat; r->mis = L.mis + Rr.mis; r->indel = L.ins + L.del + Rr.ins + Rr.del;
        r->ident = (double)r->mat / (double)(r->mat + r->mis + r->indel);
        r->cons = malloc(p); r->dead = 0; r->anchor = a; r->aph = (int)(L.qlen % p); r->core = 0; r->xd = xd; r->d0 = c->d; r->nested = 0; r->imported = 0; r->nrej = r->mrej = 0; r->rej = NULL;
        int64_t ph = mod(-L.qlen, p);
        for (int i = 0; i < p; i++) r->cons[i] = P[(ph + i) % p];
        grid_add((int)nR++);
    }
    XD = xd0;
    /* The consensus can drift away from the seed unit when the extension runs into a neighbouring
     * repeat with a related motif (a pure (TTTC)12 next to a long (TC)n ends up as TC). If the
     * period changed, test the untouched seed unit with strict scoring and keep a clean hit. */
    if (p != d0 && d0 <= 2000) {
        int prim = 1;
        for (int q = 1; q <= d0 / 2 && prim; q++) {
            if (d0 % q) continue;
            int same = 1;
            for (int i = 0; i < d0 && same; i++) same = S[i] == S[(i + q) % d0];
            if (same) prim = 0;
        }
        if (prim) {
            int mi2 = MI, gp2 = GP;
            MI = 7; GP = 7; XD = 20;
            ext_t R3 = extend(s, n, a, 1, S, d0, NULL), L3 = extend(s, n, a, -1, S, d0, NULL);
            MI = mi2; GP = gp2; XD = xd0;
            int64_t mat = R3.mat + L3.mat, mis = R3.mis + L3.mis, ind = R3.ins + R3.del + L3.ins + L3.del;
            int64_t sc = (int64_t)MA * mat - (int64_t)MI * mis - (int64_t)GP * ind;
            double cp = (double)(R3.qlen + L3.qlen) / d0;
            if ((cp >= 3.0 || (cp >= MINCOPIES && mis + ind == 0)) && sc >= MSC(d0) && sc - (int64_t)MA * d0 >= EVM(d0) && mat * 10 >= (mat + mis + ind) * 9) {
                if (nR == mR) { mR = mR ? mR * 2 : 1024; Rs = realloc(Rs, mR * sizeof(res_t)); }
                res_t *r = &Rs[nR];
                r->start = a - L3.tlen; r->end = a + R3.tlen; r->score = sc; r->p = d0; r->copies = cp;
                r->mat = mat; r->mis = mis; r->indel = ind; r->ident = (double)mat / (double)(mat + mis + ind);
                r->cons = malloc(d0); r->dead = 0; r->anchor = a; r->aph = (int)(L3.qlen % d0); r->core = 1; r->xd = 20; r->d0 = c->d; r->nested = 0; r->imported = 0;
                r->nrej = r->mrej = 0; r->rej = NULL;
                int64_t ph = mod(-L3.qlen, d0);
                for (int i = 0; i < d0; i++) r->cons[i] = S[(ph + i) % d0];
                grid_add((int)nR++);
            }
        }
    }
    free(P); free(P2); free(S);
}

/* A call about to be discarded may still contain a short, much purer stretch of its own motif
 * (e.g. a perfect (TTTC)12 inside a long TC-rich region better described by another motif).
 * Re-extend from the anchor with strict scoring; on success the call is trimmed to that core. */
static int core_rescue(const uint8_t *s, int64_t n, res_t *r)
{
    int p = r->p, mi = MI, gp = GP, xd = XD;
    uint8_t *P = malloc(p);
    for (int i = 0; i < p; i++) P[i] = r->cons[(r->aph + i) % p];
    MI = 7; GP = 7; XD = 20;
    ext_t R = extend(s, n, r->anchor, 1, P, p, NULL), L = extend(s, n, r->anchor, -1, P, p, NULL);
    MI = mi; GP = gp; XD = xd;
    int64_t len = R.tlen + L.tlen, mat = R.mat + L.mat, mis = R.mis + L.mis, ind = R.ins + R.del + L.ins + L.del;
    int64_t sc = (int64_t)MA * mat - (int64_t)MI * mis - (int64_t)GP * ind;
    double copies = (double)(R.qlen + L.qlen) / p;
    if (len * 10 > (r->end - r->start) * 8 || sc < MSC(p) || sc - (int64_t)MA * p < EVM(p) || copies < MINCOPIES || (copies < 3.0 && mis + ind > 0)) { free(P); return 0; }
    r->start = r->anchor - L.tlen; r->end = r->anchor + R.tlen; r->score = sc; r->mat = mat; r->mis = mis; r->indel = ind;
    r->copies = copies; r->ident = (double)mat / (double)(mat + mis + ind); r->aph = (int)(L.qlen % p); r->core = 1;
    int64_t ph = mod(-L.qlen, p);
    for (int i = 0; i < p; i++) r->cons[i] = P[(ph + i) % p];
    free(P);
    return 1;
}

/* score per base that call q's motif achieves on [st, en): for short motifs measured directly
 * (best rotation, no gaps), otherwise q's overall density */
static double local_density(const res_t *q, const uint8_t *s, int64_t st, int64_t en)
{
    if (q->p > 64 || en - st < 8 || en - st > 2000) return DENSITY(q); /* longer stretches drift out of phase at indels */
    int64_t best = 0;
    for (int rot = 0; rot < q->p; rot++) {
        int64_t m = 0; int k = rot;
        for (int64_t i = st; i < en; i++) { m += s[i] == q->cons[k]; if (++k == q->p) k = 0; }
        if (m > best) best = m;
    }
    double id = (double)best / (double)(en - st);
    return MA * id - MI * (1.0 - id);
}

static int cmp_hits(const void *a, const void *b)
{
    const cand_t *x = a, *y = b;
    if (x->hits != y->hits) return y->hits - x->hits;
    return x->d - y->d;
}
static int cmp_score(const void *a, const void *b)
{
    const res_t *x = &Rs[*(const int *)a], *y = &Rs[*(const int *)b];
    if (x->score != y->score) return x->score > y->score ? -1 : 1;
    return x->p - y->p;
}
static int cmp_start(const void *a, const void *b)
{
    const res_t *x = &Rs[*(const int *)a], *y = &Rs[*(const int *)b];
    return x->start < y->start ? -1 : x->start > y->start;
}

/* A candidate anchored just outside a big, clean array with a compatible period would re-align
 * that whole array. Stop its extension a little inside the existing call instead. */
static void set_barriers(const cand_t *c, int64_t nb)
{
    int64_t a = c->anchor < 0 ? 0 : c->anchor, e = a + c->d;
    int64_t b0 = a / BIN - 8, b1 = e / BIN + 8;
    LIMLO = 0; LIMHI = INT64_MAX;
    if (b0 < 0) b0 = 0;
    if (b1 >= nb) b1 = nb - 1;
    for (int64_t b = b0; b <= b1; b++)
        for (int x = head[b]; x >= 0; x = nodes[x].next) {
            res_t *r = &Rs[nodes[x].r];
            if (r->end - r->start < 50000 || r->ident < 0.85) continue;
            if (r->start <= a && r->end >= e) continue;                 /* anchor inside r: handled by covered() */
            int64_t lo = c->d < r->p ? c->d : r->p, hi = c->d < r->p ? r->p : c->d;
            if (!(hi * 2 <= lo * 3 || near_multiple((int)hi, (int)lo) >= 1)) continue;
            /* a long candidate that is a multiple of r's period may be the higher-order repeat of r's array */
            if (JUMPS && c->d >= JMIN && c->d > r->p && c->run >= 24 && near_multiple(c->d, r->p) >= 2) continue;
            int64_t pad = 2 * hi;
            if (r->end <= e && r->end - pad > LIMLO) LIMLO = r->end - pad;
            if (r->start >= a && r->start + pad < LIMHI) LIMHI = r->start + pad;
        }
    if (LIMLO > a) LIMLO = a;
    if (LIMHI < e) LIMHI = e;
}

/* Shared board (workers only): large calls are published as soon as a worker finishes them, and
 * the other workers load them as read-only context. A worker whose chunk lies inside an array that
 * is already called then treats its own candidates as it would in a single process, instead of
 * aligning the whole array again. */
#define BOARDMAX 8192
#define BIGCALL 200000
typedef struct { int64_t start, end, score; int p, d0, owner; double ident; uint8_t cons[64]; } bent_t;
typedef struct { volatile int lock; volatile int n; bent_t e[BOARDMAX]; } board_t;
static board_t *BOARD; static int board_seen;
static void board_lock(void) { while (__sync_lock_test_and_set(&BOARD->lock, 1)) usleep(50); }
static void board_unlock(void) { __sync_lock_release(&BOARD->lock); }
static void board_publish(const res_t *r)
{
    if (!BOARD || r->end - r->start < BIGCALL) return;
    board_lock();
    if (BOARD->n < BOARDMAX) {
        bent_t *b = &BOARD->e[BOARD->n];
        b->start = r->start; b->end = r->end; b->score = r->score; b->p = r->p; b->d0 = r->d0; b->owner = (int)getpid(); b->ident = r->ident;
        if (r->p <= 64) memcpy(b->cons, r->cons, r->p);
        __sync_synchronize();
        BOARD->n++;
    }
    board_unlock();
}
static void board_import(void)
{
    if (!BOARD) return;
    int n = BOARD->n, me = (int)getpid();
    for (; board_seen < n; board_seen++) {
        bent_t *b = &BOARD->e[board_seen];
        if (b->owner == me) continue;
        if (nR == mR) { mR = mR ? mR * 2 : 1024; Rs = realloc(Rs, mR * sizeof(res_t)); }
        res_t *r = &Rs[nR];
        memset(r, 0, sizeof *r);
        r->start = b->start; r->end = b->end; r->score = b->score; r->p = b->p; r->d0 = b->d0; r->ident = b->ident;
        r->copies = (double)(b->end - b->start) / b->p; r->imported = 1;
        r->cons = malloc(b->p <= 64 ? b->p : 1);
        if (b->p <= 64) memcpy(r->cons, b->cons, b->p);
        grid_add((int)nR++);
    }
}

/* densest existing call that contains the candidate's anchor unit, or -1 */
static int best_cover(const cand_t *c)
{
    int64_t a = c->anchor < 0 ? 0 : c->anchor; int best = -1;
    for (int x = head[a / BIN]; x >= 0; x = nodes[x].next) {
        res_t *r = &Rs[nodes[x].r];
        if (r->start > a || r->end < a + c->d) continue;
        if (best < 0 || DENSITY(r) > DENSITY(&Rs[best])) best = nodes[x].r;
    }
    return best;
}

/* does an earlier call containing position a already have (about) period p? */
static int cover_has_period(int64_t a, int d, int p, int skip)
{
    for (int x = head[a / BIN]; x >= 0; x = nodes[x].next) {
        res_t *r = &Rs[nodes[x].r];
        if ((JUMPS ? nodes[x].r >= skip : nodes[x].r == skip) || r->start > a || r->end < a + d) continue; /* skip: results of this candidate itself */
        if (near_multiple(p, r->p) == 1) return 1;
    }
    return 0;
}

/* undo the most recent result (its grid nodes are the heads of their bins) */
static void drop_last(void)
{
    res_t *r = &Rs[nR - 1];
    for (int64_t b = r->start / BIN; b <= (r->end - 1) / BIN; b++) head[b] = nodes[head[b]].next;
    free(r->cons); free(r->rej); nR--;
}

static int cmp_anchor(const void *a, const void *b)
{
    const cand_t *x = a, *y = b;
    if (x->anchor != y->anchor) return x->anchor < y->anchor ? -1 : 1;
    if (x->d != y->d) return x->d - y->d;
    return y->hits - x->hits;
}

/* seeding + refinement for k-mer starts in [cs, ce); raw calls are appended to Rs */
static void find_raw(const uint8_t *s, int64_t n, tier_t *T, int nt, int64_t base, int64_t stamp, int64_t cs, int64_t ce)
{
    nC = nN = 0;
    int64_t nb = n / BIN + 2;
    head = malloc(nb * sizeof(int));
    fhead = malloc(nb * sizeof(int)); nF = 0;
    for (int64_t i = 0; i < nb; i++) head[i] = fhead[i] = -1;
    redun = calloc(n / REDBIN + 2, sizeof(int));
    for (int t = 0; t < nt; t++) seed(&T[t], s, n, base, stamp, cs, ce);
    if (nC) { /* exact and coarse chains often emit the same (anchor, period): keep the strongest */
        qsort(C, nC, sizeof(cand_t), cmp_anchor);
        size_t k = 0;
        for (size_t i = 0; i < nC; i++)
            if (k == 0 || C[i].anchor != C[k - 1].anchor || C[i].d != C[k - 1].d) C[k++] = C[i];
        nC = k;
        qsort(C, nC, sizeof(cand_t), cmp_hits);
    }
    board_seen = 0;
    for (size_t i = 0; i < nC; i++) {
        board_import();
        {
            int64_t ax = C[i].anchor < 0 ? 0 : C[i].anchor; int skip = 0;
            if (redun[ax / REDBIN] >= REDMAX)
                for (int x = head[ax / BIN]; x >= 0 && !skip; x = nodes[x].next) {
                    res_t *r = &Rs[nodes[x].r];
                    if (r->start <= ax && r->end >= ax + C[i].d && r->end - r->start >= REDCOVER && near_multiple(C[i].d, r->p) >= 1) skip = 1;
                }
            if (skip) continue;
        }
        int pr = -1, cv = covered(&C[i], nb, &pr);
        if (cv == 1 || (cv >= 2 && !probe_ok(s, n, &C[i], &Rs[pr], cv == 3))) continue;
        set_barriers(&C[i], nb);
        size_t r0 = nR;
        /* Inside a large existing call, first refine within a window. A candidate that fills the
         * window is one more description of the same array: it earns an end-to-end alignment only
         * if it is a clearly shorter period that fits about as well (or is a decent repeat in its
         * own right), or any period that fits clearly better. Otherwise it is dropped and its period is not tried there again. */
        int bc = best_cover(&C[i]);
        int64_t a0 = C[i].anchor < 0 ? 0 : C[i].anchor, W = 32 * (int64_t)C[i].d > 10000 ? 32 * (int64_t)C[i].d : 10000;
        if (bc >= 0 && Rs[bc].end - Rs[bc].start <= 4 * W) {
            /* the densest cover is a small local call, but the anchor may still lie inside a huge
             * array: then that array decides, or every such candidate is aligned end to end */
            int big = -1;
            for (int x = head[a0 / BIN]; x >= 0; x = nodes[x].next) {
                res_t *r = &Rs[nodes[x].r];
                if (r->start > a0 || r->end < a0 + C[i].d || r->end - r->start <= 4 * W) continue;
                if (big < 0 || r->end - r->start > Rs[big].end - Rs[big].start) big = nodes[x].r;
            }
            if (big >= 0) bc = big;
        }
        if (bc >= 0 && Rs[bc].end - Rs[bc].start > 4 * W) {
            int64_t lo = a0 - W > LIMLO ? a0 - W : LIMLO, hi = a0 + C[i].d + W < LIMHI ? a0 + C[i].d + W : LIMHI;
            int64_t klo = LIMLO, khi = LIMHI;
            LIMLO = lo; LIMHI = hi;
            refine(s, n, &C[i]);
            LIMLO = klo; LIMHI = khi;
            if (nR > r0) {
                res_t *r = &Rs[nR - 1], *q = &Rs[bc];
                int fills = (r->start <= lo + C[i].d && lo > klo) || (r->end >= hi - C[i].d && hi < khi);
                if (fills) {
                    int shorter = (int64_t)r->p * 3 <= (int64_t)q->p * 2;
                    /* a shorter period is another level of structure: keep it if it is a decent repeat in its own right */
                    int better = shorter ? (DENSITY(r) >= 0.95 * DENSITY(q) || DENSITY(r) >= 0.5) : DENSITY(r) > 1.15 * DENSITY(q);
                    if (cover_has_period(a0, C[i].d, r->p, JUMPS ? (int)r0 : (int)nR - 1)) better = 0; /* reduced to a period that is already called */
                    if (JUMPS) while (nR > r0) drop_last(); else drop_last();
                    if (better) refine(s, n, &C[i]);
                    else {
                        redun[a0 / REDBIN]++;
                        if (q->nrej == q->mrej) { q->mrej = q->mrej ? q->mrej * 2 : 16; q->rej = realloc(q->rej, q->mrej * sizeof(int)); }
                        q->rej[q->nrej++] = C[i].d; LIMLO = 0; LIMHI = INT64_MAX; continue;
                    }
                } else if (Rs[bc].end - Rs[bc].start >= REDCOVER && cover_has_period(a0, C[i].d, r->p, JUMPS ? (int)r0 : (int)nR - 1)) redun[a0 / REDBIN]++;
            }
        } else refine(s, n, &C[i]);
        LIMLO = 0; LIMHI = INT64_MAX;
        for (size_t k = r0; k < nR; k++) board_publish(&Rs[k]);
        if (nR == r0) fail_add(C[i].anchor < 0 ? 0 : C[i].anchor, C[i].d);
    }
    free(head); free(fhead); free(redun);
}

/* per-copy coordinates: redo the final alignment of a call and record where each unit starts */
static void write_units(const char *name, const uint8_t *s, int64_t n, const res_t *r)
{
    int p = r->p, mi = MI, gp = GP, xd = XD;
    uint8_t *P = malloc(p);
    for (int i = 0; i < p; i++) P[i] = r->cons[(r->aph + i) % p];
    if (r->core) { MI = 7; GP = 7; XD = 20; } else XD = r->xd;
    nB = 0; BON = 1; BPH = mod(-(int64_t)r->aph, p);
    extend(s, n, r->anchor, 1, P, p, NULL); extend(s, n, r->anchor, -1, P, p, NULL);
    BON = 0; MI = mi; GP = gp; XD = xd; free(P);
    bnd_push(r->start);
    for (size_t i = 1; i < nB; i++) { /* insertion sort is fine: two nearly sorted runs */
        int64_t v = BND[i]; size_t j = i;
        while (j > 0 && BND[j - 1] > v) { BND[j] = BND[j - 1]; j--; }
        BND[j] = v;
    }
    int64_t prev = -1; int idx = 0;
    for (size_t i = 0; i <= nB; i++) {
        int64_t b = i < nB ? BND[i] : r->end;
        if (b < r->start || b > r->end || b == prev) continue;
        if (prev >= 0) fprintf(UNITS, "%s\t%ld\t%ld\t%s:%ld-%ld:%d\t%d\n", name, (long)prev, (long)b, name, (long)r->start, (long)r->end, p, ++idx);
        prev = b;
    }
}

/* E-value calibration. tandex was run on i.i.d. random sequences of six base compositions (AT 50-75%,
 * 30 Mb each); the rate of calls with evidence >= x fell off log-linearly, with intercept and slope
 * depending on pm = sum of squared base frequencies (the chance that two random bases match):
 *   ln(rate per bp) = (2.96 - 21.7 pm) + (-1.419 + 3.938 pm) x
 * The E-value of a call is that rate times the sequence length: the number of calls with at least
 * this much evidence expected in a random sequence of the same length and composition. Real
 * sequence is not i.i.d., so treat it as a ranking aid, not an exact probability. */
static void set_null(const uint8_t *s, int64_t n)
{
    double cnt[4] = {0, 0, 0, 0}, tot = 0, pm = 0;
    for (int64_t i = 0; i < n; i++) if (s[i] < 4) cnt[s[i]]++;
    for (int b = 0; b < 4; b++) tot += cnt[b];
    if (tot < 1) tot = 1;
    for (int b = 0; b < 4; b++) pm += (cnt[b] / tot) * (cnt[b] / tot);
    if (pm < 0.25) pm = 0.25;
    if (pm > 0.33) pm = 0.33;
    EA = 2.96 - 21.7 * pm; EB = -1.419 + 3.938 * pm; ENULL = tot;
}
static double evalue(int64_t ev) { return ENULL * exp(EA + EB * (double)ev); }

/* redundancy removal over all raw calls of a sequence, then output */
static void finish(const char *name, const uint8_t *s, int64_t n)
{
    static const char *nt4 = "ACGT";
    int64_t nb = n / BIN + 2;
    head = malloc(nb * sizeof(int));
    for (int64_t i = 0; i < nb; i++) head[i] = -1;
    nN = 0;
    int *ord = malloc((nR + 1) * sizeof(int)), no = 0;
    for (size_t i = 0; i < nR; i++) ord[i] = (int)i;
    if (!KEEPALL) {
        if (nR) qsort(ord, nR, sizeof(int), cmp_score);
        for (size_t i = 0; i < nR; i++) {
            res_t *r = &Rs[ord[i]]; int drop = 0, qd = -1;
            int64_t len = r->end - r->start;
            for (int64_t b = r->start / BIN; b <= (r->end - 1) / BIN && !drop; b++)
                for (int x = head[b]; x >= 0; x = nodes[x].next) {
                    res_t *q = &Rs[nodes[x].r];
                    if (q->dead) continue;
                    int64_t ov = (r->end < q->end ? r->end : q->end) - (r->start > q->start ? r->start : q->start);
                    /* a longer period must have clearly more evidence than a shorter one over the same region
                     * (evidence, not raw score: a long consensus fits its own first copy for free);
                     * an intermediate level between two others needs a larger margin */
                    if (r->p * 3 <= q->p * 2 && ov * 10 >= (q->end - q->start) * 8 &&
                        (double)(q->score - (int64_t)MA * q->p) < (q->nested ? 1.5 : LONGER) * (double)(r->score - (int64_t)MA * r->p) * (len > q->end - q->start ? (double)(q->end - q->start) / (double)len : 1.0)) {
                        if (q->core) continue;
                        if (!(core_rescue(s, n, q) && DENSITY(q) > COREF * local_density(r, s, q->start, q->end))) q->dead = 1;
                        continue;
                    }
                    if (ov * 10 < len * 9) continue;
                    /* a clearly shorter period inside a longer-period call is another level of structure
                     * (motif inside a unit, monomer inside a higher-order repeat): keep it */
                    if (r->p * 3 <= q->p * 2) { r->nested = 1; continue; }
                    drop = 1; qd = nodes[x].r; break;
                }
            /* a core must be a different motif: never the displacing call's own period, and a multiple of
             * it only when nearly perfect (noisy multiples of a diverged repeat are not new motifs) */
            int mq = drop ? near_multiple(r->p, Rs[qd].p) : 0;
            if (drop && mq != 1 && (r->core || core_rescue(s, n, r)) && (mq == 0 || r->ident >= 0.97) && DENSITY(r) > COREF * local_density(&Rs[qd], s, r->start, r->end)) {
                drop = 0;
                int64_t l2 = r->end - r->start;
                for (int64_t b = r->start / BIN; b <= (r->end - 1) / BIN && !drop; b++)
                    for (int x = head[b]; x >= 0; x = nodes[x].next) {
                        res_t *q = &Rs[nodes[x].r];
                        if (q->dead) continue;
                        int64_t ov = (r->end < q->end ? r->end : q->end) - (r->start > q->start ? r->start : q->start);
                        int m2 = near_multiple(r->p, q->p);
                        if (ov * 10 >= l2 * 9 && (m2 == 1 || (m2 > 1 && r->ident < 0.97) || local_density(q, s, r->start, r->end) * COREF >= DENSITY(r))) { drop = 1; break; }
                    }
            }
            if (drop) r->dead = 1; else grid_add(ord[i]);
        }
    }
    for (size_t i = 0; i < nR; i++) if (!Rs[i].dead) ord[no++] = (int)i;
    if (no) qsort(ord, no, sizeof(int), cmp_start);
    for (int i = 0; i < no; i++) {
        res_t *r = &Rs[ord[i]];
        if (i > 0) { /* exact duplicates can arise from different chunks or anchors */
            res_t *q = &Rs[ord[i - 1]];
            if (q->start == r->start && q->end == r->end && q->p == r->p) continue;
        }
        long ev = (long)(r->score - (int64_t)MA * r->p);
        double E = evalue(ev);
        if (MAXEVALUE >= 0 && E > MAXEVALUE) continue;
        if (BEDOUT) {
            printf("%s\t%ld\t%ld\t", name, (long)r->start, (long)r->end);
            if (r->p <= 30) for (int j = 0; j < r->p; j++) putchar(nt4[r->cons[j]]); else printf("p%d", r->p);
            printf("\t%ld\t.\t%d\t%.1f\t%.1f\n", (long)(r->score > 1000 ? 1000 : r->score), r->p, r->copies, 100.0 * r->ident);
        } else {
            printf("%s\t%ld\t%ld\t%d\t%.1f\t%ld\t%.1f\t%ld\t%ld\t%ld\t", name, (long)r->start, (long)r->end, r->p,
                   r->copies, (long)r->score, 100.0 * r->ident, (long)r->mat, (long)r->mis, (long)r->indel);
            for (int j = 0; j < r->p; j++) putchar(nt4[r->cons[j]]);
            printf("\t%ld\t%s\t%.2g\n", ev, r->core ? "core" : "full", E);
        }
        if (UNITS) write_units(name, s, n, r);
    }
    for (size_t i = 0; i < nR; i++) { free(Rs[i].cons); free(Rs[i].rej); }
    free(ord); free(head);
}

/* TANDEX_PROFILE=1 in the environment prints where the wall time goes (stderr) */
static int PROF; static double PT[3];
static double wall(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }

static void process(const char *name, const uint8_t *s, int64_t n, tier_t *T, int nt, int64_t base, int64_t stamp)
{
    nR = 0;
    double t0 = wall(), t1;
    if (BOARD) { BOARD->n = 0; BOARD->lock = 0; } /* the board is per sequence */
    if (THREADS <= 1 || n < 2 * CHUNK) {
        find_raw(s, n, T, nt, base, stamp, 0, n);
    } else {
        /* worker processes share the sequence copy-on-write; each seeds one chunk but may extend a
         * repeat anywhere, so arrays that span chunks are found whole and merged in finish() */
        int nch = (int)((n + CHUNK - 1) / CHUNK), running = 0;

        FILE **tmp = malloc(nch * sizeof(FILE *));
        fflush(stdout); if (UNITS) fflush(UNITS);
        for (int c = 0; c < nch; c++) {
            tmp[c] = tmpfile();
            if (!tmp[c]) { perror("tmpfile"); exit(1); }
            if (running == THREADS) { wait(NULL); running--; }
            pid_t pid = fork();
            if (pid < 0) { perror("fork"); exit(1); }
            if (pid == 0) {
                int64_t cs = (int64_t)c * CHUNK, ce = cs + CHUNK < n ? cs + CHUNK : n;
                double w0 = wall();
                find_raw(s, n, T, nt, base, stamp, cs, ce);
                if (PROF && wall() - w0 >= 5) fprintf(stderr, "profile chunk\t%s\t%ld\t%.1f s\n", name, (long)cs, wall() - w0);
                for (size_t i = 0; i < nR; i++) if (!Rs[i].imported) { fwrite(&Rs[i], sizeof(res_t), 1, tmp[c]); fwrite(Rs[i].cons, 1, Rs[i].p, tmp[c]); }
                fflush(tmp[c]);
                _exit(0);
            }
            running++;
        }
        int st, bad = 0;
        while (wait(&st) > 0) if (!WIFEXITED(st) || WEXITSTATUS(st)) bad = 1;
        if (bad) { fprintf(stderr, "tandex: a worker failed on %s\n", name); exit(1); }
        for (int c = 0; c < nch; c++) {
            res_t r;
            rewind(tmp[c]);
            while (fread(&r, sizeof(res_t), 1, tmp[c]) == 1) {
                r.cons = malloc(r.p); r.rej = NULL; r.nrej = r.mrej = 0;
                if (fread(r.cons, 1, r.p, tmp[c]) != (size_t)r.p) { fprintf(stderr, "tandex: truncated worker output\n"); exit(1); }
                if (nR == mR) { mR = mR ? mR * 2 : 1024; Rs = realloc(Rs, mR * sizeof(res_t)); }
                Rs[nR++] = r;
            }
            fclose(tmp[c]);
        }
        free(tmp);
    }
    t1 = wall(); PT[0] += t1 - t0;
    set_null(s, n);
    finish(name, s, n);
    PT[1] += wall() - t1;
    if (PROF) fprintf(stderr, "profile seq\t%s\t%ld bp\tsearch %.1f s\tfinish %.1f s\n", name, (long)n, t1 - t0, wall() - t1);
}

static void usage(void)
{
    fprintf(stderr,
        "tandex %s - de novo tandem repeat detection\n"
        "Usage: tandex [options] <in.fa[.gz] | ->  > repeats.tsv\n"
        "  -p INT   max period [20000]\n"
        "  -s INT   min alignment score [40]\n"
        "  -e INT   min evidence = score - match*period, i.e. score beyond the first copy [30]\n"
        "  -E NUM   max E-value: chance calls with at least this evidence expected in a random\n"
        "           sequence of the same length and base composition [no limit]\n"
        "  -S       sensitive for short repeats: periods of 4 bp or less need only score 22 and\n"
        "           evidence 14 (finds microsatellites of about 11 bp and up; more chance calls)\n"
        "  -H       higher-order repeat mode: follow copies that lack or repeat whole sub-units, so that\n"
        "           the full repeat unit of satellite arrays is reported (about 2 times slower)\n"
        "  -c NUM   min copies [1.8]\n"
        "  -m/-x/-g INT  match/mismatch/gap [2/3/5]     -X INT  X-drop [40]     -b INT  band [10]\n"
        "  -t INT   worker processes for sequences longer than 4 Mb [1]\n"
        "  -f STR   output format: tsv or bed [tsv]\n"
        "  -u FILE  also write per-copy (unit) coordinates as BED\n"
        "  -a       keep all redundant calls\n"
        "  -v       version\n"
        "TSV: seq start end period copies score %%identity matches mismatches indels consensus evidence type evalue\n"
        "BED: seq start end motif|pPERIOD score(<=1000) . period copies %%identity\n", VERSION);
}

int main(int argc, char **argv)
{
    int o;
    const char *ufile = NULL;
    double tstart = wall(); PROF = getenv("TANDEX_PROFILE") != NULL;
    while ((o = getopt(argc, argv, "m:x:g:X:s:e:E:c:p:b:t:f:u:aSHvh")) != -1) {
        switch (o) {
        case 'm': MA = atoi(optarg); break;
        case 'x': MI = atoi(optarg); break;
        case 'g': GP = atoi(optarg); break;
        case 'X': XD = atoi(optarg); break;
        case 's': MINSCORE = atoi(optarg); break;
        case 'e': EVMIN = atoi(optarg); break;
        case 'E': MAXEVALUE = atof(optarg); break;
        case 'c': MINCOPIES = atof(optarg); break;
        case 'p': MAXP = atoi(optarg); break;
        case 'b': BAND = atoi(optarg); if (BAND > 31) BAND = 31; if (BAND < 2) BAND = 2; break;
        case 't': THREADS = atoi(optarg); if (THREADS < 1) THREADS = 1; break;
        case 'f': if (!strcmp(optarg, "bed")) BEDOUT = 1; else if (strcmp(optarg, "tsv")) { fprintf(stderr, "tandex: unknown format %s\n", optarg); return 1; } break;
        case 'u': ufile = optarg; break;
        case 'a': KEEPALL = 1; break;
        case 'S': SHORT = 1; break;
        case 'H': JUMPS = 1; break;
        case 'v': puts(VERSION); return 0;
        default: usage(); return o == 'h' ? 0 : 1;
        }
    }
    if (optind >= argc) { usage(); return 1; }
    if (MAXP < 1 || MA < 1 || MI < 1 || GP < 1) { fprintf(stderr, "tandex: scores and max period must be positive\n"); return 1; }
    gzFile f = strcmp(argv[optind], "-") ? gzopen(argv[optind], "r") : gzdopen(0, "r");
    if (!f) { perror(argv[optind]); return 1; }
    if (ufile && !(UNITS = fopen(ufile, "w"))) { perror(ufile); return 1; }
    if (THREADS > 1) {
        BOARD = mmap(NULL, sizeof(board_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
        if (BOARD == MAP_FAILED) BOARD = NULL; /* workers then simply do not share */
    }
    tier_t T[4] = {{4, 1, 10, 6, 0}, {6, 5, 64, 5, 0}, {8, 30, 600, 4, 8}, {11, 300, MAXP, 5, 32}};
    int nt = MAXP > 600 ? 4 : MAXP > 64 ? 3 : MAXP > 10 ? 2 : 1;
    for (int t = 0; t < nt; t++) {
        if (T[t].dmax > MAXP) T[t].dmax = MAXP;
        T[t].tab = calloc((size_t)1 << (2 * T[t].k), sizeof(int64_t));
        for (int k = 0; k < 3; k++) {
            T[t].ch[k] = calloc(T[t].dmax + 2, sizeof(chain_t));
            T[t].act[k] = malloc((T[t].dmax + 2) * sizeof(int));
        }
    }
    uint8_t code[256]; memset(code, 4, 256);
    code['A'] = code['a'] = 0; code['C'] = code['c'] = 1; code['G'] = code['g'] = 2; code['T'] = code['t'] = 3;
    code['U'] = code['u'] = 3;
    char name[1024] = ""; int nl = 0;
    uint8_t *s = NULL; size_t n = 0, m = 0; int64_t base = 0, stamp = 0;
    int have = 0, state = 0; /* state 0: line start, 1: header name, 2: rest of header, 3: sequence line */
    static char buf[1 << 16];
    if (!BEDOUT) puts("#seq\tstart\tend\tperiod\tcopies\tscore\tidentity\tmatches\tmismatches\tindels\tconsensus\tevidence\ttype\tevalue");
    for (;;) {
        int got = gzread(f, buf, sizeof buf);
        if (got < 0) { fprintf(stderr, "tandex: read error\n"); return 1; }
        if (got == 0) break;
        if (n + got > m) { m = (n + got) * 2 + 1024; s = realloc(s, m); }
        for (int i = 0; i < got; i++) {
            char ch = buf[i];
            if (ch == '\n') { state = 0; continue; }
            if (state == 0) {
                if (ch == '>') {
                    if (have && n) { process(name, s, n, T, nt, base, ++stamp); base += n + 1; }
                    n = 0; nl = 0; name[0] = 0; have = 1; state = 1; continue;
                }
                state = 3;
            }
            if (state == 1) {
                if (ch == ' ' || ch == '\t' || ch == '\r') state = 2;
                else if (nl < 1023) { name[nl++] = ch; name[nl] = 0; }
            } else if (state == 3) {
                if (!have) { have = 1; strcpy(name, "seq"); }
                if (ch > ' ') s[n++] = code[(uint8_t)ch];
            }
        }
    }
    if (have && n) process(name, s, n, T, nt, base, ++stamp);
    gzclose(f);
    if (UNITS) fclose(UNITS);
    if (PROF) fprintf(stderr, "profile total\t%.1f s\tsearch %.1f s\tfinish %.1f s\tread and other %.1f s\n", wall() - tstart, PT[0], PT[1], wall() - tstart - PT[0] - PT[1]);
    return 0;
}
