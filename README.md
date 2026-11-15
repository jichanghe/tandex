# tandex

Fast de novo tandem repeat detection in DNA sequences. FASTA (plain or gzip) in, TSV or BED out.
One C file, no dependencies beyond zlib.

Repository: https://github.com/jichanghe/tandex

[![tests](https://github.com/jichanghe/tandex/actions/workflows/test.yml/badge.svg)](https://github.com/jichanghe/tandex/actions/workflows/test.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24256789.svg)](https://doi.org/10.5281/zenodo.24256789)

```
git clone https://github.com/jichanghe/tandex && cd tandex
make            # or: gcc -O3 -o tandex src/tandex.c -lm -lz
make test       # 14 regression tests
./tandex genome.fa.gz > repeats.tsv
./tandex -t 8 -f bed -u units.bed genome.fa > repeats.bed
```

## Options

| Option | Meaning | Default |
|---|---|---|
| `-p INT` | maximum period (use 50000 to capture human rDNA arrays) | 20000 |
| `-s INT` | minimum alignment score | 40 |
| `-e INT` | minimum evidence (score beyond the first copy, see below) | 30 |
| `-E NUM` | maximum E-value (see Method, step 5) | no limit |
| `-S` | sensitive for short repeats: periods of 4 bp or less need only score 22 and evidence 14 | off |
| `-c NUM` | minimum copy number | 1.8 |
| `-m/-x/-g INT` | match / mismatch / gap scores | 2 / 3 / 5 |
| `-X INT`, `-b INT` | X-drop (grows with period), band half-width | 40, 10 |
| `-t INT` | worker processes (used for sequences longer than 4 Mb) | 1 |
| `-f tsv\|bed` | output format | tsv |
| `-u FILE` | also write the coordinates of every repeat copy (BED) | off |
| `-a` | keep all redundant calls | off |

**TSV columns** (0-based, half-open): `seq start end period copies score %identity matches
mismatches indels consensus evidence type evalue`. `type` is `full`, or `core` for a call trimmed
to its high-purity part (see Method, step 3).

**BED columns**: `seq start end name score . period copies %identity`; `name` is the motif when
it is 30 bp or shorter, otherwise `p<period>`; score is capped at 1000.

**Units file** (`-u`): `seq unit_start unit_end call_id copy_number`, where `call_id` is
`seq:start-end:period`. Units tile each call exactly.

## Method

1. **Seeding.** For four k-mer sizes (4, 6, 8, 11) each position votes for the distance to the
   previous occurrence of its k-mer. Runs of votes for one distance, or for a small window of
   distances so that indels do not split the signal, become candidates. Linear time, about 1-3%
   of the run time.
2. **Refinement.** A seed unit is cut from the cleanest stretch of the candidate. Both flanks are
   aligned to the endlessly repeated unit with an adaptive-band X-drop extension. The traceback
   votes a new consensus (substitutions, insertions, deletions) and the alignment is repeated
   until the consensus is stable. Cost depends on repeat content, not on the period limit.
3. **Selection.** Periods are reduced to their primitive divisor. When several periods explain
   the same region, they are compared by evidence per base, and a longer period must be clearly
   better. A clearly shorter period inside a longer one is kept as a second level of structure
   (motif inside a unit, monomer inside a higher-order repeat). A call that loses this
   comparison is re-extended with strict scoring; if a much purer stretch remains, it is
   reported as a `core` (for example a perfect (TTTC)12 inside a long TC-rich region).
   To stay fast on large satellite arrays, a candidate inside or next to an array that is
   already called cleanly is first probed over a few units, and its extension stops at the
   existing call instead of re-aligning it. Inside a large array, an alternative period is first
   aligned in a window; it earns an end-to-end alignment only if it is a clearly shorter period
   that fits about as well, or any period that fits clearly better.
   When the consensus drifts away from the seed unit because the extension ran into a
   neighbouring repeat with a related motif, the untouched seed unit is tested with strict
   scoring and a clean hit is kept as a `core` too.
4. **Filtering by evidence.** A repeat unit always matches itself, so the first copy says
   nothing about repetition. `evidence = score - match x period` counts only what further copies
   add. Requiring evidence of 30 removes almost all chance calls (see below) without losing
   sensitivity on the simulation.

5. **E-value.** Each call gets the number of calls with at least its evidence that would be
   expected in a random sequence of the same length and base composition. The model was fitted
   to tandex runs on random sequences of six compositions (AT 50-75%) and checked on fresh ones:
   at E-value cut-offs of 10, 30, 100 and 300 the observed counts were 6-11, 22-41, 78-116 and
   227-310. Real sequence is not random, so use it to rank calls, not as an exact probability.

`-t N` forks worker processes that share the sequence. Each seeds one 2 Mb chunk but may extend
a repeat anywhere, so arrays that span chunks are found whole and merged afterwards. Workers
publish calls longer than 200 kb to a shared board and load each other's as read-only context. Results
differ from a single-process run by about one call in 10,000.

## Benchmarks

Machine: 2 cores, 7 GB RAM. All timings single-process unless marked. Results are for version
0.6.0 unless a table says otherwise; 0.6.0 differs from 0.5.0 only when `-t` is used, and 0.7.0
differs from 0.6.0 only when `-S` is used (default output on chr21 is identical). Tool versions: TRF 4.10,
ULTRA 1.2.2, tantan 49, longdust (git, Nov 2026), SRF (git), mreps 2.6, pytrf 1.5.0.
Scripts are in `bench/`.

### 1. Simulation with known truth (held-out set)

15.2 Mb random sequence, 2,000 implanted repeats: periods 1 bp to 6 kb, 2.5 to 50 copies, 0 to
25% divergence. "Found" = at least half the locus covered; "period right" = a covering call has
the true period within 5%.

| Tool (settings) | Time | Peak RAM | Found | Found, period right | Truth bases covered | Calls per found locus |
|---|---|---|---|---|---|---|
| **tandex** | **7.1 s** | 59 MB | **89.5%** | **88.7%** | **99.7%** | **1.08** |
| TRF (2 7 7 80 10 50 2000) | 30.8 s | 102 MB | 58.2% | 57.4% | 24.1% | 1.34 |
| TRF (2 5 7 80 10 50 2000) | 33.9 s | 110 MB | 63.6% | 62.6% | 25.5% | 1.34 |
| ULTRA (-p 500, 2 threads) | 934 s | 4260 MB | 63.9% | 55.9% | 11.0% | 1.51 |
| tantan (-w 2000) | 238 s | 203 MB | 53.3% | 42.5% | 21.6% | 2.24 |
| tantan (default) | 8.8 s | 35 MB | 39.4% | 30.1% | 2.9% | 2.65 |
| mreps (-res 3) | 558 s | 406 MB | 50.4% | 41.4% | 19.4% | 5.18 |
| longdust | 22.4 s | 32 MB | 38.5% | reports no period | 10.4% | 1.13 |
| pytrf (ATR + GTR, motif up to 500) | 52 s | 65 MB | 29.8% | 29.2% | 3.2% | 2.48 |

Restricted to loci each tool can represent (found, period right): periods up to 2,000: tandex
87.4%, TRF 64.0% / 69.9%, ULTRA 62.4%, tantan -w 2000 47.5%. Periods up to 500: tandex 86.0%,
TRF 64.2% / 70.4%, ULTRA 70.4%, tantan -w 2000 45.5%.

tandex by divergence: 100 / 100 / 99 / 92 / 86 / 80 / 62% at 0 / 2 / 5 / 10 / 15 / 20 / 25%.

### 2. Chance calls on 20 Mb of random DNA

| tandex | TRF 2/7/7 | TRF 2/5/7 | tantan | ULTRA (-p 500) |
|---|---|---|---|---|
| **3** | 16 | 54 | 788 | 38 in the first 3 Mb |

Before the evidence filter tandex made 121.

### 3. Human chr21 (GRCh38) against the Broad tandem repeat catalog v1.0

66,298 catalog loci on chr21. "Exact" = a covering call has the catalog's motif length.
"Compatible" = exact, or one motif is the other repeated with at least 80% identity (for example
CCTTCCTTCCTTCGTT against CCTT).

| Tool | Time | Calls | 25-49 bp loci (n=3,223): found / exact / compatible | 50-99 bp (n=539) | 100 bp and longer (n=182) | 9-24 bp found |
|---|---|---|---|---|---|---|
| tandex | 36 s | 21,650 | 99.3 / 93.9 / 96.5% | 100 / 85.9 / 91.8% | 100 / 84.6 / 89.0% | 36.1% |
| TRF default | 90 s | 14,860 | 99.6 / **94.3** / **97.7**% | 99.8 / **89.4** / **93.7**% | 100 / 86.8 / 91.8% | 20.3% |
| TRF lenient (2 5 7, score 30) | 152 s | 64,694 | 99.9 / 93.6 / 97.7% | 99.8 / 87.2 / 92.4% | 100 / 86.3 / 90.7% | 45.2% |
| tantan | **25.6 s** | 45,343 | 98.4 / 85.9 / 93.9% | 99.1 / 85.7 / 93.7% | 99.5 / **89.6** / **92.9**% | 59.6% |
| ULTRA -p 100 | 688 s | 42,753 | 99.9 / 75.6 / 92.3% | 99.8 / 77.0 / 86.3% | 99.5 / 86.8 / 89.0% | **59.8%** |
| longdust | 11.5 s | 17,999 | 98.0% found | 95.4% found | 96.7% found | 40.7% |

Every tool finds nearly all catalog loci of 25 bp or more. On motif agreement **TRF is still
ahead of tandex**, now by 0.4 to 3.5 points (exact); it was 2.9 to 11.1 points in version 0.2 and
1.9 to 7.1 in version 0.4. The gain comes from keeping clean seed-unit repeats (`core` calls,
about 2,000 of the 21,650). There is a trade-off: allowing more cores raises agreement with the
catalog further (to 94.2 / 87.4 / 88.5%, ahead of TRF on the longest loci) but adds redundant
calls on the simulation, so the stricter setting was kept. The remaining disagreements are mostly
two-copy 25-40 bp motifs inside AT-rich regions, which TRF and tantan cover by reporting several
overlapping periods per locus. The catalog was built partly from TRF-style annotations.

### 4. Are tandex's extra calls on chr21 real?

9,853 tandex calls (1.04 Mb) are not in TRF's default output. Each was tested without tandex
code (`bench/validate.py`): consecutive repeat units must be more similar than in 20 matched
control regions (chance pass rate about 5%); "any support" also counts coverage by TRF lenient,
by ULTRA, or by the catalog.

| tandex-only calls | n | Pass unit test | Any support | No support |
|---|---|---|---|---|
| calibration: calls TRF also makes (sample) | 2,000 | 98.5% | 100% | 0 |
| score 40-49 | 5,043 | 92.1% | 99.4% | 30 |
| score 50-99 | 3,058 | 94.7% | 99.3% | 21 |
| score 100 and up | 1,752 | 99.4% | 99.7% | 5 |

56 calls (5.9 kb) have no support, down from 245 in version 0.2. The cost of the evidence
filter: tandex now misses 4.9% of TRF's records (3.0% before); all 735 have a TRF score
below 100. Use `-e 0` to get them back.

### 5. Arabidopsis Col-CEN v1.2: centromere satellite

Truth: a tool-independent mask of the 178-bp satellite from direct autocorrelation (1-kb windows
where at least 60% of bases equal the base 176-180 bp downstream; background 29%). 10.9 Mb.

Whole genome (132 Mb):

| Tool | Time | Mask covered | With period 176-180 | Calls in mask |
|---|---|---|---|---|
| tandex | **123 s** (69 s with `-t 2`) | 99.9% | 99.5% | 1,498 |
| tantan (-w 200) | 146 s | 99.7% | 99.6% | 1,499 |
| TRF default | not finished after 55 min | - | - | - |

Version 0.3 needed 288 s here; 0.4 no longer re-aligns a large array once per stray candidate.

Chr1 centromere (3 Mb extract, 2.1 Mb of satellite):

| Tool | Time | Peak RAM | Mask covered with period n x 178 | Calls | Longest call |
|---|---|---|---|---|---|
| tandex | 4.7 s | 73 MB | 100% | 34 | 1.55 Mb |
| tantan (-w 200) | **4.0 s** | 41 MB | 99.9% | 273 | 546 kb |
| ULTRA (-p 200) | 102 s | 1690 MB | 100% | 51 | 176 kb |
| SRF + minimap2 | 38 s | - | 99.5% | 5 satellite units (178, 356, 534 bp) | - |
| longdust | 5.2 s | 13 MB | 100% covered, no period | 3 | 2.14 Mb |
| TRF default | not finished after 560 s | 1239 MB | - | - | - |

All finishing tools recover the satellite. tandex and tantan are equally fast on this extract;
tandex gives the least fragmented annotation.

### 6. T2T-CHM13 v2.0 against the consortium's satellite annotation

Truth: the censat v2.1 track (alpha-satellite HOR arrays, divergent HORs, monomeric alpha,
HSat1A/1B/2/3, beta and gamma satellite, other centromeric satellites, rDNA). The tables give the
share of annotated satellite bases covered by each tool's calls.

Whole genome (24 chromosomes, 3.1 Gb, 240.5 Mb of annotated satellite); tandex with
`-p 50000 -t 2`, tantan with `-f 4 -w 200`:

| Class | Annotated Mb | tandex | tantan |
|---|---|---|---|
| alpha HOR arrays | 70.33 | 100.0% | 8.7% |
| divergent HORs | 1.85 | 93.7% | 25.8% |
| monomeric alpha | 13.49 | 85.7% | 45.8% |
| HSat1A / HSat1B | 13.39 / 15.32 | 100.0% / 100.0% | 67.7% / 11.0% |
| HSat2 | 28.71 | 100.0% | 31.7% |
| HSat3 | 69.33 | 100.0% | 39.2% |
| beta satellite | 8.61 | 81.8% | 20.6% |
| gamma satellite | 0.65 | 95.7% | 6.1% |
| other centromeric satellite | 8.91 | 74.6% | 20.4% |
| rDNA | 9.93 | 98.7% | 16.6% |
| **all classes** | 240.52 | **97.5%** | 27.1% |
| alpha satellite with period n x 171 | 85.67 | 97.5% | 14.8% |

tandex makes 1.61 million calls and puts 381 Mb of the genome in tandem repeats. tantan's
figures are for its repeat-unit output, which is what reports a period; its masking mode covers
more (68.5% on chr21 against 32.8%) but gives no period.

**The uncovered part of the weaker classes is mostly not satellite.** On five chromosomes the
500-bp windows of "monomeric alpha" annotation that tandex leaves uncovered were tested for
171-bp periodicity directly (share of bases equal to the base one monomer downstream): 0-1% of
them pass, against 83-98% of covered windows. The same holds for beta satellite on chr21 and
chr22 (0% against 86-92%). The gaps are 4-8 kb blocks, the size of transposon insertions, that
the annotation includes in the satellite interval. The low "other centromeric satellite" figure
comes mainly from chrY, where the annotation covers composite TSPY/RBMY gene arrays. No change
was made to tandex for these classes.

**Speed.** The whole genome takes 47 minutes on 2 cores (`-p 50000 -t 2`, clean run, peak memory
337 MB); the slowest chromosome is chr1 at 5.0 minutes and the median is 2.0 minutes. Version 0.4
needed 82 minutes, with 14 and 20 minutes for chr1 and chr16, because every alternative period was
aligned across their 13 Mb HSat2 arrays, once per worker. Version 0.5 aligns alternative periods
in a window first, and version 0.6 lets workers publish large calls to a shared board so that
other workers do not align the same array again (chr16: 20 min, then 4.0 min, now 2.2 min).
Coverage is identical in all three versions.

chr21 (45 Mb, 7.79 Mb of annotated satellite):

| Class | Mb | tandex | tandex `-p 50000` | TRF | ULTRA -p 200 | tantan -w 200 |
|---|---|---|---|---|---|---|
| alpha HOR arrays | 0.43 | 100% | 100% | 100% | 100% | 14.5% |
| monomeric alpha | 0.65 | 88.4% | 88.4% | 89.4% | 86.4% | 66.8% |
| HSat1A | 0.91 | 100% | 100% | 100% | 100% | 61.1% |
| HSat1B | 0.07 | 100% | 100% | 51.9% | 60.1% | 18.4% |
| HSat3 | 1.60 | 99.8% | 99.8% | 98.7% | 97.6% | 42.6% |
| beta satellite | 1.03 | 93.3% | 93.3% | 81.3% | 80.9% | 28.5% |
| gamma satellite | 0.03 | 94.8% | 94.8% | 75.6% | 59.5% | 16.6% |
| other centromeric satellite | 0.55 | 90.7% | 90.7% | 61.0% | 55.5% | 19.7% |
| rDNA (45 kb unit) | 2.50 | 35.4% | **98.0%** | 24.0% | 32.3% | 15.8% |
| **all classes** | 7.79 | 76.7% | **96.8%** | 68.6% | 70.4% | 32.8% |
| alpha satellite with period n x 171 | 1.09 | 92.8% | 92.8% | 93.5% | 91.5% | 45.5% |
| **time** | | 94 s | 100 s (80 s in v0.5) | 904 s | 1,976 s | 52 s |

Other chromosomes (all satellite classes):

| Chromosome | Annotated satellite | tandex | tantan -w 200 | tandex time | tantan time |
|---|---|---|---|---|---|
| chr22 | 12.6 Mb | 92.7% (97.5% with `-p 50000`) | 27.1% | 138 s | 60 s |
| chr9 (28 Mb of HSat3) | 31.6 Mb | 99.5% | 41.6% | 280 s | 176 s |
| chrX | 3.9 Mb | 97.5% | not run | 215 s | - |

(Times in these two tables are single-process runs.)

Reading: tandex covers the most satellite sequence in every class but monomeric alpha, where
TRF is 1 point ahead, and it is 10 to 20 times faster than TRF and ULTRA. The rDNA arrays have a
45 kb unit, above the default period limit; `-p 50000` finds them (44.8 kb units at 99.7-100%
identity) at almost no extra cost. tantan is about twice as fast as tandex here but its window
of 200 bp cannot see higher-order repeat units of 1-3 kb, and larger windows are much slower.
TRF finished chr21 in 15 minutes but did not finish chr22 within 45 minutes. (The chr21 and
chr22 columns were measured with version 0.4; version 0.5 gives the same coverage.)

### 7. Short microsatellites: the `-S` option (version 0.7.0)

By default tandex needs a score of 40, so it skips repeats shorter than about 20 bp. `-S` lowers
the bar for periods of 4 bp or less only. Catalog loci of 9 to 24 bp on GRCh38 chr21, against
chance calls on 20 Mb of random sequence:

| Tool | 9-24 bp loci found | Chance calls, uniform random | Chance calls, 70% AT random |
|---|---|---|---|
| tandex default | 36.1% | 3 | 126 |
| **tandex -S** | **59.3%** | 291 | 2,385 |
| tantan | 59.6% | 788 | 4,663 |
| ULTRA -p 100 | 59.8% | about 250 (38 in the first 3 Mb, at -p 500) | not run |

With `-S` tandex finds as many short loci as tantan and ULTRA, with fewer chance calls than
tantan and about as many as ULTRA. Results for loci of 25 bp and longer do not change, and the
run time is the same (39 s on chr21). On the simulation `-S` raises loci found from 89.5% to
90.8%, with 1.20 calls per locus instead of 1.08. The thresholds (period 4, score 22, evidence
14) were chosen on this same chr21 catalog and random sequence, so the table is not a held-out
test. A sweep of uniform thresholds is in `bench/short_sweep.txt`.

### 8. Scaling

Per-chromosome runs on T2T-CHM13 (`-p 50000`, 2 workers) fit

    time  = 0.61 s per Mb of sequence + 4.0 s per Mb of annotated satellite   (R2 = 0.86)
    memory = 81 MB + 1.06 MB per Mb of the longest sequence                   (R2 = 0.98)

so run time grows linearly with genome size, plus a cost for the satellite arrays, which are
aligned base by base. One against two workers (the machine has only 2 cores, so nothing beyond
2 could be measured):

| Chromosome | 1 worker | 2 workers | Speed-up |
|---|---|---|---|
| chr21 (45 Mb) | 107.9 s | 56.2 s | 1.92 |
| chr18 (81 Mb) | 99.9 s | 49.2 s | 2.03 |
| chr10 (135 Mb) | 172.0 s | 89.3 s | 1.93 |
| chr9 (151 Mb, 32 Mb satellite) | 307.3 s | 165.5 s | 1.86 |

## Where tandex stands

- **Better than the others tested:** sensitivity and period accuracy on the simulation (all
  period and divergence classes), coverage of annotated human satellites on CHM13, long periods
  (above 2 kb, up to 45 kb rDNA units), fewest chance calls, least redundant output, low memory,
  and it is far faster than TRF and ULTRA on large satellite arrays.
- **Not better:** TRF agrees slightly more often with a curated human catalog's motifs (0.4 to
  3.5 points), as does tantan on loci of 100 bp or more; tantan is faster on human chromosomes in
  a single process (about 1.4 to 2 times). On repeats shorter than 25 bp tandex is level with
  ULTRA and tantan only with `-S`; the default setting finds fewer.
- **Untested:** TRF and ULTRA on CHM13 beyond chr21 (TRF timed out on chr22), TRASH (needs an R
  stack), ULTRA with period limits above 500 (needs more than 7 GB).
- The simulation is the author's own design, and parameters were tuned on a sibling simulation
  and on the chr21 catalog. Arabidopsis was not used for tuning.
- tandex annotates assembled sequence. It does not genotype repeats from reads (TRGT, LongTR,
  Straglr, ExpansionHunter, GangSTR, HipSTR do that).

## Comparisons still to run

Two comparisons need a larger machine than the one used here. Scripts are provided but were not
run by the author:

- `bench/run_trash.sh` runs TRASH (needs R 4.4.2 or newer) and converts its per-monomer output
  with `bench/trash2tsv.py`; the converter was checked on TRASH's bundled example output.
- `bench/run_ultra_long.sh` runs ULTRA with a long period limit. ULTRA's own estimate for
  `-p 2000 -t 2` is 16.9 GB of memory.

Both print the commands to score the result with the evaluation scripts used for the other tools.

## Limitations

- The E-value assumes random sequence with the chromosome's overall base composition; locally
  AT-rich or low-complexity regions make it optimistic.
- No alignment text output.
- `-t` is not used for sequences under 4 Mb, so assemblies of many small contigs run on one
  core.
- With `-t`, workers that start on the same large array at the same moment still each align it
  once; later workers reuse the published call.
- Arrays with fewer than about 3 copies at more than 20% divergence are often missed, as are
  microsatellites under about 20 bp.
- Compound repeats are reported as separate, possibly overlapping calls.

## Data sources

chr21: nf-core/test-datasets (GRCh38). Catalog: broadinstitute/tandem-repeat-catalog release
v1.0. Arabidopsis: schatzlab/Col-CEN v1.2. CHM13 v2.0 and censat v2.1: marbl/CHM13.

## Citing

If you use tandex, please cite the archived release:

> tandex 0.6.0. Zenodo. https://doi.org/10.5281/zenodo.24256789

## Licence

MIT.
