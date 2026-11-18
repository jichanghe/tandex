# Benchmarks

Scripts and saved results behind the numbers in the main README. Paths inside the shell scripts
point at the machine the benchmarks were run on (2 cores, 7 GB); edit them before running.

## Simulation with known truth

| File | What it does |
|---|---|
| `simulate.py` | Makes a random genome with implanted tandem repeats; writes FASTA and a truth table |
| `evaluate.py` | Scores a call set against the truth: base and locus recall, precision, period accuracy |
| `run_all.sh` | Runs tandex, TRF, ULTRA, mreps and pytrf on one FASTA and records time and memory |
| `run_extra.sh` | Same for tantan, longdust and SRF |
| `run_pytrf.py` | Wrapper that runs pytrf and writes calls in the common format |
| `tm.py` | Times one command and appends wall time and peak memory to a file |

## Human GRCh38 chr21

| File | What it does |
|---|---|
| `catalog_eval.py` | Recall and motif agreement against the Broad tandem repeat catalog v1.0 |
| `validate.py` | Independent check (unit against next unit, with matched controls) of calls TRF does not report |
| `conc.py` | Overlap between two call sets |
| `incompat.py` | Lists catalog loci where a tool's motif disagrees with the catalog |

## Arabidopsis Col-CEN

| File | What it does |
|---|---|
| `cen_mask.py` | Builds a tool-independent mask of the 178-bp centromere satellite by autocorrelation |
| `cen_eval.py` | Coverage of that mask by each tool |

## Human T2T-CHM13

| File | What it does |
|---|---|
| `censat_eval.py` | Coverage of the censat v2.1 satellite annotation, by class |
| `gap_check.py` | Tests whether annotated satellite that tandex leaves uncovered is actually periodic |
| `run_group.sh` | Runs tandex chromosome by chromosome on the whole genome |

## Seven ape genomes

| File | What it does |
|---|---|
| `run_apes.sh` | Downloads the six T2T ape assemblies from NCBI, runs tandex, writes per-genome summaries |
| `classify_periods.py` | Bases in repeats by period family, from a tandex TSV |
| `run_next.sh` | Lists ape annotation files on GenomeArk, measures scaling with workers, runs ULTRA (long periods) and TRASH on CHM13 chr21 |
| `ape_censat_eval.py` | Coverage of a T2T CenSat annotation by calls, pairing sequence names by length |
| `scaling_14threads.tsv` | Wall time with 1 to 12 workers on a 14-thread desktop |
| `apes_seven.tsv` | Results for human CHM13 and the six apes |

## Higher-order repeats (`-H`)

| File | What it does |
|---|---|
| `hor_check.py` | Compares the HOR unit tandex reports per array with an autocorrelation estimate |
| `hor_chm13.tsv` | Result for the 28 active HOR arrays of T2T-CHM13, tandex 0.8.0 `-H` |
| `run_hprc_pilot.sh`, `run_hprc_hor.sh` | Pilot on 20 HPRC release 2 haplotypes: download and run, then re-run with `-H` |
| `hprc_units.py`, `run_hprc_all.sh`, `hprc_all.tsv` | All 462 HPRC haplotypes: cut out the arrays, run tandex `-H`, unit by tandex and by self-comparison |
| `hprc_verify.py` | The self-comparison check for the 20-haplotype pilot |
| `hprc_cen_hor.py` | Per active centromere array: length and the HOR unit tandex reports |

## ULTRA with long periods and TRASH

Run on CHM13 chr21 by `run_next.sh` (part C) and `run_fix.sh`; results in `chr21_rivals.txt`.
`run_trash.sh` and `run_ultra_long.sh` are the earlier stand-alone drafts and were not used.

| File | What it does |
|---|---|
| `run_trash.sh`, `trash2tsv.py` | Run TRASH and convert its output for the evaluation scripts |
| `run_ultra_long.sh` | Run ULTRA with a period limit above 500 |

## Saved results

`*.eval.txt` and `*.out` hold evaluation output; `*.times` hold wall time and peak memory per tool.

| Files | Data set |
|---|---|
| `sim2.*` | Simulation |
| `chr21.*`, `catalog.out`, `validate.out` | GRCh38 chr21 |
| `cen1.*`, `ath.times` | Arabidopsis Col-CEN |
| `chr9`, `chr22`, `chrX`, `wg.*`, `chm13*.times` | T2T-CHM13 (`chm13_workers.times`: 1 against 2 workers, version 0.7.0) |
| `short_sweep.txt` | Short-repeat recall against chance calls for uniform thresholds |
