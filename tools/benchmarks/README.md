# Reproducible performance experiments

These experiments answer three separate questions: whether a randomized sketch
accelerates the package's penalized sparse PCA objective; how the unpenalized
special case compares with a standard external PCA implementation; and how
preserving sparse input storage changes memory use. They also include a flat
spectrum where a small sketch captures less of the leading subspace.

## Reproduce

Run from the repository root on Linux with R, `rsvd`, `Matrix`, Python 3, NumPy,
Matplotlib, and GNU `time` installed:

```sh
python3 tools/benchmarks/run.py --data-dir /path/to/temporary/benchmark-inputs
python3 tools/benchmarks/summarize.py
```

Set `R_LIBS_USER` if dependencies are in a nondefault library. Pass `--rscript`
to select an R version. `--cases spca-small --repetitions 1` provides a pilot;
save pilots to a separate `--results-dir` so they cannot be mixed with the
three-repetition study. Input files are reproducible and large, so keep the
chosen data directory outside the Git repository.

The harness snapshots the current R sources for the duration of the experiment,
records their SHA-256 hashes, and refuses to mix source revisions in one results
directory. It runs one fresh R process at a time and alternates method order
between repetitions. Each completed run has a CSV measurement, stdout/stderr log,
and GNU time report. The combined `results/measurements.csv` is the raw dataset;
`results/summary.csv`, `results/tables.md`, and the figures are derived from it.
Changing a caption or plot never requires repeating the numerical experiment.

The first small-case measurements overlapped a short package test run. Those
six measurements are preserved in `results/superseded-small/` and excluded
from the summary. The complete small-case group was then repeated using the
same frozen source, input, parameters, seeds, and run order after the remaining
workloads finished. No other measurements were discarded.

## Workloads and comparisons

All workloads request ten components, center the input, and retain its original
column scales. Data generation uses seed 1729. Independent fit seeds 5001, 5002,
and 5003 are shared by the methods within each repetition.

| Case | Matrix dimensions | Data and comparison |
| --- | --- | --- |
| `spca-small` | 2,000 × 600 | Ten disjoint loading groups plus Gaussian noise; `spca` versus `rspca`. |
| `spca-medium` | 5,000 × 1,500 | The same generating model at a larger size. |
| `spca-large` | 10,000 × 3,000 | The same generating model at a larger size. |
| `pca-lowrank` | 10,000 × 3,000 | Ordinary PCA on the approximately low-rank input: `stats::prcomp` versus unpenalized `rspca`. |
| `pca-flat` | 2,000 × 600 | Independent Gaussian entries, with no small set of dominant components: `prcomp` versus unpenalized `rspca`. |
| `sparse-input` | 10,000 × 3,000 | Ten nonoverlapping observation/feature groups with small noise in stored entries; 10% stored entries. Compare `rspca(as.matrix(X))` with `rspca(X)`. |

The regularized fits use `alpha = 0.001`, `beta = 0.0001`, the elastic-net
penalty, `tol = 1e-5`, and `max_iter = 500`. Randomized fits use `o = 20` and
`q = 2`. These are relative penalties: the deterministic solver scales them by
the squared spectral norm of the preprocessed input, and the randomized solver
uses the squared norm of its sketch. In addition to each fit's own full-input
objective, the summary evaluates both methods with the deterministic spectral
scale. This prevents small differences in penalty scaling from being mistaken
for differences in objective quality.

The external comparison uses **zero penalties for both methods** and
`prcomp(rank. = 10, center = TRUE, scale. = FALSE)`. Both return loadings and
scores. The reference is an established deterministic implementation of the
same PCA problem. R's `prcomp` calls the ordinary SVD and then returns the
requested number of vectors; `rank.` does not turn it into a randomized or
iterative truncated SVD. Its additional factorization work is relevant to
interpreting the comparison. See the official
[prcomp documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html)
and [R source](https://svn.r-project.org/R/trunk/src/library/stats/R/prcomp.R).

The SPCA-versus-RSPCA comparison measures acceleration within this package. It
does not establish superiority over every sparse PCA algorithm. In particular,
these experiments do not treat the normalized loadings returned by `elasticnet`
as identical to this package's unnormalized penalized coefficients.

## Timing and memory definitions

Elapsed time surrounds the public fitting call, including preprocessing,
factorization or sketch construction, optimization, and result construction.
Input generation, input-file loading, package loading, a small 80 × 30 warm-up,
and post-fit quality calculations are excluded. The dense-conversion method in
the sparse-input experiment includes `as.matrix(X)` inside the timer. Garbage
collection is requested before the timer, and garbage collection that occurs
inside the fitting call remains included. Filesystem caches are not flushed.

All BLAS/OpenMP thread environment variables are set to one before starting R.
The supplied results were collected in a shared virtualized environment with
other package builds paused during measurements. They are measurements of this
machine and software stack, not hardware-independent performance guarantees.

The primary memory measure, `peak_rss_mib`, reads Linux's `VmHWM` immediately
after fitting. It is the **whole R process resident high-water mark through fit
completion**, in MiB (2^20 bytes). It includes the R runtime, loaded namespaces,
input loading and storage, warm-up, preprocessing, temporary allocations, and
the fitted result. It excludes any later memory peak due to quality validation.
It is neither incremental memory, cumulative allocation, nor output size.
Every worker preloads `rsvd` and `Matrix`, including the `prcomp` reference,
so the comparisons share the same namespace setup.

GNU time independently records `process_peak_rss_mib` for the complete worker,
including post-fit validation. Quality is computed in bounded row blocks so
that checking sparse input does not allocate a full dense reference matrix.
The two resident-memory measures may differ and are both retained. The sparse
comparison begins with the same compressed-column sparse input in both workers;
the dense branch explicitly converts that input as part of the fitted operation.

## Quality and aggregation

Reconstruction accuracy is `norm(Xc - scores %*% t(transform), "F") /
norm(Xc, "F")`, where `Xc` is the identically centered input. This is an
actual full-input residual, including the part omitted by a randomized sketch.
The raw data also record fitted energy, adjusted score energy, loading
nonzeros, effective penalties, iteration counts, convergence status, warnings,
and the full reconstruction-plus-penalty objective.
The separately saved `logged_objective` for `rspca` is evaluated on its sketch;
it must not be compared directly with the full-input `spca` objective. The
full-input objective columns and the common-penalty summary provide that
comparison.

The input contains noise, so even the best rank-ten fit has a nonzero residual.
On the flat-spectrum case, a relative reconstruction error near one is expected:
ten components explain little of a matrix with hundreds of similarly sized
singular values. Compare captured/fitted energy as well as residual error there;
a small difference between two large residuals can conceal a meaningful loss
of leading-subspace quality.

Reported elapsed times and resident memory are medians of three isolated runs.
Time error bars show the observed minimum and maximum, not confidence intervals.
Speedup is the **reference median elapsed time divided by the candidate median**.
Its error bars show the observed range of paired reference/candidate ratios.
No fragile wall-clock assertions are part of the package's unit tests.

## Provenance

The supplied environment is R 4.6.1 on Ubuntu 24.04.3 LTS, `rsvd` 1.0.6 and
`Matrix` 1.7-6, with OpenBLAS 0.3.26 and LAPACK 3.12.0. The virtual machine
exposes nine Intel Xeon Platinum 8573C CPUs and an 8 GiB cgroup memory limit;
each numerical fit uses one BLAS thread. Exact session, hardware, thread,
source-hash, and process reports are saved under `results/`.

These runs execute a frozen snapshot of candidate R source. The release review
must compare the benchmarked implementations with the final source and explain
any transferred evidence. They should not be described as a fresh benchmark of
a later archive merely because its version number is the same.
