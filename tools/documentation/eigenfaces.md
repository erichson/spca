# Reproduce the face application

The README face figures use all 2,410 cropped images in the prepared
Extended Yale Face Database B matrix. Each image has 96 rows and 84
columns. The fitted observations are the rows of `t(as.matrix(faces))`,
giving a 2,410 by 8,064 matrix. Every fit centers each pixel and divides
it by that pixel's sample standard deviation. The resulting components
describe variation in standardized pixel intensity.

The data preparation and the four displayed parameter settings follow
the package's [original application](https://github.com/erichson/spca/blob/fa896289a68801a4798f53948ff6c2b8f7c54e1c/README.md).
The figures are recomputed through the public APIs using the current
implementation. An additional matched nonrobust fit makes the robust
comparison easier to interpret.

## Obtain the input and run the fits

Download the prepared R data file to a directory outside the repository.
The immutable source URL is:

```text
https://raw.githubusercontent.com/erichson/data/fca24620977d21f5230ca9cc7a6107b25847d6fa/R/faces.RData
```

The expected file is 17,539,054 bytes. Its SHA-256 is:

```text
82674d0785f38c4c57eb688992eaf0fb2453e8c0409ef1f25b40da8de6a979a5
```

Its Git blob is `884cf6fa86db4f26adbeb22eeb019101f99585c8`. The script
checks the SHA-256, object name, dimensions, finite values, and intensity
range before fitting. The R file contains one object, `faces`, a data
frame with pixels in rows and images in columns. Values range from 0 to
255. No source photographs are copied into the repository or CRAN archive.

Run from the repository root with the development packages `pkgload` and
`digest` installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tools/documentation/eigenfaces.R --fit \
  /path/to/faces.RData /path/to/eigenface-checkpoints
```

The first path names the downloaded input. The second names a directory
outside the repository for intermediate coefficient/transform/score checkpoints.
Fitting uses `pkgload::load_all()` so it evaluates the current checkout.
The script records the package version, Git revision, hashes of executable
R files, dependency versions, BLAS/LAPACK, seed, and thread settings.

The full input occupies approximately 149 MiB as a dense double matrix.
Robust fitting also needs several dense work arrays. These applications
are development experiments and are excluded from routine unit tests,
help examples, and the CRAN source archive.

Individual fits can be executed and checkpointed separately:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tools/documentation/eigenfaces.R --fit-model robust \
  /path/to/faces.RData /path/to/eigenface-checkpoints
```

Supported model names are `ordinary`, `sparse`, `stronger`, `robust`,
and `robust_control`. After all five checkpoints exist, assemble the
compact results and figures:

```sh
Rscript tools/documentation/eigenfaces.R --assemble /path/to/eigenface-checkpoints
```

## What the comparisons measure

Every model has 25 components. Both randomized calls use `o = 20`,
`q = 2`, and `set.seed(42)`, so they use the same random starting draws.
The iterative fits use the original example's `tol = 1e-4` and
`max_iter = 1000`.

| Model | Public routine | `alpha` | `beta` | `gamma` |
| --- | --- | ---: | ---: | ---: |
| Ordinary | `stats::prcomp` | 0 | 0 | — |
| Sparse | `rspca` | 0.0001 | 0.1 | — |
| Stronger penalties | `rspca` | 0.0002 | 0.2 | — |
| Robust | `robspca` | 0.0001 | 0.01 | 0.9 |
| Matched nonrobust control | `spca` | 0.0001 | 0.01 | — |

The ordinary reference uses `rank. = 25`, `center = TRUE`, and
`scale. = TRUE`. Its loading matrix is `rotation`. For the other methods,
the displayed coefficients are the returned `loadings` matrix, B.
SPCA's iterative tolerance does not apply to the direct PCA calculation.

The first sparse comparison shows how these regularized settings reduce
the support of the coefficient maps. The stronger setting increases
both penalties. It therefore measures their joint effect. It does not
isolate a change in `alpha` alone. The robust setting also changes
`beta` relative to those randomized examples. The additional deterministic
nonrobust control uses the same `alpha`, `beta`, rank, data, preprocessing,
and stopping tolerance as the robust fit. This is the comparison that
isolates the effect of adding the sparse error term at the stated settings.

The script evaluates full-data reconstruction error for every fit:

```r
norm(processed - tcrossprod(fit$scores, fit$transform), "F") /
  norm(processed, "F")
```

For `prcomp`, the corresponding factors are `x` and `rotation`.
The implementation evaluates this formula in 100-row blocks to bound
temporary memory. It also evaluates the robust complete-model residual,
including the fitted sparse matrix S. A smaller error after adding S
describes how well the larger model represents the observed images; it
does not measure recovery of an independently known clean face. The data
contain no supplied pixelwise corruption labels.

Recorded nonzero counts use exact zeros in B and S. `active_components`
counts nonzero columns of B. The script verifies objective descent up to
floating-point roundoff and checks the final deterministic and robust
objectives by independently accumulating their loss and penalty terms.
It also verifies the robust soft-threshold update against the full-data
residual in every row block. The reported convergence flag is the solver's
relative objective stopping criterion.

The saved elapsed times record one run per model. They are execution
records for reproducing the application, not a timing comparison. Input
loading and the additional diagnostic calculations are outside those
times. The existing performance study remains in `tools/benchmarks/`.

### How stopping tolerance affects selected pixels

The sparse setting also has a separate completed run at `tol = 1e-5`.
With the data, penalties, seed, and implementation held fixed, the
nonzero coefficient fraction changes from 51.2183% at `tol = 1e-4`
to 40.3681% at `tol = 1e-5`. The full-data relative reconstruction
errors are 0.4732872 and 0.4724738, respectively. Similar reconstruction
errors can therefore accompany appreciably different support counts.
If conclusions depend on which pixels are selected, assess support
stability as well as the change in objective.

The main figures use the original tolerance, `1e-4`. The stricter sparse
fit is retained separately in `results/eigenfaces-sparse-stricter.rds`,
with its own parameters and objective history. The histories agree
through the shorter run; only the stopping tolerance differs. Data,
penalties, seed, runtime source hashes, dependency session, and all
other controls were checked for equality. To reproduce this additional
comparison in a separate checkpoint directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tools/documentation/eigenfaces.R --fit-model sparse --tol=1e-5 \
  /path/to/faces.RData /path/to/stricter-checkpoint
```

To retain that sensitivity fit alongside the five-model application and
write its comparison table, run the following R code from the repository
root after replacing the two checkpoint paths:

```r
regular <- readRDS("/path/to/eigenface-checkpoints/sparse.rds")
stricter <- readRDS("/path/to/stricter-checkpoint/sparse.rds")
stopifnot(identical(
  stricter$objective[seq_along(regular$objective)], regular$objective
))
saveRDS(stricter, "tools/documentation/results/eigenfaces-sparse-stricter.rds",
        compress = "xz", version = 3)
comparison <- rbind(
  cbind(tol = regular$metadata$common$tol, regular$metrics),
  cbind(tol = stricter$metadata$common$tol, stricter$metrics)
)
write.csv(comparison,
          "tools/documentation/results/eigenfaces-tolerance-comparison.csv",
          row.names = FALSE)
```

## Read and redraw the figures

For display, each coefficient vector is divided by its largest absolute
entry. This preserves its zero entries and makes its spatial pattern
visible across components with different amounts of ridge shrinkage.
The four coefficient panels share one symmetric scale from -1 to 1.
Color shows coefficient size relative to the largest magnitude in that
component. Raw fitted coefficients, their original norms, and their
maximum magnitudes are saved.
White denotes zero, blue negative values, and red positive values.
The displayed support fraction is the proportion of pixels with a
nonzero coefficient in that component.

Each loading column is oriented so its largest-magnitude entry is
positive. These sign changes affect only the display; original factors
remain in the saved checkpoints and compact result file. Components
keep the solver's order. Matching panel numbers across two fits does not
establish that their components have identical meanings.

The input vectorization puts 84 horizontal pixels before the next image
row. The display reshapes a vector with `matrix(v, nrow = 84, ncol = 96)`
and reverses its second index in `image()`, keeping the faces upright.
This convention is specific to the supplied matrix; an image flattened
by R's column-major `as.vector()` requires its corresponding inverse.

The sparse-error panel contains the first 25 rows of S, reshaped with
that same orientation. The row choice follows the input order, without
selection by fit quality or residual magnitude. All 25 errors share a
single symmetric color scale in standardized pixel units. Their signs
are preserved. These are fitted additive residuals for individual
observations, whereas the other panels contain component coefficients.

The files in `tools/documentation/results/` have the `eigenfaces-` prefix:

| File | Contents |
| --- | --- |
| `eigenfaces-data.rds` | Exact raw loadings, transforms, scores, means, scales, objective histories, and display signs for all five fits; selected sparse-error rows, numerical summaries, settings, and provenance. No input photographs or full dense S. |
| `eigenfaces-summary.csv` | Nonzero counts, full-data residual errors, objectives, convergence, and single-run elapsed times for all five models. |
| `eigenfaces-components.csv` | Per-component support counts, loading norms, maximum magnitudes, and display signs for all five fits. |
| `eigenfaces-objectives.csv` | Every recorded objective value for the iterative fits. The randomized and full-data objectives concern different optimization problems. |
| `eigenfaces-session-info.txt` | R, packages, platform, locale, BLAS/LAPACK, and thread settings. |
| `eigenfaces-verification.txt` | Output of the independent saved-factor verification. |
| `eigenfaces-figure-checks.json` | SHA-256 values for the ten final figures and the result of redrawing them from the saved factors. |
| `eigenfaces-sparse-stricter.rds` | Separate completed sparse fit at `tol = 1e-5`, including its factors, objective history, and original generation metadata. |
| `eigenfaces-tolerance-comparison.csv` | The two sparse fits' measurements, labeled by their stopping tolerance. |

To redraw all five PNGs and SVGs without accessing the source images,
loading the package, or fitting a model:

```sh
Rscript tools/documentation/eigenfaces.R --plot-only
```

To verify the reported reconstruction errors, support counts, factor
relationships, and objectives directly from the saved factors and input:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tools/documentation/check-eigenfaces.R /path/to/faces.RData
```

This check uses the recorded `rsvd` dependency and seed to reconstruct
the 45-row randomized sketch. It does not rerun any optimization.
The robust S entries are recovered by soft thresholding the residual,
then compared with the saved aggregate measurements and selected rows.
The separate stricter sparse fit is also verified when its saved RDS file
is present, including the common prefix of the two objective histories.

The figure files are `man/figures/eigenfaces-ordinary`,
`eigenfaces-sparse`, `eigenfaces-stronger`, `eigenfaces-robust`, and
`eigenfaces-errors`, each with `.png` and `.svg` extensions. PNGs use
1,600 by 1,904 pixels. SVGs retain vector text and embed the coefficient
maps as raster images, avoiding millions of small rectangles.

## Dataset acknowledgement

This experiment uses the **Extended Yale Face Database B** and acknowledges
the **Yale Face Database B**. Dataset provenance is documented by the
[UCSD dataset listing](https://vision.ucsd.edu/datasets/extended-yale-face-database-b-b)
and the [Yale dataset description and research-use terms](https://www.cs.yale.edu/cvc/projects/yalefacesB/yalefacesB.html).
The prepared cropped matrix was deposited by N. Benjamin Erichson in
[`erichson/data`](https://github.com/erichson/data/commit/fca24620977d21f5230ca9cc7a6107b25847d6fa).

The Yale database description requests acknowledgement and citation of
A. S. Georghiades, P. N. Belhumeur, and D. J. Kriegman, "From Few to Many:
Illumination Cone Models for Face Recognition under Variable Lighting and
Pose," *IEEE Transactions on Pattern Analysis and Machine Intelligence*
23(6), 643–660 (2001), [doi:10.1109/34.927464](https://doi.org/10.1109/34.927464).
The source images retain their dataset terms. The package software license
does not replace those terms.
