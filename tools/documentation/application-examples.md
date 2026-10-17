# Application examples and figure provenance

The README combines two freshly computed applications with the original
face illustrations. The fresh examples require no external dataset
download. The historical pictures are retained illustrations, not new
version 1.0.0 fits.

## Reproduce the new figures

From the repository root, with the development dependency `pkgload` installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tools/documentation/application-examples.R
```

This loads the current source checkout, fits both examples through the public
`spca()`, `robspca()`, and `stats::prcomp()` APIs, saves the complete data and
fitted objects, and writes these figures:

- `man/figures/iris-application.png` and `.svg`;
- `man/figures/robust-sensors.png` and `.svg`.

To redraw the figures from the recorded data without loading the package,
refitting, or changing the saved measurements:

```sh
Rscript tools/documentation/application-examples.R --plot-only
```

Both commands accept the repository path as an optional final argument.
The plotting code uses base R graphics. RDS files and CSV files in `results/`
have the `application-` prefix:

| File | Contents |
| --- | --- |
| `application-data.rds` | Complete inputs, fits, reconstructions, display coordinates, parameters, generation revision, RNG kind, and thread settings. |
| `application-summary.csv` | Numerical results in application/method/metric/value form. `objective_converged = 1` means the stated objective stopping criterion was met. |
| `application-iris-scores.csv` | Original and displayed scores, observation indices, and species labels. |
| `application-iris-loadings.csv` | Original and displayed loading coefficients. |
| `application-sensor-traces.csv` | Every sensor observation, known clean value, fitted reconstruction, injected error, and estimated sparse error. |
| `application-sensor-corruptions.csv` | All locations with an injected error or a nonzero fitted sparse error, including the two additional small estimates. |
| `application-session-info.txt` | R, package, BLAS/LAPACK, platform, locale, and thread environment. |

These are small explanatory experiments, not performance benchmarks.
The README's existing performance measurements and figures are maintained
separately. All development scripts and README-only assets remain excluded
from the CRAN source archive by `.Rbuildignore`.

## Iris: interpret the score map through its coefficients

The example uses the 150 flowers and four measurements in R's
[`datasets::iris`](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/iris.html).
The dataset documentation credits Anderson's measurements and Fisher's
published data. Species labels color the figure only; neither fitting
call receives them.

Both methods center each measurement and divide it by its sample standard
deviation. The ordinary reference is `prcomp(..., rank. = 2, center = TRUE,
scale. = TRUE)`. The sparse fit uses two components, `alpha = 0.01`,
`beta = 0.001`, `tol = 1e-8`, and `max_iter = 3000`.

The sparse fit has three nonzero coefficients in its first component and
one in its second, compared with four in each ordinary component. The first
sparse score combines sepal length, petal length, and petal width; the second
uses sepal width. This is a descriptive projection of the observed data,
not a trained species classifier or an estimate of predictive accuracy.

Each displayed component is oriented so its largest-magnitude coefficient
is positive. The same sign change is applied to its scores. The original
signs remain in the fit snapshot and CSV files. The numerical coefficient
panels are rounded to three decimal places; exact zeros are shown as `0`.
The two score plots use common axis limits, and both coefficient panels
use the same color scale. The sparse coefficients are the returned matrix
`B`; they are not additionally normalized, and they are not interpreted as
variable-score correlations.

The recorded reconstruction errors use the same centered and scaled input
for both fits: `norm(reconstruction - processed, "F") / norm(processed, "F")`.
The retained ordinary factors give 0.2046165 and the sparse factors give
0.2108651. The comparison uses `ordinary$x %*% t(ordinary$rotation)` for PCA
and `sparse$scores %*% t(sparse$transform)` for SPCA. The extra error
diagnostics were calculated from the saved fits without refitting them.

## Sensors: separate large entrywise errors from shared variation

The script generates 120 observations of eight synthetic sensors using
`set.seed(17)`. Sensors 1--4 share a sine signal, and sensors 5--8 share a
cosine signal over two periods. Independent Gaussian measurement noise
has standard deviation 0.08. Four individual entries receive additive
errors of `+7`, `-7`, `+8`, and `-8` at the recorded sample/sensor pairs.

Both sparse fits use `k = 2`, `alpha = 0.005`, `beta = 0.001`,
`center = FALSE`, `scale = FALSE`, `tol = 1e-8`, and `max_iter = 3000`.
The robust fit additionally uses `gamma = 0.3`. Centering and scaling are
disabled so the generating signal, the reconstruction, the sparse errors,
and `gamma` retain the same sensor units. The known truth is available
because the data are simulated; it is not supplied to either optimizer.

The reconstruction is `tcrossprod(fit$scores, fit$transform)`, namely
`X B A'` here. In the robust model, `X` in those scores is still the observed
input. The separate `fit$sparse` object is `S`; it is not itself a loading
matrix or the reconstructed signal.

For this fixed simulation, relative clean-signal error is
`norm(reconstruction - truth, "F") / norm(truth, "F")`, measured across all
120 x 8 entries. The ordinary sparse fit gives 0.42500559, and the robust
fit gives 0.09090952. The robust fit identifies all four injected locations
with nonzero `S` and produces two other small entries of approximately
-0.02135 and -0.02144. The ordinary, robust, and iris sparse fits meet their
objective stopping criterion after 102, 622, and 296 iterations, respectively.
These results illustrate the selected data and penalties; they are not
an outlier detection guarantee or a proof of global optimization.

The top sensor panel shows the observed readings and clean truth for
Sensor 1. The middle panel shows the two reconstructions using an expanded,
explicitly labeled vertical scale. In particular, a corruption in Sensor 3
at sample 75 affects the ordinary reconstruction of Sensor 1, illustrating
how a bad entry can influence a shared component. The final panel shows
`S` across all eight sensors and the known injected errors. Its two very
small additional estimates occur at the same sample index, 34, in different
sensors; their markers overlap at this display scale. No error values are
removed from the recorded data.

## Historical face illustrations

The original GitHub walkthrough is preserved at
[`fa896289a68801a4798f53948ff6c2b8f7c54e1c`](https://github.com/erichson/spca/blob/fa896289a68801a4798f53948ff6c2b8f7c54e1c/README.md)
(April 7, 2018). The first three panels were added in
[`b1e34c0d79e8dec2d4c2538c923fe48b495fc025`](https://github.com/erichson/spca/commit/b1e34c0d79e8dec2d4c2538c923fe48b495fc025),
and the robust loading/error panels in
[`0c75adf458baa6fe3762fe244232947b0324dcdb`](https://github.com/erichson/spca/commit/0c75adf458baa6fe3762fe244232947b0324dcdb)
(both March 26, 2018). All five image files are unchanged from that history.

| Historical asset | Meaning in the original walkthrough |
| --- | --- |
| `plots/eigenfaces.png` | Twenty-five ordinary loading coefficient maps from a fit with `alpha = beta = 0`. |
| `plots/sparseeigenfaces.png` | Twenty-five randomized sparse loading maps, `alpha = 1e-4`, `beta = 1e-1`. |
| `plots/sparsereigenfaces.png` | Stronger-penalty randomized sparse maps, `alpha = 2e-4`, `beta = 2e-1`. Both penalties changed. |
| `plots/robustsparseeigenfaces.png` | Robust sparse loading maps, `alpha = 1e-4`, `beta = 1e-2`, `gamma = 0.9`. |
| `plots/sparsecomp.png` | Sparse errors in the first 25 observed faces: rows of `S`, reshaped into images. |

The original walkthrough describes 2,410 cropped/downsampled faces of 38
people from the Extended Yale Face Database B, with 96 x 84 pixels per
image. It uses `t(faces)` so observations are rows and pixel positions are
columns. Those dimensions are recorded from the historical explanation;
the current documentation refresh did not download or refit the full face
data. The face input is dense, even though the fitted loadings can be sparse.

Individual historical `image()` calls choose their own color range.
Compare the localization of coefficients and errors, not the absolute
gray level across components or panels. The robust images were created
before the objective correction in 1.0.0; no claim is made that they show
newly converged 1.0.0 results. Historical variance percentages and timings
are not used as current numerical evidence.

The original data link points to `erichson/data/R/faces.RData`. The most
recent modifying file commit is
[`fca24620977d21f5230ca9cc7a6107b25847d6fa`](https://github.com/erichson/data/commit/fca24620977d21f5230ca9cc7a6107b25847d6fa)
(August 23, 2016, "cropped yale faces"). The immutable raw URL is:

```text
https://raw.githubusercontent.com/erichson/data/fca24620977d21f5230ca9cc7a6107b25847d6fa/R/faces.RData
```

The file metadata reports 17,539,054 bytes and Git blob
`884cf6fa86db4f26adbeb22eeb019101f99585c8`. No new copies of the source face
photographs are included in this documentation change.

Dataset provenance remains with the
[Extended Yale Face Database B](https://vision.ucsd.edu/datasets/extended-yale-face-database-b-b)
and [Yale Face Database B](https://www.cs.yale.edu/cvc/projects/yalefacesB/yalefacesB.html).
The original Yale page describes research use, acknowledgement/citation
requirements, and restrictions on incorporating source images into a larger
publicly distributed database. Retaining the existing research illustrations
does not grant additional rights to the source photographs through the
package's software license. The old UCSD archive link is unavailable; the
current dataset listing is used for provenance instead.

The published 0.1.2 CRAN source archive contains the three function help
pages, source, and tests, but no README, face data, or image assets. Its
robust help example is a Gaussian factor simulation without explicit gross
contamination. The new sensor application makes that distinction visible.
