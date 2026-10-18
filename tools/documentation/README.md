# Documentation examples and illustrations

The README is an executable introduction to the public package API. Its code
is checked against an **installed** sparsepca package; generated help examples
are extracted from the source `man/` pages. This keeps documentation validation
separate from `pkgload`, which exposes development-only namespace objects.

## Verify the examples and HTML help

Install the intended candidate and its dependencies separately. From the
repository root, run:

```sh
Rscript --vanilla tools/documentation/check-examples.R . /tmp/sparsepca-documentation-run
```

The output directory must be new or empty. A fresh directory prevents stale
plots or logs from being mistaken for results of the current README. If the
candidate is in a separate R library, either set `R_LIBS_USER` or pass that
library as the third argument:

```sh
Rscript --vanilla tools/documentation/check-examples.R . /tmp/sparsepca-documentation-run /path/to/candidate-library
```

The optional library argument can contain several paths separated by the
platform's normal R library separator (`:` on Unix, `;` on Windows). The
runner requires the installed package version to match source DESCRIPTION.
It records the package location and installed code-file checksums; a matching
version string alone does not establish that two source archives are identical.
The runner never installs packages or downloads data.

The check:

- Parses all labeled README R chunks, recording their stable IDs and source
  line ranges. Duplicate IDs, missing labels, malformed markers, and unclosed
  fences fail validation.
- Runs README chunks in document order in one shared example environment.
  Chunks after a failure or unexpected warning are marked blocked, because
  they may depend on unfinished work.
- Runs every help topic in its own environment. An object created by one help
  example is not available to another topic.
- Exposes the installed package's public exports and normal attached R
  packages, without making the runner's helpers or package internals available
  as unqualified functions.
- Captures printed output, messages, warnings, errors, and plots. Any unexpected
  warning makes the run fail; a warning is not treated as a successful example.
- Renders all source help topics to HTML, including build/install-stage Rd
  macros and local alias links. Rendering warnings and errors also fail.
- Copies this R installation's `R.css` stylesheet and creates a local topic
  index. The explicitly referenced `stats::prcomp` page is rendered from its
  installed help source, without running its examples or fetching other pages.
- Verifies that the source files did not change during execution.

`example-results.csv` identifies each stage, its status, source lines, and any
explicit skip reason. `code/` contains the exact extracted R code; `logs/`,
`plots/`, and `html/` retain the associated output. PNG plots are used when the
R build supports them, otherwise PDF. `source-manifest.csv` and
`installed-package-manifest.csv` identify the source and installed files with
MD5 content checksums. `run-manifest.txt`, its RDS counterpart, and
`session-info.txt` record Git state, package/dependency versions, library paths,
execution flags, and timestamps. The metadata labels the hash algorithm; these
checksums identify inputs and are not the release archive's SHA-256 receipt.

Open `html/00Index.html` to browse the package's eight help topics. Package
links to `stats::prcomp` resolve to `html/stats-prcomp.html` within this preview.
Further help-topic references on that external page are rendered as plain text;
the runner does not copy or crawl the rest of R's help system. Only that exact
external reference is localized, using parsed Rd rather than rewriting HTML.
The logs record its installed package version and location. Mathematical
formatting still uses the R HTML converter's configured renderer and assets.

This is a source-documentation check against the identified installed package.
Release validation still builds a source archive and runs `R CMD check
--as-cran` on those exact bytes, including installed-package examples and the
complete PDF manual. HTML generation alone does not establish correct browser
layout; inspect rendered mathematics, links, figures, and the PDF separately.

## Inspect the actual GitHub rendering

GitHub's mathematics renderer can reject a macro even when a local TeX
renderer accepts it. The optional browser check verifies the public README
at an immutable commit, including the output inside each equation container.
It requires Node.js, `playwright-core`, and Chromium as development tools;
they are separate from the R package dependencies. Set `NODE_PATH` if the
Node module is installed outside the normal search path, and set
`SPARSEPCA_CHROMIUM_EXECUTABLE` to an existing browser executable.

First run the offline fixtures that exercise the previously missed macro
error, incomplete renderings, and the corrected soft-threshold formula:

```sh
node tools/documentation/test-github-rendering.cjs
```

After pushing a candidate, replace `COMMIT` below with its full 40-character
SHA. The output directory must be absent or empty.

```sh
git show COMMIT:README.md > /tmp/sparsepca-readme.md
node tools/documentation/check-github-rendering.cjs \
  COMMIT /tmp/sparsepca-readme.md /tmp/sparsepca-github-rendering
```

The check first requires the snapshot to match GitHub's raw source byte for
byte. It then checks rendered mathematics, loaded images, local anchors,
browser errors, and mobile overflow. The output contains the inspected
source identity, rendered HTML, screenshots, and structured findings.
Exit status 0 means these checks passed, 1 indicates findings, and 2 means
the inspection could not complete. Inspect the screenshots as well to judge
equation layout and figure readability.

The browser uses a fresh unauthenticated context and normally verifies
certificates. If an isolated public-browser environment cannot use its
proxy's certificate authority, the explicit environment setting
`SPARSEPCA_PUBLIC_BROWSER_IGNORE_HTTPS_ERRORS=1` permits that inspection
and records the exception in its results.

Keep formulas within GitHub's supported
[mathematical-expression syntax](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
The soft-threshold formula uses `\mathrm{soft}` and `\mathrm{sign}` for
plain function labels. The browser check detects GitHub's macro-rejection
messages as failures, even when they appear inside a normal equation wrapper.

## Mark README code explicitly

Place an HTML comment immediately before each fenced R block. These comments
are invisible in the rendered article. IDs are unique lowercase words and
numbers separated by hyphens; use meaningful names that survive section edits.

````markdown
<!-- example: iris-fit -->
```r
X <- as.matrix(iris[, 1:4])
fit <- spca(X, k = 2, scale = TRUE, verbose = FALSE)
summary(fit)
```
````

The default is `eval: true`. An installation block is intentionally left for
the reader to execute:

````markdown
<!-- example: install-development; eval: false -->
```r
remotes::install_github("erichson/spca")
```
````

For other deliberately unevaluated examples, provide a concrete reason:

````markdown
<!-- example: external-images; eval: false; reason: Requires separately obtained image files -->
```r
images <- readRDS("my-authorized-image-data.rds")
```
````

The `example`, `eval`, and `reason` fields are separated by semicolons, so do
not use semicolons inside a reason. An `install-` ID with `eval: false` receives
an explicit installation reason automatically; other unevaluated examples
must supply their own. Only ordinary `r` fences are executable, with either
backticks or tildes. Math and other language fences are not executed.

Skips are never inferred from a comment or string mentioning an installation
function. Accidentally executable installation/download calls are rejected
and reported as errors; they are not silently converted into skips. External
data illustrations need an honest execution boundary and provenance. A
self-contained simulation or built-in dataset should supply the runnable
teaching example when historical data are unavailable.

Execution checks complement mathematical validation. Application examples
should also establish the claims they illustrate: use identical preprocessing
for PCA comparisons; compare reconstruction or subspaces rather than vector
signs; calculate prediction from training preprocessing; assess robust models
in processed units; and use genuine compression when illustrating randomized
sketches. Do not assert timing thresholds, monotone support along every
nonconvex path, or general support-recovery guarantees from one simulation.

## Iris, sensors, and face images

The [feature-count example](target-features.md) uses `spca_target()` to
search for a representation based on three of the four iris measurements.
Its script saves every evaluated penalty, the chosen coefficients, and
session information, then draws the figure from that saved snapshot.

The [application reproduction and provenance guide](application-examples.md)
documents the iris and sensor examples and their retained data and numerical
results. Its figure script can redraw those figures from saved results with
`--plot-only`, without refitting. The [face reproduction guide](eigenfaces.md)
provides the corresponding workflow for the full cropped Yale face dataset,
with saved coefficients and error maps for all five application figures.

## Retained loading illustration

The existing loading illustration uses simulated data generated entirely by
`illustrate-loadings.R`. Twelve variables form three groups, with four noisy
measurements per latent signal. The latent signals are correlated, so ordinary
PCA mixes the groups. The illustration compares three-component ordinary PCA
with elastic net and hard-threshold sparse fits on the same centered and scaled
data.

Run from the repository root, with `pkgload` installed:

```sh
Rscript tools/documentation/illustrate-loadings.R
```

This existing illustration script regenerates both the fits and figures. It
saves raw loading coefficients, model parameters, nonzero counts, relative
Frobenius reconstruction errors, convergence diagnostics, and R session
information under `results/`. Its SVG and PNG figures are in `man/figures/`.
Signs are flipped for display only so each component's largest-magnitude
coefficient is positive; sparsity and reconstruction do not change. Original
signed coefficients are retained in CSV. Colors use a common coefficient
scale, with white for zero.

The different penalties illustrate interpretation and reconstruction
tradeoffs. They are not mathematically equivalent penalty strengths or a speed
benchmark. All data in this illustration are synthetic.

The benchmark figures and raw performance measurements are maintained
separately under `tools/benchmarks/` and `inst/figures/`. Documentation example
checks do not rerun or overwrite those experiments or plots. README-only
figure assets and these development scripts are excluded from the CRAN source
archive. Authoritative help text remains in `R/documentation.R`.
