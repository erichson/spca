# Reproducing the release checks

The scripts in this directory are development tools and are excluded from the
CRAN source archive. Their output is evidence from the selected versions and
environment; a completed run should be retained with the release archive.

## Build and check one preserved archive

Commit the reviewed source, then build to a new directory outside the checkout:

```sh
R_BUILD=/path/to/R bash tools/release/build-source.sh /path/to/release-candidate
Rscript tools/release/check-source.R \
  /path/to/release-candidate/sparsepca_1.0.0.tar.gz /path/to/release-check
```

The build script refuses to overwrite an existing archive or build an uncommitted
checkout. It saves source commit/tree identifiers, the build environment and log,
an archive listing, SHA-256, and per-file hashes. The check script verifies the
archive bytes before and after a complete `--as-cran` check and retains all
diagnostics. It requires `rcmdcheck`, `testthat`, `V8`, HTML Tidy, and a working
TeX installation for the manual. Run it in a fresh process for each R version,
selecting suitable dependency libraries and single-thread BLAS settings before
starting R. No NOTE is hidden; errors and warnings return a failure status.

## Baseline and candidate libraries

Install the published baseline and the candidate into separate R libraries.
Start a fresh R process for each comparison, putting the selected library before
the shared dependency library in `R_LIBS_USER`. Do not replace an already loaded
namespace within a process.

```sh
R CMD INSTALL --library=baseline-library sparsepca_0.1.2.tar.gz
R CMD INSTALL --library=candidate-library sparsepca_1.0.0.tar.gz
```

Run the commands with the required dependency library also available. Record
the source archive SHA-256 hashes, the candidate source commit and tree, R and
package versions, operating system, and BLAS/LAPACK implementation. For the
release candidate, check the archive produced by `R CMD build`, and retain those
same bytes for any subsequent checks or upload.

## Downstream archives

Download the current source archives for the four direct CRAN imports:
`categoryEncodings`, `clusterability`, `robCompositions`, and `SparseBiplots`.
Install their required and suggested dependencies and the development tool
`rcmdcheck`. Run the following once for each selected sparsepca library:

```sh
R_LIBS_USER=baseline-library:dependencies-library \
  Rscript tools/release/check-downstream.R archives baseline-checks 0.1.2
R_LIBS_USER=candidate-library:dependencies-library \
  Rscript tools/release/check-downstream.R archives candidate-checks 1.0.0
```

These downstream checks run installation, loading, code analysis, examples, and
the downstream packages' test suites. They use `--no-manual`: PDF typesetting of
unchanged downstream documentation is outside this compatibility check. Incoming
submission checks are disabled because the downstream packages are already
published. Suggested packages are required, and the packages' own CRAN skips
are respected. All errors, warnings, NOTES, check output, installed package
versions, and settings are retained. Compare the two `summary.csv` files and
inspect individual diagnostics rather than interpreting a baseline failure as
a regression in the candidate.

## Focused integrations

Install `parameters` and Bioconductor's `scPCA` as well as the four CRAN packages.
The integration script exercises the public entry points that use sparsepca,
including the robust branch of `parameters::principal_components()` and both
variable-projection branches of `scPCA::scPCA()`. It also checks the score-only
path inside `clusterability` with all four centering/scaling combinations.

```sh
R_LIBS_USER=baseline-library:dependencies-library \
  Rscript tools/release/integration-probes.R baseline-probes 0.1.2
R_LIBS_USER=candidate-library:dependencies-library \
  Rscript tools/release/integration-probes.R candidate-probes 1.0.0
```

On Windows, separate library paths with `;` instead of `:`. The scripts use
small, seeded synthetic datasets. Results, warnings, errors, versions, and full
returned objects are saved. Numerical values can change when an objective or
variance calculation is corrected; the probes check that real consumer code
continues to run and receives the documented structures. They do not replace a
complete check of the large optional integration package `parameters` or the
complete Bioconductor dependency ecosystem.

## Continuous integration

`R-CMD-check.yaml` covers Linux at R 4.0.0, R 4.4.0, current R release, and R-devel,
as well as current R on Windows and macOS. Its R 4.0.0 job uses Matrix 1.3-0 and
rsvd 1.0.5, with compatible test tools from a dated CRAN snapshot. Other jobs use
current compatible dependencies. These routine checks explicitly omit the PDF
manual and preserve all remaining diagnostics.

`CRAN-source-check.yaml` builds one archive under current R release, records its
SHA-256 and source identity, and checks the identical archive under release and
R-devel. It requires the declared suggested packages, runs incoming checks,
installs HTML and math validators, and builds the full PDF manual. It saves the
archive, logs, environment, manual, and result objects. Errors and warnings fail
the workflow; NOTES remain visible and must be assessed before submission.

## Policy references

- [CRAN repository policy](https://cran.r-project.org/web/packages/policies.html)
- [CRAN submission checklist](https://cran.r-project.org/web/packages/submission_checklist.html)
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-devel/R-exts.html)
- [Current sparsepca reverse dependencies](https://cran.r-project.org/package=sparsepca)
- [Bioconductor scPCA](https://bioconductor.org/packages/scPCA/)

Before submission, check the current CRAN package/check pages and public incoming
queues again. A clear preflight result does not replace CRAN's incoming review
or the maintainer's confirmation of the submission.
