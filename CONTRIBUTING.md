# Contributing to sparsepca

Bug reports, reproducible numerical examples, and focused improvements are
welcome. Open an [issue](https://github.com/erichson/spca/issues) describing
the result you expected and what happened. Include a small runnable example,
the relevant parameter settings, `sessionInfo()`, and the package version.
For randomized fits, include `set.seed()` and a description of approximation
quality. Avoid attaching private data; a simulated example is usually easier
to investigate and share.

## Development setup

```r
install.packages(c("Matrix", "testthat", "roxygen2"))
```

The repository is named `spca`; the installed R package is `sparsepca`.
From a checkout, run the package tests with:

```r
testthat::test_local()
```

Build and check the source archive from the repository's parent directory:

```sh
R CMD build spca
R CMD check --as-cran sparsepca_1.0.0.tar.gz
```

Suggested packages and documentation tools must be available for a full
check. A check of an installed package or working directory is useful
during development; release evidence should identify the exact source
archive, its checksum, source commit, R version, dependency versions, and
any checks the environment could not run.

## Numerical changes

Start with a small regression or mathematical reference test. Useful
assertions compare objective values, residual norms, subspaces, sparsity,
and public result shapes. Singular-vector signs, ordering within ties, and
bases in repeated singular subspaces can vary across numerical libraries.
Use justified tolerances and compare the represented model where possible.
Avoid assertions on a signed sum of residual entries: cancellation can hide
a large error.

Preserve public arguments and defaults unless a documented correction is
necessary. Explain how a change affects loadings, preprocessing, variance
statistics, prediction, or downstream use. Sparse loadings and sparse input
storage are separate concerns. When adding sparse operations, compare their
products with an explicit dense reference before testing an entire fit.

Large benchmarks belong under `tools/benchmarks`, outside the routine test
suite and the CRAN source archive. Save raw measurements and environment
details; generate figures from those saved measurements. Use comparable
objectives, preprocessing, and requested outputs, and report approximation
quality beside time and memory. Do not use fragile timing assertions in
unit tests.

## Documentation and release notes

The authoritative help text is maintained with roxygen comments in
`R/documentation.R`; exported functions and S3 registrations are declared
beside their implementations. Regenerate help pages and `NAMESPACE` with:

```r
roxygen2::roxygenise()
```

Update `NEWS.md` for user-visible changes. Run examples and inspect generated
HTML and PDF documentation. The README should teach the public API with
small examples that can run without external data downloads; any larger
visual demonstration must record its data provenance and reuse conditions.

For a release candidate, build once and check the identical archive on the
supported environments. Keep check findings, including NOTES and environment
limitations, in the release record. An unexecuted workflow is configuration,
not evidence of a successful check.
