# Documentation examples and illustrations

The loading illustration uses simulated data generated entirely by
`illustrate-loadings.R`. Twelve variables form three groups, with four
noisy measurements per latent signal. The latent signals are correlated,
so ordinary PCA mixes the groups. The illustration compares three-component
ordinary PCA with elastic net and hard-threshold sparse fits on the same
centered and scaled data.

Run from the repository root, with `pkgload` installed:

```sh
Rscript tools/documentation/illustrate-loadings.R
```

The script saves raw loading coefficients, model parameters, nonzero counts,
relative Frobenius reconstruction errors, convergence diagnostics, and R
session information under `results/`. It generates SVG and PNG figures in
`man/figures/`. Signs are flipped for display only so each component's
largest-magnitude coefficient is positive; this does not change sparsity
or reconstruction. The original signed coefficients are retained in the
CSV file. Colors use one common coefficient scale, and white denotes zero.

The fits use different penalties to illustrate interpretation and
reconstruction tradeoffs. This is not a speed benchmark or a claim that
the penalty strengths are mathematically equivalent. All illustration data
are synthetic; no external images or human-subject dataset is needed.

The scripts and README-only figure assets are excluded from the CRAN source
archive. The help examples are generated from `R/documentation.R` and run
by `R CMD check`.

## Verify the examples and HTML help

The following development check loads the source package with `pkgload`,
executes the examples in every help topic, and executes the README's R
blocks in document order. It skips the installation block explicitly.
It saves example output, any warnings or errors, plots, a results table,
HTML help, and session information in the chosen output directory:

```sh
Rscript tools/documentation/check-examples.R . /tmp/sparsepca-documentation
```

This is a source-checkout documentation check. Release validation also
builds a source archive and runs `R CMD check --as-cran` on that archive,
including installed-package examples and the complete PDF manual.
