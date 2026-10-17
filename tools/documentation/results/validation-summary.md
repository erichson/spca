# Documentation validation for sparsepca 1.0.0

## Scope and source identity

This record concerns the source checkout used to prepare the release candidate.
The source files and rendered manual are identified by SHA-256 in
`documentation-source-manifest.json`. The recorded Git commit may precede
uncommitted release documentation; file hashes identify the inspected content.
This is not a substitute for checking the exact source archive. Release archive
checks, operating-system coverage, and installed-package checks are recorded
separately by the release workflow.

## Environment

The examples ran under R 4.6.1 (2026-06-24), x86_64-pc-linux-gnu, Ubuntu
24.04.3 LTS, with Matrix 1.7-6, rsvd 1.0.6, pkgload 1.5.3, and testthat 3.3.2.
Documentation was generated with roxygen2 8.1.1. OpenBLAS 0.3.26 supplied BLAS
and LAPACK 3.12.0; OPENBLAS_NUM_THREADS and OMP_NUM_THREADS were both one.
`example-session-info.txt` records the full R session, locale, and time zone.

## Executed examples

The command was:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
R_LIBS_USER=/workspace/scratch/4668db2fbf57/validation/R-library \
/opt/R/4.6.1/bin/Rscript tools/documentation/check-examples.R . \
  /workspace/scratch/4668db2fbf57/validation/documentation/examples
```

All six help topics containing executable examples passed: `spca`, `rspca`,
`robspca`, `predict.spca`, `spca_path`, and `spca-methods`. The package overview
has no executable examples. All eight executable README R blocks passed in
document order. The first README block, which gives package-installation
commands, was explicitly skipped; no network installation was performed by
this documentation check. There were no example warnings or errors.

`example-results.csv` records each block's status. The adjacent `help-*.log`
and `README-*.log` files retain numerical and printed output. Benchmark prose,
tables, and images were added after the final example run; no executable R
blocks changed in that step.

## Generated help and PDF manual

All seven HTML help topics were generated with `tools::Rd2HTML` without
warnings. Local browser inspection is recorded separately after rendering.

The full indexed manual was built with:

```sh
TEXMFHOME=/workspace/scratch/4668db2fbf57/validation/texmf \
TEXMFCONFIG=/workspace/scratch/4668db2fbf57/validation/texmf-config \
TEXMFVAR=/workspace/scratch/4668db2fbf57/validation/texmf-var \
/opt/R/4.6.1/bin/R CMD Rd2pdf --no-preview --force \
  --output=/workspace/scratch/4668db2fbf57/validation/documentation/sparsepca-manual.pdf \
  /workspace/scratch/4668db2fbf57/spca
```

The complete manual contains 19 pages, including its index. All pages were
rendered with Poppler and visually inspected; equations, argument tables,
examples, references, and index entries are readable without clipping or
overlap. After the last wording clarification, the affected path pages were
rendered and inspected again. The final TeX log contains no overfull boxes.
The manual PDF and page renders remain in the separate validation directory;
they are not shipped inside the CRAN source archive.

The successful build used texi2dvi emulation because the standalone texi2dvi
program was unavailable. TeX reported an unavailable italic Inconsolata font
shape in the index and substituted its normal shape. These environment
messages are preserved in `manual-build.log`; neither prevented creation of
the complete indexed manual. Intermediate rerunfilecheck messages resolved
during normal repeated LaTeX passes. They must not be described as zero
messages from TeX or confused with the separately recorded R CMD check status.

## Figures and their provenance

The loading illustration was generated from the synthetic data in
`tools/documentation/illustrate-loadings.R`. Both SVG and PNG versions were
visually checked. The raw coefficients, parameter values, reconstruction
errors, and convergence diagnostics are retained in `loadings-*.csv`.
Signs are changed only for display. The common legend, group labels, and
sparsity/error annotations are readable and do not overlap.

Benchmark figures and raw performance measurements have their own scripts,
source hashes, and methodology under `tools/benchmarks`. Historical eigenface
illustrations are retained with an explicit provenance and reuse caveat;
they are not presented as newly computed 1.0.0 results.
