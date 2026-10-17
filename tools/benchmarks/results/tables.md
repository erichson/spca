# Measured benchmark results

Three isolated processes per method and workload; medians are reported. Brackets show the observed elapsed-time range. Speedup is the reference median divided by the candidate median. Full methodology is in `tools/benchmarks/README.md`.

## Sparse loadings on approximately low-rank dense data

| Matrix | SPCA seconds [range] | RSPCA seconds [range] | Speedup | Relative error, SPCA / RSPCA | Zero loadings |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2,000 × 600 | 1.462 [0.812, 1.631] | 0.192 [0.176, 0.531] | 7.61× | 0.07282 / 0.07282 | 90.0% |
| 5,000 × 1,500 | 7.446 [6.682, 7.672] | 1.624 [1.196, 2.313] | 4.58× | 0.07735 / 0.07739 | 90.0% |
| 10,000 × 3,000 | 45.386 [41.321, 52.049] | 4.850 [3.108, 6.392] | 9.36× | 0.08365 / 0.08365 | 90.0% |

All fits use rank 10, `alpha = 0.001`, `beta = 0.0001`, centering, no scaling, `tol = 1e-5`, and at most 500 iterations. RSPCA uses `o = 20`, `q = 2`. The two methods use the same relative penalties; the complete CSV also evaluates both fits using the exact deterministic spectral scale.

## External PCA reference and a flat-spectrum limitation

| Workload | Method | Seconds | Speedup vs prcomp | Relative error | Fitted energy / input energy |
| --- | --- | ---: | ---: | ---: | ---: |
| Low-rank signal, 10,000 × 3,000 | stats::prcomp | 36.939 | 1.00× | 0.07107 | 99.495% |
| Low-rank signal, 10,000 × 3,000 | Randomized SPCA | 7.283 | 5.07× | 0.07107 | 99.495% |
| Gaussian noise, 2,000 × 600 | stats::prcomp | 0.967 | 1.00× | 0.98073 | 3.816% |
| Gaussian noise, 2,000 × 600 | Randomized SPCA | 0.256 | 3.78× | 0.98261 | 3.448% |

The external comparison sets `alpha = beta = 0`, so both methods target ordinary rank-10 PCA with the same preprocessing and returned scores. `prcomp(rank. = 10)` uses R's deterministic SVD; the randomized fit is approximate. This does not compare sparsity-penalized fitting with ordinary PCA.

## Peak resident memory

| Workload | Method | Peak RSS at fit completion (MiB) | Seconds | Relative error |
| --- | --- | ---: | ---: | ---: |
| Dense input, 10,000 × 3,000 | SPCA | 1716.2 | 45.386 | 0.08365 |
| Dense input, 10,000 × 3,000 | Randomized SPCA | 1122.2 | 4.850 | 0.08365 |
| 10% stored input, 10,000 × 3,000 | Dense conversion + RSPCA | 1354.3 | 7.065 | 0.05559 |
| 10% stored input, 10,000 × 3,000 | Sparse input + RSPCA | 342.2 | 2.733 | 0.05559 |

RSS is the whole R process high-water mark, including runtime, namespaces, input loading, a small warm-up, preprocessing, fitting, and result construction. It is read immediately after fitting and before blocked quality validation. It is neither incremental memory nor cumulative allocation. GNU time's additional whole-process peak (including validation) is also retained in the raw data.
