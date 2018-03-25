Sparse Principal Component Analysis (SPCA) via Variable Projection
*******************************************************************

Sparse principal component analysis is a modern variant of PCA. Specifically, SPCA attempts to find sparse 
weight vectors (loadings), i.e., a weight vector with only a few `active' (nonzero) values. This approach 
leads to an improved interpretability of the model, because the principal components are formed as a 
linear combination of only a few of the original variables. Further, SCPA avoids overfitting in a 
high-dimensional data setting where the number of variables is greater than the number of observations.

This package provides robust and randomized accelerated SPCA routines in R:
 
* Sparse PCA: ``spca(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, scale=FALSE,...)``.
* Randomized accelerated sparse PCA: ``rspca(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, scale=FALSE,...)``.
* Robust sparse PCA: ``robspca(X, k=NULL, alpha=1e-4, beta=1e-4, gamma=100, center=TRUE, scale=FALSE,...)``.

All the implementations are using variable projection as optimization strategy. 


Installation
************

Install the spca package via CRAN
```R
install.packages("spca")
```

You can also install the development version from GitHub using [devtools](https://cran.r-project.org/package=devtools):

```r
devtools::install_github("erichson/spca")
```

The source packge can be obtained here: [CRAN: rsvd](https://cran.r-project.org/web/packages/rsvd/index.html).

Example
*******


References
*************
* [N. Benjamin Erichson, et al. Randomized Matrix Decompositions using R. (2016)](http://arxiv.org/abs/1608.02148)
* [Sergey Voronin, Per-Gunnar Martinsson. RSVDPACK: Subroutines for computing partial singular value decompositions via randomized sampling on single core, multi core, and GPU architectures. (2015)](https://arxiv.org/abs/1502.05366)
* [Nathan Halko, et al. Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. (2011)](https://arxiv.org/abs/0909.4061)
