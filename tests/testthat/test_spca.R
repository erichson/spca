#devtools::use_package("testthat")

context("Sparse PCA")

#Load rsvd library
library(spca)

#Set seed
set.seed(1234)

#Accuray
atol_float64 <- 1e-8

#--------------------------------------------------------------------
# Generate Some Data in R
#--------------------------------------------------------------------
m = 1000
V1 = rnorm(m, 0, 290)
V2 = rnorm(m, 0, 300)
V3 = -0.1*V1 + 0.1*V2 + rnorm(m,0,100)

X = cbind(V1,V1,V1,V1, V2,V2,V2,V2, V3,V3)
X = X + matrix(rnorm(length(X),0,1), ncol = ncol(X), nrow = nrow(X))


#*************************************************************************************
# Test: SPARSE PCA - center = TRUE, scale. = TRUE
#*************************************************************************************

pca_out <- prcomp(X, center = TRUE, scale. = FALSE) #Deterministic PCA
spca_out <- spca(X, k=3, alpha=0, beta=0, center = TRUE, scale = FALSE, verbose=0) #Sparse PCA

#Test1: SPCA recovers PCA for alpha = 0 and beta = 0.
testthat::test_that("Test 1: cov; alpha = 0 and beta = 0", {
  testthat::expect_equal(pca_out$sdev[1:3], spca_out$sdev[1:3])
  testthat::expect_equal(sum(diag(1,3,3) - t(spca_out$loadings)%*%spca_out$loadings), 0 )
})


pca_out <- prcomp(X, center = TRUE, scale. = TRUE) #Deterministic PCA
spca_out <- spca(X, k=3, alpha=0, beta=0, center = TRUE, scale = TRUE, verbose=0) #Sparse PCA

#Test1: SPCA recovers PCA for alpha = 0 and beta = 0.
  testthat::test_that("Test 2: cor; alpha = 0 and beta = 0", {
    testthat::expect_equal(pca_out$sdev[1:3], spca_out$sdev[1:3])
    testthat::expect_equal(sum(diag(1,3,3) - t(spca_out$loadings)%*%spca_out$loadings), 0 )
  })





#*************************************************************************************
# Test: Sparse PCA - Reconstruction error
#*************************************************************************************
Re <- pca_out$x %*% t(pca_out$rotation)
Re2 <- spca_out$scores %*% t(spca_out$transform)

testthat::test_that("Test 3: Sparse PCA reconstruction error", {
  testthat::expect_lt(sum(Re-Re2), atol_float64)
})


#*************************************************************************************
# Test: Sparse PCA - Sparse Loadings
#*************************************************************************************
spca_out <- spca(X, k=3, alpha=1e-3, beta=0, center = TRUE, scale = FALSE, verbose=0) #Sparse PCA

testthat::test_that("Test 4: Sparse PCA; sparse loadings with alpha = 1e-3", {
  testthat::expect_equal(colSums(abs(spca_out$loadings)>0.0), c(4,4,2))
})


spca_out <- spca(X, k=3, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0) #Sparse PCA

testthat::test_that("Test 5: Sparse PCA; sparse loadings with alpha = 1e-3 and beta = 1e-3", {
  testthat::expect_equal(colSums(abs(spca_out$loadings)>0.0), c(4,4,2))
})


#*************************************************************************************
# Test: Sparse PCA - Reconstruction error with sparse loadings
#*************************************************************************************
Re3 <- spca_out$scores %*% t(spca_out$transform)

testthat::test_that("Test 6: Sparse PCA reconstruction error II", {
  testthat::expect_lt(sum(Re-Re3), atol_float64)
})

