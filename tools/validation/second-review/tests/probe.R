library(sparsepca)
C <- diag(c(8, 4, 2, 1))
A <- diag(4)
for (fraction in c(0.25, 1 - 1e-9)) {
  B <- fraction * A
  CA <- C %*% A
  CB <- C %*% B
  exact <- 0.5 * sum((C - CB %*% t(A))^2)
  expanded <- 0.5 * (sum(C^2) - 2 * sum(CA * CB) + sum(CB^2))
  got <- sparsepca:::.spca_loss(C, CA, CB, A, sum(C^2))
  print(c(fraction=fraction, direct=exact, expanded=expanded, got=got,
          observed_ratio=got/exact, expanded_ratio=expanded/exact))
}
X <- rbind(diag(5), -diag(5))
B <- cbind(c(0,0,0,0,0), c(1,0,0,0,0), c(3,0,0,0,0),
           c(1,2,0,0,0), c(1,1,3,0,0))
p <- sparsepca:::.spca_preprocess(X, center=FALSE, scale=FALSE)
f <- list(B=B, A=diag(5), objective=0, converged=TRUE, improvement=0)
r <- sparsepca:::.spca_result(f,p,'spca',0,0,'elasticnet')
print(r$eigenvalues)
print(r$score_variance)
print(list(qr_pivot=qr(r$scores)$pivot, rank=qr(r$scores)$rank))
progress <- sparsepca:::.spca_progress
print(progress(2,1,tol=.1,evaluation_scale=1000))
print(progress(1,1+1e-12,tol=0,evaluation_scale=1000))
print(progress(1,1.001,tol=1e-5,evaluation_scale=1000))
print(progress(0,0,tol=0,evaluation_scale=0))
print(sessionInfo())
