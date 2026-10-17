.libPaths(c('/workspace/scratch/4668db2fbf57/validation/R-library', .libPaths()))
repo <- '/workspace/scratch/4668db2fbf57/spca'
for (file in list.files(file.path(repo, 'R'), full.names = TRUE)) source(file)
x <- c(-5, -2, -1, .5, 1, 3, 20)
X <- matrix(x, ncol = 1)
huber <- function(r, gamma) ifelse(abs(r) <= gamma, 0.5 * r^2, gamma * abs(r) - 0.5 * gamma^2)
scalar_reference <- function(alpha, beta, gamma) {
  spectral2 <- sum(x^2)
  score <- function(b) -sum(abs(x) * pmin(abs(x) * (1 - b), gamma)) + alpha * spectral2 + beta * spectral2 * b
  if (score(0) >= 0) return(0)
  uniroot(score, c(0, 1), tol = 1e-14)$root
}
scalar_rows <- list()
for (settings in list(c(.01, .05, .5), c(.01, .001, .1), c(.1, .1, 0), c(0, 0, 0), c(.001, 0, .5))) {
 alpha <- settings[1]; beta <- settings[2]; gamma <- settings[3]
 fit <- robspca(X, 1, alpha = alpha, beta = beta, gamma = gamma, center = FALSE, tol = 1e-13, max_iter = 20000, verbose = FALSE)
 reference <- if (alpha + beta == 0 && gamma == 0) NA_real_ else scalar_reference(alpha, beta, gamma)
 scalar_rows[[length(scalar_rows) + 1]] <- data.frame(alpha, beta, gamma, coefficient = abs(fit$loadings[1]), reference, iterations=fit$iterations, converged=fit$converged, objective=tail(fit$objective,1), reconstruction_difference=max(abs(X-fit$scores%*%t(fit$transform)-fit$sparse)))
}
scalar <- do.call(rbind, scalar_rows)
write.csv(scalar, '/workspace/scratch/4668db2fbf57/validation/second-review/solver/scalar-robust-results.csv', row.names=FALSE)
print(scalar)
rows <- list()
for(seed in 1:10) {
 set.seed(seed)
 X <- matrix(rnorm(60*8),60,8)
 X[,5:8] <- X[,5:8] + .8*X[,1:4]
 X[1,1] <- X[1,1] + 15
 xc <- sweep(X,2,colMeans(X),'-')
 for(method in c('spca','robspca')) for(tol in c(1e-5,1e-10)) {
  args <- list(X=X,k=2,alpha=.02,beta=.01,tol=tol,max_iter=5000,verbose=FALSE)
  if (method == 'robspca') args$gamma <- .5
  fit <- do.call(get(method), args)
  XS <- if(method == 'robspca') xc-fit$sparse else xc
  M <- crossprod(XS, xc%*%fit$loadings)
  decomp <- svd(M,nu=2,nv=2)
  bestA <- tcrossprod(decomp$u,decomp$v)
  current_loss <- sum((XS-fit$scores%*%t(fit$transform))^2)/2
  polished_loss <- sum((XS-fit$scores%*%t(bestA))^2)/2
  bestAres <- norm(M-fit$transform%*%((crossprod(fit$transform,M)+crossprod(M,fit$transform))/2),'F') / max(norm(M,'F'),1)
  rows[[length(rows)+1]] <- data.frame(seed,method,tol,iterations=fit$iterations,converged=fit$converged,objective=tail(fit$objective,1),polar_gap=current_loss-polished_loss,polar_relative=(current_loss-polished_loss)/max(tail(fit$objective,1),1),polar_stationarity=bestAres)
 }
}
polars <- do.call(rbind, rows)
write.csv(polars, '/workspace/scratch/4668db2fbf57/validation/second-review/solver/final-polar-results.csv',row.names=FALSE)
print(aggregate(cbind(polar_relative,polar_stationarity,iterations)~method+tol,polars,max))
