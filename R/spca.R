#' @title  Sparse Principal Component Analysis (spca).
#
#' @description Implementation of SPCA, using variable projection as an optimization strategy.
#
#' @details
#' Sparse principal component analysis is a modern variant of PCA. Specifically, SPCA attempts to find sparse
#' weight vectors (loadings), i.e., a weight vector with only a few 'active' (nonzero) values. This approach
#' leads to an improved interpretability of the model, because the principal components are formed as a
#' linear combination of only a few of the original variables. Further, SPCA avoids overfitting in a
#' high-dimensional data setting where the number of variables \eqn{p} is greater than the number of
#' observations \eqn{n}.
#'
#' Such a parsimonious model is obtained by introducing prior information like sparsity promoting regularizers.
#' More concreatly, given an \eqn{(n,p)} data matrix \eqn{X}, SPCA attemps to minimize the following
#' objective function:
#'
#' \deqn{ f(A,B) = \frac{1}{2} \| X - X B A^\top \|^2_F + \psi(B) }
#'
#' where \eqn{B} is the sparse weight (loadings) matrix and \eqn{A} is an orthonormal matrix.
#' \eqn{\psi} denotes a sparsity inducing regularizer such as the LASSO (\eqn{\ell_1}{l1} norm) or the elastic net
#' (a combination of the \eqn{\ell_1}{l1} and \eqn{\ell_2}{l2} norm). The principal components \eqn{Z} are formed as
#'
#' \deqn{ Z = X B }{Z = X * B}
#'
#' and the data can be approximately rotated back as
#'
#' \deqn{ \tilde{X} = Z A^\top }{ X = Z t(A)}
#'
#' The print and summary method can be used to present the results in a nice format.
#'
#' @param X       array_like; \cr
#'                a real \eqn{(n, p)} input matrix (or data frame) to be decomposed.
#'
#' @param k       integer; \cr
#'                specifies the target rank, i.e., the number of components to be computed.
#'
#' @param alpha   float; \cr
#'                Sparsity controlling parameter. Higher values lead to sparser components.
#'
#' @param beta    float; \cr
#'                Amount of ridge shrinkage to apply in order to improve conditioning.
#'
#' @param center  bool; \cr
#'                logical value which indicates whether the variables should be
#'                shifted to be zero centered (TRUE by default).
#'
#' @param scale   bool; \cr
#'                logical value which indicates whether the variables should
#'                be scaled to have unit variance (FALSE by default).
#'
#' @param max_iter integer; \cr
#'                 maximum number of iterations to perform before exiting.
#'
#' @param tol float; \cr
#'            stopping tolerance for the convergence criterion.
#'
#' @param verbose bool; \cr
#'                logical value which indicates whether progress is printed.
#'
#'
#'
#'@return \code{spca} returns a list containing the following three components:
#'\item{loadings}{  array_like; \cr
#'           sparse loadings (weight) vector;  \eqn{(p, k)} dimensional array.
#'}
#'
#'\item{transform}{  array_like; \cr
#'           the approximated inverse transform; \eqn{(p, k)} dimensional array.
#'}
#'
#'\item{scores}{  array_like; \cr
#'           the principal component scores; \eqn{(n, k)} dimensional array.
#'}
#'
#'\item{eigenvalues}{  array_like; \cr
#'          the approximated eigenvalues; \eqn{(k)} dimensional array.
#'}
#'
#'\item{center, scale}{  array_like; \cr
#'                     the centering and scaling used.
#'}
#'
#'
#' @references
#' \itemize{
#'
#'  \item [1] N. B. Erichson, P. Zheng, K. Manohar, S. Brunton, J. N. Kutz, A. Y. Aravkin.
#'  "Sparse Principal Component Analysis via Variable Projection."
#'  Submitted to IEEE Journal of Selected Topics on Signal Processing (2018).
#'  (available at `arXiv \url{https://arxiv.org/abs/1804.00341}).
#' }
#'
#' @author N. Benjamin Erichson, Peng Zheng, and Sasha Aravkin
#'
#' @seealso \code{\link{rspca}}, \code{\link{robspca}}
#'
#' @examples
#'
#' # Create artifical data
#' m <- 10000
#' V1 <- rnorm(m, 0, 290)
#' V2 <- rnorm(m, 0, 300)
#' V3 <- -0.1*V1 + 0.1*V2 + rnorm(m,0,100)
#'
#' X <- cbind(V1,V1,V1,V1, V2,V2,V2,V2, V3,V3)
#' X <- X + matrix(rnorm(length(X),0,1), ncol = ncol(X), nrow = nrow(X))
#'
#' # Compute SPCA
#' out <- spca(X, k=3, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
#' print(out)
#' summary(out)
#'
#'


#' @export
spca <- function(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, scale=FALSE, max_iter=1000, tol=1e-5, verbose=TRUE) UseMethod("spca")

#' @export
spca.default <- function(X, k=NULL, alpha=1e-4, beta=1e-4, center=TRUE, scale=FALSE, max_iter=1000, tol=1e-5, verbose=TRUE) {

  X <- as.matrix(X)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checks
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (any(is.na(X))) {
    warning("Missing values are omitted: na.omit(X).")
    X <- stats::na.omit(X)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Init rpca object
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  spcaObj = list(loadings = NULL,
                 transform = NULL,
                 scores = NULL,
                 eigenvalues = NULL,
                 center = center,
                 scale = scale)

  n <- nrow(X)
  p <- ncol(X)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set target rank
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(k)) k <- min(n,p)
  if(k > min(n,p)) k <- min(n,p)
  if(k<1) stop("Target rank is not valid!")



  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Center/Scale data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(center == TRUE) {
    spcaObj$center <- colMeans(X)
    X <- sweep(X, MARGIN = 2, STATS = spcaObj$center, FUN = "-", check.margin = TRUE)
  } else { spcaObj$center <- FALSE }

  if(scale == TRUE) {
    spcaObj$scale <- sqrt(colSums(X**2) / (n-1))
    if(is.complex(spcaObj$scale)) { spcaObj$scale[Re(spcaObj$scale) < 1e-8 ] <- 1+0i
    } else {spcaObj$scale[spcaObj$scale < 1e-8] <- 1}
    X <- sweep(X, MARGIN = 2, STATS = spcaObj$scale, FUN = "/", check.margin = TRUE)
  } else { spcaObj$scale <- FALSE }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute SVD for initialization of the Variable Projection Solver
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  svd_init <- svd(X)

  Dmax <- svd_init$d[1] # l2 norm

  A <- svd_init$v[,1:k]
  B <- svd_init$v[,1:k]

  V <- svd_init$v
  VD = sweep(V, MARGIN = 2, STATS = svd_init$d, FUN = "*", check.margin = TRUE)
  VD2 = sweep(V, MARGIN = 2, STATS = svd_init$d**2, FUN = "*", check.margin = TRUE)


  #--------------------------------------------------------------------
  #   Set Tuning Parameters
  #--------------------------------------------------------------------
  alpha <- alpha *  Dmax**2
  beta <- beta * Dmax**2

  nu <- 1.0 / (Dmax**2 + beta)
  kappa <- nu * alpha

  obj <- c()
  improvement <- Inf

  #--------------------------------------------------------------------
  #   Apply Variable Projection Solver
  #--------------------------------------------------------------------
  noi <- 1
  while (noi <= max_iter && improvement > tol) {

        # Update A:  X'XB = UDV'
        Z <- VD2 %*% (t(V) %*% B)
        svd_update <- svd(Z)
        A <- svd_update$u %*% t(svd_update$v)


        # Proximal Gradient Descent to Update B:
        grad <- VD2 %*% (t(V) %*% (A - B)) - beta * B
        B_temp <- B + nu * grad

        # l1 soft-threshold
        idxH <- which(B_temp > kappa)
        idxL <- which(B_temp <= -kappa)

        B = matrix(0, nrow = nrow(B_temp), ncol = ncol(B_temp))
        B[idxH] <- B_temp[idxH] - kappa
        B[idxL] <- B_temp[idxL] + kappa



      # compute residual
      R <- t(VD) - (t(VD) %*% B) %*% t(A)

      # compute objective function
      obj <- c(obj, 0.5 * sum(R**2) + alpha * sum(abs(B)) + 0.5 * beta * sum(B**2))

      # Break if obj is not improving anymore
      if(noi > 1){
        improvement <- (obj[noi-1] - obj[noi]) / obj[noi]
      }

      # Trace
      if(verbose > 0 && (noi-1) %% 10 == 0) {
        print(sprintf("Iteration: %4d, Objective: %1.5e, Relative improvement %1.5e", noi, obj[noi], improvement))
      }


      # Next iter
      noi <- noi + 1


  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Update spca object
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  spcaObj$loadings <- B
  spcaObj$transform <- A
  spcaObj$scores <- X %*% B
  spcaObj$eigenvalues <- svd_update$d / (n - 1)
  spcaObj$objective <- obj


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Explained variance and explained variance ratio
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  spcaObj$sdev <-  sqrt( spcaObj$eigenvalues )
  spcaObj$var <- sum( apply( Re(X) , 2, stats::var ) )
  if(is.complex(X)) spcaObj$var <- Re(spcaObj$var + sum( apply( Im(X) , 2, stats::var ) ))


  class(spcaObj) <- "spca"
  return( spcaObj )

}


#' @export
print.spca <- function(x , ...) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Standard deviations:\n")
  print(round(x$sdev, 3))
  cat("\nEigenvalues:\n")
  print(round(x$eigenvalues, 3))
  cat("\nSparse loadings:\n")
  print(round(x$loadings, 3))
}


#' @export
summary.spca <- function( object , ... )
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Summary spca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  variance = object$sdev**2
  explained_variance_ratio = variance / object$var
  cum_explained_variance_ratio = cumsum( explained_variance_ratio )

  x <- t(data.frame( var = round(variance, 3),
                     sdev = round(object$sdev, 3),
                     prob = round(explained_variance_ratio, 3),
                     cum = round(cum_explained_variance_ratio, 3)))

  rownames( x ) <- c( 'Explained variance',
                      'Standard deviations',
                      'Proportion of variance',
                      'Cumulative proportion')

  colnames( x ) <- paste(rep('PC', length(object$sdev)), 1:length(object$sdev), sep = "")

  x <- as.matrix(x)

  return( x )
}
