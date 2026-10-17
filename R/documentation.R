#' Sparse principal component analysis
#'
#' Fit interpretable component scores as sparse linear combinations of the
#' input variables. Variable projection alternates an exact orthogonal
#' reconstruction update with a proximal step for the loading coefficients.
#'
#' @name spca
#' @aliases spca.default
#' @param X A real numeric matrix, data frame with numeric columns, or
#'   real \pkg{Matrix} object. Rows are observations and columns are
#'   variables. Rows containing missing values are omitted with a warning;
#'   at least two complete observations and one variable are required.
#'   Logical inputs are treated as zero/one. Infinite and complex values
#'   are not supported. \code{spca} converts
#'   \pkg{Matrix} input to a dense matrix; see \code{\link{rspca}} for
#'   computation that preserves sparse input storage.
#' @param k Positive integer number of components. \code{NULL} uses
#'   \code{min(nrow(X), ncol(X))} after incomplete rows are omitted. Larger
#'   values are reduced to that limit. A component can have zero scores
#'   when the data are rank deficient or its loadings are all zero.
#' @param alpha Finite nonnegative sparsity penalty, relative to the squared
#'   largest singular value of the processed data. The effective penalty is
#'   \code{alpha * d1^2}; see Details. Increasing \code{alpha} encourages
#'   fewer nonzero coefficients but does not impose an exact feature count.
#' @param beta Finite nonnegative ridge penalty, using the same relative
#'   scaling as \code{alpha}. The default is \code{1e-4}. Ridge shrinkage
#'   can improve conditioning for strongly correlated variables.
#' @param center Logical; subtract column means before fitting (default
#'   \code{TRUE}). The fitted means are stored for prediction.
#' @param scale Logical; divide each processed column by its scale factor
#'   (default \code{FALSE}). Factors are square roots of column sums of
#'   squares divided by \eqn{n-1}, after any centering. This is the sample
#'   standard deviation for centered data and a root-mean-square quantity
#'   otherwise. Factors below \code{1e-8} are replaced by one.
#' @param max_iter Positive integer limit on optimization iterations.
#'   Reaching the limit does not establish convergence.
#' @param tol Finite nonnegative tolerance for relative objective change.
#'   This controls an optimization stopping rule, not a reconstruction-error
#'   bound or a statistical accuracy guarantee. An objective increase larger
#'   than a roundoff allowance stops with a warning and does not count as
#'   convergence; see Details.
#' @param verbose Logical, or the equivalent numeric zero/one; print
#'   optimization progress when enabled.
#' @param penalty Sparsity penalty: \code{"elasticnet"} (the default)
#'   uses an entrywise absolute-value penalty and soft thresholding;
#'   \code{"l0"} counts nonzero coefficients and uses hard thresholding.
#'   With \code{penalty = "l0"} and \code{beta > 0}, the combination of
#'   cardinality and ridge penalties is also called ZeroNet.
#' @param B_init Optional real finite \eqn{p} by \eqn{k} matrix of initial
#'   loading coefficients. \code{NULL} initializes from leading right
#'   singular vectors. A previous fit's \code{loadings} can supply a warm
#'   start when variable order, component count, and preprocessing match.
#'   Initial loadings are not normalized.
#'
#' @details
#' Let \eqn{X} denote the retained, centered and/or scaled data, with
#' \eqn{n} observations and \eqn{p} variables. The default method minimizes
#' \deqn{\frac{1}{2}\|X-XBA^T\|_F^2 +
#'       \lambda_1\|B\|_1 + \frac{\lambda_2}{2}\|B\|_F^2}
#' subject to \eqn{A^T A=I}. Here \eqn{A} and \eqn{B} are \eqn{p} by
#' \eqn{k}; \eqn{\|B\|_1} is the sum of absolute entries, and
#' \eqn{\|\cdot\|_F} is the Frobenius norm. The effective coefficients
#' are \eqn{\lambda_1=\mathtt{alpha}\,d_1^2} and
#' \eqn{\lambda_2=\mathtt{beta}\,d_1^2}, where \eqn{d_1} is the largest
#' singular value of the processed matrix. For \code{penalty = "l0"},
#' replace \eqn{\|B\|_1} with the number of nonzero entries of \eqn{B}.
#' The same numeric \code{alpha} does not imply equivalent sparsity under
#' the two penalties.
#'
#' For fixed \eqn{B}, an SVD of \eqn{X^T X B} gives the best orthogonal
#' reconstruction factor \eqn{A}. A gradient step followed by thresholding
#' updates \eqn{B}. The joint problem is nonconvex, so the solution can
#' depend on initialization. Small objective change is a convergence
#' diagnostic, not proof of a global minimum. With \code{alpha = beta = 0},
#' deterministic fitting recovers ordinary PCA with matching preprocessing,
#' up to signs and bases within repeated singular subspaces.
#'
#' Relative improvement divides the decrease in objective by the magnitude
#' of the current objective, with a positive machine denominator at zero.
#' Objective increases within \code{64 * .Machine$double.eps} times the
#' largest magnitude of the consecutive objectives or squared data norm
#' are treated as numerical stagnation. Larger increases stop with a
#' warning and \code{converged = FALSE}. For \code{rspca}, the squared
#' data norm in this allowance is that of the sketch.
#'
#' Scores are \eqn{Z=XB}; reconstruction in processed coordinates is
#' \eqn{ZA^T}. The loading matrix \eqn{B} retains its penalized scale.
#' It is generally neither normalized nor orthogonal. The reconstruction
#' factor \eqn{A} is not a matrix inverse of \eqn{B}. To reconstruct in
#' input coordinates, multiply columns of \eqn{ZA^T} by the stored scales
#' and then add the stored means, omitting operations stored as
#' \code{FALSE}.
#'
#' @section Variance and loading interpretation:
#' Sparse component scores can be correlated. The historical
#' \code{eigenvalues} field contains sequentially adjusted score energies:
#' project each score column onto the span of preceding columns, subtract
#' that projection, and divide the squared residual norm by \eqn{n-1}.
#' This is equivalent to squared diagonal entries in an unpivoted QR
#' decomposition of the scores, with dependent columns assigned zero.
#' The adjustment follows Zou, Hastie and Tibshirani (2006), applied to
#' the returned, unnormalized loading coefficients. It depends on component
#' order and loading scale; values are not necessarily decreasing.
#'
#' For centered data these are adjusted sample variances. They are not
#' generally eigenvalues of a covariance matrix. \code{score_variance}
#' contains the marginal squared score norms divided by \eqn{n-1},
#' without adjustment for overlapping components. If centering is disabled,
#' both quantities are uncentered second moments with that denominator,
#' rather than centered sample variances.
#'
#' \code{var} is the total squared energy of the entire processed input
#' divided by \eqn{n-1}; for centered data it is the exact sum of input
#' sample variances. Summary proportions use this denominator, including
#' variation in omitted directions. They need not sum to one. For zero
#' total energy, variance proportions are mathematically undefined; the
#' summary displays zero by convention.
#'
#' Loadings are coefficients, not variable--component correlations.
#' For complete training data, correlations can be computed directly as
#' \code{cor(as.matrix(X), fit$scores)}. Multiplying loadings by
#' \code{sdev} does not generally produce these correlations for sparse
#' components.
#'
#' @return A list of class \code{"spca"} containing:
#' \describe{
#'   \item{loadings}{A dense \eqn{p} by \eqn{k} matrix \eqn{B} of fitted
#'     sparse coefficients. Zeros are stored in an ordinary R matrix.}
#'   \item{transform}{A dense \eqn{p} by \eqn{k} reconstruction factor
#'     \eqn{A} with orthonormal columns.}
#'   \item{scores}{A dense \eqn{n} by \eqn{k} matrix \eqn{XB}. Rows
#'     correspond to complete observations retained from the input.}
#'   \item{eigenvalues}{Numeric vector of length \eqn{k} containing
#'     sequentially adjusted score energies divided by \eqn{n-1}; see
#'     the variance discussion. The field name is retained for compatibility.}
#'   \item{sdev}{Square roots of \code{eigenvalues}, not generally the
#'     marginal standard deviations of correlated score columns.}
#'   \item{score_variance}{Numeric vector of marginal score squared norms
#'     divided by \eqn{n-1}; sample variances when centering is enabled.}
#'   \item{var}{Total squared energy of the processed input divided by
#'     \eqn{n-1}, including directions omitted from the model.}
#'   \item{center, scale}{Training means and scale factors, respectively,
#'     or \code{FALSE} when the corresponding operation is disabled.}
#'   \item{objective}{Numeric vector with one penalized objective value
#'     per completed optimization iteration.}
#'   \item{iterations}{Number of completed iterations.}
#'   \item{converged}{Logical indicator of objective-based convergence.
#'     This is not a global-optimality certificate.}
#'   \item{improvement}{Relative change used by the final convergence check.}
#'   \item{alpha, beta}{Effective penalty coefficients after multiplication
#'     by the squared leading singular value.}
#'   \item{penalty}{The selected sparsity penalty.}
#'   \item{feature_names}{Training column names used for prediction, or
#'     \code{NULL} for unnamed input.}
#' }
#' The \code{summary} method returns a numeric matrix with four rows:
#' adjusted variance, its square root, proportion of total input energy,
#' and cumulative proportion. The historical row labels and three-decimal
#' rounding are retained. The \code{print} method displays the fitted
#' statistics and loadings.
#'
#' @note Version 1.0.0 corrects the variance statistics: older versions used
#' singular values from an auxiliary Procrustes update, which were not the
#' variances of the returned scores. Loading signs and bases within tied
#' singular subspaces may differ between valid numerical implementations.
#'
#' @references
#' Erichson, N. B., Zheng, P., Manohar, K., Brunton, S. L., Kutz, J. N.,
#' and Aravkin, A. Y. (2020). Sparse Principal Component Analysis via
#' Variable Projection. \emph{SIAM Journal on Applied Mathematics},
#' 80(2), 977--1002. \doi{10.1137/18M1211350}.
#'
#' Zou, H., Hastie, T., and Tibshirani, R. (2006). Sparse Principal
#' Component Analysis. \emph{Journal of Computational and Graphical
#' Statistics}, 15(2), 265--286. \doi{10.1198/106186006X113430}.
#'
#' @seealso \code{\link{rspca}} for randomized fitting,
#'   \code{\link{robspca}} for sparse additive errors,
#'   \code{\link{spca_path}} for a penalty path,
#'   \code{\link{predict.spca}} for new observations,
#'   \code{\link[stats]{prcomp}} for ordinary PCA.
#'
#' @examples
#' X <- as.matrix(iris[, 1:4])
#' fit <- spca(X, k = 2, alpha = 0.01, beta = 0.001,
#'             scale = TRUE, verbose = FALSE)
#' colSums(fit$loadings != 0)
#' summary(fit)
#'
#' # Reconstruct in original measurement coordinates.
#' reconstructed <- tcrossprod(fit$scores, fit$transform)
#' reconstructed <- sweep(reconstructed, 2, fit$scale, "*")
#' reconstructed <- sweep(reconstructed, 2, fit$center, "+")
#'
#' # Hard thresholding with a ridge penalty (ZeroNet).
#' hard <- spca(X, k = 2, alpha = 0.005, beta = 0.001,
#'              penalty = "l0", scale = TRUE, verbose = FALSE)
#' colSums(hard$loadings != 0)
#'
#' # Ordinary PCA with matching preprocessing is the unpenalized case.
#' ordinary <- spca(X, k = 2, alpha = 0, beta = 0,
#'                  scale = TRUE, verbose = FALSE)
#' reference <- stats::prcomp(X, rank. = 2, scale. = TRUE)
#' cbind(ordinary$sdev, reference$sdev[1:2])
NULL

#' Randomized sparse principal component analysis
#'
#' Approximate a sparse principal component fit by compressing the rows
#' into a randomized subspace before optimizing the loading coefficients.
#' The final component scores use the complete processed data.
#'
#' @name rspca
#' @aliases rspca.default
#' @inheritParams spca
#' @param o Nonnegative integer oversampling parameter (default \code{20}).
#'   The sketch width is \code{min(k + o, nrow(X), ncol(X))} after
#'   incomplete rows are omitted. Additional columns can improve the
#'   approximation while increasing memory and computation.
#' @param q Nonnegative integer number of power iterations (default
#'   \code{2}). Each iteration requires additional products with the
#'   processed input and its transpose, with intermediate orthogonalization.
#'   Increasing \code{q} can help for a slowly decaying singular spectrum.
#' @param X A real numeric matrix, numeric data frame, or real
#'   \pkg{Matrix} object. Observations are rows and variables are columns.
#'   Sparse \pkg{Matrix} inputs retain sparse storage for matrix products
#'   and implicit centering/scaling. Missing rows are omitted with a warning;
#'   at least two complete rows are required. Infinite and complex values
#'   are unsupported.
#'
#' @details
#' The algorithm forms an orthonormal basis \eqn{Q} of width at most
#' \eqn{k+o}, using random combinations of the processed input columns.
#' Power iterations refine this basis. It then replaces \eqn{X} with
#' \eqn{C=Q^T X} in the sparse fitting objective described by
#' \code{\link{spca}}. The coefficients still correspond to the original
#' variables. Returned scores are \eqn{XB}, and reconstruction is
#' \eqn{XBA^T} in processed coordinates.
#'
#' Both relative penalties are scaled by the squared leading singular value
#' of the sketch \eqn{C}. The \code{objective} history refers to the
#' sketched optimization problem, not the full-data reconstruction loss.
#' Variance statistics are calculated from the full-data scores using the
#' conventions documented for \code{\link{spca}}. Use full-data residuals
#' or another common reference when comparing randomized and deterministic
#' fit quality.
#'
#' Randomized compression is most useful when \eqn{k+o} is much smaller
#' than the matrix dimensions and the data have a useful low-rank
#' approximation. Small matrices, a large requested rank, or a flat singular
#' spectrum reduce the advantage. A full-rank noisy matrix can still admit
#' an accurate low-rank approximation; the relevant quantity is the energy
#' outside the retained subspace. Omitting \code{k} requests the full
#' available component count and offers no compression advantage.
#'
#' Call \code{set.seed} before fitting to reproduce random draws. The
#' function uses R's current random-number stream and does not reset it.
#' Numerical libraries can still choose different signs or bases for tied
#' subspaces. A warm start reuses loading coefficients but a separate call
#' creates its own sketch; \code{\link{spca_path}} reuses one sketch across
#' the whole path.
#'
#' @section Sparse input and memory:
#' For sparse input, products with the centered/scaled matrix are evaluated
#' without creating the entire centered dense matrix. For example,
#' \eqn{(X-\mathbf{1}\mu^T)W=XW-\mathbf{1}(\mu^T W)}. This preserves
#' sparse storage of \eqn{X}; the mean correction uses smaller dense
#' factors. The basis \eqn{Q}, projected matrix \eqn{C}, loading factors,
#' and score matrix are generally dense, with dimensions \eqn{n} by
#' \eqn{\ell}, \eqn{\ell} by \eqn{p}, \eqn{p} by \eqn{k}, and
#' \eqn{n} by \eqn{k}, respectively, where
#' \eqn{\ell=\min(k+o,n,p)}. A full reconstruction is dense.
#'
#' Sparse input and sparse loading coefficients are separate concepts:
#' \code{rspca} supports both, but zero coefficients in the returned
#' loadings are stored in an ordinary matrix. See \code{\link{spca}} and
#' \code{\link{robspca}} for methods that use dense input workspace.
#'
#' @return A list of class \code{"rspca"} with the same components and
#'   dimensions as a \code{\link{spca}} result. \code{objective},
#'   effective \code{alpha}, and effective \code{beta} refer to the sketch;
#'   \code{scores}, \code{var}, \code{score_variance},
#'   \code{eigenvalues}, and \code{sdev} use the full processed data.
#'   The \code{print} and \code{summary} methods follow the corresponding
#'   \code{spca} methods.
#'
#' @references
#' Erichson, N. B., Zheng, P., Manohar, K., Brunton, S. L., Kutz, J. N.,
#' and Aravkin, A. Y. (2020). Sparse Principal Component Analysis via
#' Variable Projection. \emph{SIAM Journal on Applied Mathematics},
#' 80(2), 977--1002. \doi{10.1137/18M1211350}.
#'
#' Erichson, N. B., Voronin, S., Brunton, S. L., and Kutz, J. N. (2019).
#' Randomized Matrix Decompositions Using R. \emph{Journal of Statistical
#' Software}, 89(11), 1--48. \doi{10.18637/jss.v089.i11}.
#'
#' @seealso \code{\link{spca}}, \code{\link{robspca}},
#'   \code{\link{spca_path}}, \code{\link{predict.spca}}.
#'
#' @examples
#' # A small example illustrating sparse storage and sparse coefficients.
#' X <- Matrix::bdiag(outer(1:20, 1:6), outer(1:20, 6:1),
#'                    outer(1:20, 1:6))
#' set.seed(42)
#' fit <- rspca(X, k = 3, alpha = 0.005, scale = TRUE,
#'              o = 5, q = 2, verbose = FALSE)
#' colSums(fit$loadings != 0)
#' summary(fit)
#' predict(fit, X[1:3, , drop = FALSE])
NULL

#' Robust sparse principal component analysis
#'
#' Fit sparse component loadings and a sparse additive error matrix to
#' reduce the influence of large individual reconstruction residuals.
#' The method combines variable projection with an entrywise Huber loss.
#'
#' @name robspca
#' @aliases robspca.default
#' @inheritParams spca
#' @param gamma Finite nonnegative penalty on the absolute entries of the
#'   sparse error matrix (default \code{100}). This is a threshold in the
#'   centered/scaled input units and is not scaled by a singular value.
#'   Smaller values allow more residual entries into the sparse component.
#'   A sufficiently large value leaves that component zero and recovers
#'   the ordinary squared reconstruction loss.
#' @param X A real numeric matrix, numeric data frame, or real
#'   \pkg{Matrix} object, with observations in rows and variables in columns.
#'   Input is converted to a dense matrix. Missing rows are omitted with
#'   a warning; at least two complete rows are required. Infinite and
#'   complex values are unsupported.
#'
#' @details
#' For processed data \eqn{X}, robust SPCA minimizes
#' \deqn{\frac{1}{2}\|X-XBA^T-S\|_F^2 +
#'       \lambda_1 P(B) + \frac{\lambda_2}{2}\|B\|_F^2 +
#'       \gamma\|S\|_1}
#' subject to \eqn{A^T A=I}. The function \eqn{P} is the entrywise
#' absolute-value sum for \code{penalty = "elasticnet"} or the number
#' of nonzero entries for \code{penalty = "l0"}. As in
#' \code{\link{spca}}, \eqn{\lambda_1} and \eqn{\lambda_2} are the
#' user-supplied \code{alpha} and \code{beta} multiplied by the squared
#' leading singular value of the processed input.
#'
#' For fixed \eqn{B} and \eqn{S}, an SVD of
#' \eqn{(X-S)^T X B} gives the orthogonal factor \eqn{A}. The algorithm
#' takes a proximal gradient step for \eqn{B}, then soft-thresholds
#' \eqn{X-XBA^T} at \code{gamma} to update \eqn{S}. The objective must
#' include \eqn{S} in the squared reconstruction residual as well as
#' charging its penalty. Eliminating \eqn{S} gives the entrywise Huber
#' loss: quadratic near zero and linear beyond the threshold.
#'
#' The predictor matrix in the score map remains the processed \eqn{X},
#' not \eqn{X-S}. Consequently \code{scores = X \%*\% loadings} and
#' the low-rank reconstruction is
#' \code{tcrossprod(scores, transform)}. Adding \code{sparse} gives the
#' full fitted model in processed coordinates. To obtain original units,
#' undo scaling for both components and add the training means only once.
#' Prediction applies the learned linear map and does not estimate a new
#' sparse-error matrix.
#'
#' Robustness here concerns an additive sparse-error model for individual
#' entries. It is not a guarantee against every form of contamination.
#' Centering and scaling use ordinary means and sums of squares rather
#' than robust location/scale estimates. \code{gamma = 0} allows the
#' sparse component to absorb the entire remaining residual and can give
#' a degenerate low-rank fit. Choose the threshold with the processed
#' measurement units and the application in mind.
#'
#' \code{robspca} retains dense input, residual, and sparse-error matrices;
#' accepting a sparse \pkg{Matrix} does not imply sparse-memory computation.
#' The optimization is nonconvex, and \code{max_iter} remains a work limit.
#' Variance diagnostics summarize the returned full-data scores; they do
#' not measure the amount of corruption removed or certify correct
#' identification of outliers.
#'
#' @return A list of class \code{"robspca"} with the components documented
#'   for \code{\link{spca}}, plus \code{sparse}: a dense \eqn{n} by
#'   \eqn{p} matrix \eqn{S} of estimated sparse errors in processed input
#'   units. Its rows correspond to retained complete observations.
#'   \code{objective} includes the complete squared residual and all
#'   three penalties. The \code{print} and \code{summary} methods follow
#'   the corresponding \code{spca} methods.
#'
#' @note Version 1.0.0 corrects an objective bookkeeping error reported by
#'   Maria Sans-Fuentes. Earlier versions charged the penalty on \eqn{S}
#'   but omitted it from the squared reconstruction residual, which could
#'   make the recorded objective rise and terminate optimization too early.
#'   Robust fits with nonzero \code{sparse} can therefore run longer and
#'   return different fitted values after the correction.
#'
#' @references
#' Erichson, N. B., Zheng, P., Manohar, K., Brunton, S. L., Kutz, J. N.,
#' and Aravkin, A. Y. (2020). Sparse Principal Component Analysis via
#' Variable Projection. \emph{SIAM Journal on Applied Mathematics},
#' 80(2), 977--1002. \doi{10.1137/18M1211350}.
#'
#' @seealso \code{\link{spca}}, \code{\link{rspca}},
#'   \code{\link{predict.spca}}.
#'
#' @examples
#' set.seed(42)
#' signal <- matrix(rnorm(120 * 2), 120, 2)
#' weights <- matrix(0, 8, 2)
#' weights[1:4, 1] <- 1
#' weights[5:8, 2] <- 1
#' X <- tcrossprod(signal, weights) + matrix(rnorm(120 * 8, sd = 0.1), 120)
#' X[cbind(c(3, 40, 99), c(1, 6, 4))] <- c(8, -10, 9)
#' fit <- robspca(X, k = 2, alpha = 0.001, gamma = 0.5,
#'                verbose = FALSE)
#' summary(fit)
#' sum(fit$sparse != 0)
#' low_rank <- tcrossprod(fit$scores, fit$transform)
#' processed_fit <- low_rank + fit$sparse
#' # Here scale = FALSE; add the training means to restore input units.
#' reconstructed <- sweep(processed_fit, 2, fit$center, "+")
NULL

#' Display and summarize sparse PCA fits
#'
#' Print sparse loading coefficients and adjusted score statistics, or
#' return a compact variance-summary matrix for a fitted model.
#'
#' @name spca-methods
#' @aliases print.spca print.rspca print.robspca summary.spca summary.rspca summary.robspca
#' @param x,object A fitted object returned by \code{\link{spca}},
#'   \code{\link{rspca}}, or \code{\link{robspca}}.
#' @param ... Additional arguments; currently ignored.
#' @details
#' All three fit classes use the variance conventions described in
#' \code{\link{spca}}. The historical labels \code{eigenvalues} and
#' ``Standard deviations'' refer to sequentially adjusted score energies
#' and their square roots. Sparse components can be correlated, and these
#' are not generally covariance eigenvalues or marginal score standard
#' deviations. The summary divides by total processed input energy,
#' including directions omitted from the model, rather than by the sum
#' of retained values. For zero input energy it displays zero proportions
#' by convention.
#'
#' With \code{center = FALSE}, the displayed quantities use uncentered
#' squared energies divided by \eqn{n-1}. For robust fits the statistics
#' summarize the score map and do not measure the amount of outlier
#' contamination removed. Inspect \code{objective}, \code{iterations},
#' and \code{converged} for optimization diagnostics.
#'
#' @return \code{summary} returns a numeric matrix with one column per
#' component and rows named \code{"Explained variance"},
#' \code{"Standard deviations"}, \code{"Proportion of variance"}, and
#' \code{"Cumulative proportion"}, rounded to three decimal places.
#' The print methods display standard deviations, adjusted variances,
#' and sparse loadings and invisibly return the rounded loading matrix.
#' @seealso \code{\link{spca}}, \code{\link{rspca}},
#'   \code{\link{robspca}}, \code{\link{predict.spca}}.
#' @examples
#' fit <- spca(iris[, 1:4], k = 2, scale = TRUE, verbose = FALSE)
#' print(fit)
#' summary(fit)
NULL

#' Predict sparse principal component scores
#'
#' Apply a fitted sparse PCA transformation to new observations using the
#' stored training means, scales, and loading coefficients.
#'
#' @name predict.spca
#' @aliases predict.rspca predict.robspca
#' @param object A fitted object returned by \code{\link{spca}},
#'   \code{\link{rspca}}, or \code{\link{robspca}}.
#' @param newdata A real numeric matrix, numeric data frame, or real
#'   \pkg{Matrix} object with the same variables as the training data.
#'   Missing and nonfinite values are not allowed. If both training and
#'   new data have column names, names must be unique and identify the
#'   exact same variable set; columns are reordered to training order.
#'   Otherwise supply columns in training order. If omitted, return the
#'   stored training scores.
#' @param ... Additional arguments; currently ignored.
#'
#' @details
#' Prediction first applies the stored column means and scale factors,
#' skipping operations stored as \code{FALSE}, and then multiplies by
#' \code{object$loadings}. It does not estimate new preprocessing
#' parameters. Centering sparse input can be applied as a correction to
#' the score product; the returned score matrix is dense.
#'
#' Predictions for robust fits use this same linear score map. The
#' method does not infer new outliers or subtract the training
#' \code{sparse} matrix from new observations. The stored training scores
#' and new scores are in component coordinates, not original variable
#' coordinates. Multiplication by \code{t(object$transform)} reconstructs
#' the low-rank part in processed coordinates; undo column scaling and
#' centering to restore original units.
#'
#' @return A dense numeric matrix with one row per new observation and
#'   one column per fitted component. With no \code{newdata}, returns
#'   \code{object$scores}, including its retained training-row order.
#'
#' @seealso \code{\link{spca}}, \code{\link{rspca}},
#'   \code{\link{robspca}}.
#'
#' @examples
#' fit <- spca(iris[1:140, 1:4], k = 2, alpha = 0.01,
#'             scale = TRUE, verbose = FALSE)
#' scores <- predict(fit, iris[141:150, 1:4])
#' head(scores)
#' # A complete named variable set can be supplied in a different order.
#' reordered <- predict(fit, iris[141:150, 4:1])
#' all.equal(scores, reordered)
#' identical(predict(fit), fit$scores)
NULL

#' Sparse PCA regularization paths
#'
#' Fit a sequence of sparsity penalties while reusing preprocessing and
#' the initial matrix factorization or randomized sketch. Each fit provides
#' the initial loading coefficients for the next penalty in the supplied
#' order.
#'
#' @name spca_path
#' @aliases print.spca_path
#' @inheritParams spca
#' @inheritParams rspca
#' @param alpha A nonempty numeric vector of finite nonnegative relative
#'   sparsity penalties. Fits follow this order without sorting. Increasing
#'   penalties from a weakly regularized fit is a convenient starting point.
#' @param method Fit with deterministic \code{"spca"} (the default) or
#'   randomized \code{"rspca"}. All randomized fits share one sketch.
#'   \code{o} and \code{q} are used only by the randomized method.
#' @param B_init Optional \eqn{p} by \eqn{k} loading matrix for the first
#'   fit. Later fits initialize from the preceding fit's raw loadings.
#' @param x An object returned by \code{spca_path}.
#' @param ... Additional arguments; currently ignored by the print method.
#'
#' @details
#' A path saves the preprocessing and factorization work that separate
#' calls would repeat. Warm starts may also reduce optimization iterations,
#' but this is not guaranteed. The objective is nonconvex, so a path fit
#' can differ from a fresh SVD initialization at the same penalty.
#'
#' Penalties encourage sparse coefficients; they do not impose an exact
#' count. Counts can jump or be nonmonotone. \code{nonzero} reports
#' the number of active coefficients in each component, whereas
#' \code{selected_features} counts variables active in at least one
#' component. These answer different questions when components overlap.
#' Inspect reconstruction quality and stability in addition to counts.
#'
#' Centering, scaling, the component count, the ridge penalty, and the
#' chosen penalty type are shared across fits. For randomized paths,
#' \code{set.seed} controls the single initial sketch. Each fit uses the
#' result class, variance conventions, and convergence diagnostics of the
#' selected method. Predict with an individual element of \code{fits}.
#'
#' @return A list of class \code{"spca_path"} with:
#' \describe{
#'   \item{fits}{A list of fitted \code{spca} or \code{rspca} objects in
#'     the supplied penalty order.}
#'   \item{alpha}{The user-supplied penalties before singular-value scaling.}
#'   \item{nonzero}{An integer matrix with one row per penalty and one
#'     column per component, containing \code{colSums(loadings != 0)}.}
#'   \item{selected_features}{The number of variables with a nonzero
#'     loading in at least one component, for each fit.}
#'   \item{iterations}{The number of optimization iterations for each fit.}
#'   \item{method, penalty}{The selected fitting method and sparsity penalty.}
#' }
#' The print method displays a compact table of penalties, selected
#' variable counts, iteration counts, and convergence indicators;
#' per-component nonzero counts are available in \code{nonzero}.
#'
#' @seealso \code{\link{spca}}, \code{\link{rspca}},
#'   \code{\link{predict.spca}}.
#'
#' @examples
#' path <- spca_path(iris[, 1:4], k = 2,
#'                   alpha = c(0.001, 0.005, 0.01, 0.03),
#'                   scale = TRUE, verbose = FALSE)
#' print(path)
#' path$nonzero
#' predict(path$fits[[2]], iris[1:3, 1:4])
#'
#' set.seed(42)
#' randomized <- spca_path(iris[, 1:4], k = 2,
#'                         alpha = c(0.001, 0.01),
#'                         method = "rspca", scale = TRUE,
#'                         verbose = FALSE)
NULL

#' Sparse principal component analysis via variable projection
#'
#' Deterministic, randomized, and robust sparse PCA with loading penalties,
#' prediction, and shared-factorization regularization paths.
#'
#' @name sparsepca-package
#' @aliases sparsepca
#' @details
#' The package separates the sparse coefficients that create component
#' scores from the orthogonal factor that reconstructs the observations.
#' \code{\link{spca}} gives the mathematical model, parameter scaling,
#' variance conventions, and complete returned-object description.
#' \code{\link{rspca}} adds randomized compression and sparse input
#' products; \code{\link{robspca}} models sparse additive errors.
#' \code{\link{spca_path}} shares preprocessing across penalties, and
#' \code{\link{predict.spca}} projects new observations.
#'
#' Sparse loadings use a subset of variables. Sparse input storage is a
#' separate property of a matrix; only the randomized fitting path
#' preserves it during centering and scaling. All fits return dense
#' loading and score matrices. The optimization is nonconvex, and
#' sparsity does not itself guarantee stability or generalization.
#'
#' @references
#' Erichson, N. B., Zheng, P., Manohar, K., Brunton, S. L., Kutz, J. N.,
#' and Aravkin, A. Y. (2020). Sparse Principal Component Analysis via
#' Variable Projection. \emph{SIAM Journal on Applied Mathematics},
#' 80(2), 977--1002. \doi{10.1137/18M1211350}.
#' @author N. Benjamin Erichson, Peng Zheng, and Sasha Aravkin.
"_PACKAGE"
