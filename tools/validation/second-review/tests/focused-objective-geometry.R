test_that("the penalized update agrees with a finite-difference value-function gradient",
    {
        set.seed(802)
        X <- matrix(rnorm(18 * 6), 18, 6)
        X[, 2:6] <- X[, 2:6] + 0.4 * X[, 1:5]
        B <- matrix(rnorm(12), 6, 2)/3
        value <- function(B) {
            XB <- X %*% B
            0.5 * sum(X^2) + 0.5 * sum(XB^2) - sum(svd(crossprod(X,
                XB), nu = 0, nv = 0)$d)
        }
        h <- 1e-05
        gradient <- B
        for (j in seq_along(B)) {
            plus <- minus <- B
            plus[j] <- plus[j] + h
            minus[j] <- minus[j] - h
            gradient[j] <- (value(plus) - value(minus))/(2 *
                h)
        }
        spectral2 <- svd(X, nu = 0, nv = 0)$d[1]^2
        alpha <- 0.035
        beta <- 0.2
        step <- 1/(spectral2 * (1 + beta))
        z <- B - step * (gradient + beta * spectral2 * B)
        for (penalty in c("elasticnet", "l0")) {
            expected <- if (penalty == "elasticnet") {
                sign(z) * pmax(abs(z) - step * alpha * spectral2,
                  0)
            }
            else {
                z * (0.5 * z^2 > step * alpha * spectral2)
            }
            fit <- spca(X, 2, alpha = alpha, beta = beta, center = FALSE,
                penalty = penalty, B_init = B, max_iter = 1,
                verbose = FALSE)
            expect_lt(max(abs(fit$loadings - expected)), 1e-07)
        }
    })
test_that("a correlated penalized fit satisfies independent optimality conditions",
    {
        set.seed(1701)
        X <- matrix(rnorm(55 * 7), 55, 7)
        X[, 2:7] <- X[, 2:7] + 0.65 * X[, 1:6]
        fit <- spca(X, 3, alpha = 0.045, beta = 0.12, max_iter = 2000,
            tol = 1e-13, verbose = FALSE)
        D <- sweep(X, 2, fit$center, "-")
        A <- fit$transform
        B <- fit$loadings
        spectral2 <- svd(D, nu = 0, nv = 0)$d[1]^2
        alpha_effective <- 0.045 * spectral2
        beta_effective <- 0.12 * spectral2
        gradient <- crossprod(D, (D %*% B) %*% t(A) - D) %*%
            A + beta_effective * B
        active <- B != 0
        residual <- gradient
        residual[active] <- gradient[active] + alpha_effective *
            sign(B[active])
        residual[!active] <- pmax(abs(gradient[!active]) - alpha_effective,
            0)
        expect_true(fit$converged)
        expect_gt(sum(active), 0)
        expect_gt(sum(!active), 0)
        expect_lt(max(abs(residual))/spectral2, 2e-06)
        cross <- crossprod(D, D %*% B)
        best_gain <- sum(svd(cross, nu = 0, nv = 0)$d)
        expect_lt(abs(best_gain - sum(A * cross))/sum(D^2), 1e-10)
    })
