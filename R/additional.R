#' Estimate a p-th order vector autoregressive (VAR) model
#'
#' @param x Data matrix (T x n)
#' @param p Maximum lag order, i.e. VAR(p) will be estimated
#' @return Estimated parameter matrix, residuals and regression model
#' independent and dependent variables
#' @examples
#' x <- matrix(rnorm(50*2), nrow=50, ncol=2)
#' VAR(x, 2)
VAR <- function(x, p) {
    T <- nrow(x)
    Y <- x[(p + 1):T, ]
    X <- c()
    for (i in 1:p) {
        X <- cbind(X, x[(p + 1 - i):(T - i), ])
    }
    A <- solve(t(X) %*% X) %*% t(X) %*% Y
    res <- Y - X %*% A
    
    return(list(Y = Y, X = X, A = A, res = res))
}
