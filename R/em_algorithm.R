#' Computation of the expectation step in the EM-algorithm.

#' @param y Data matrix
#' @param A System state matrix
#' @param C Observation matrix
#' @param R Observation equation variance 
#' @param Q System state equation variance
#' @param W Logical matrix with dim(W) = dim(y) 
#' indicating missing observations
#' @param initx Initial value for state variable
#' @param initV Initial value for state matrix
#' @return Sufficient statistics used for M-step
Estep <- function(y, A, C, Q, R, initx, initV, W) {
    
    os <- dim(y)[1]
    T <- dim(y)[2]
    ss <- nrow(A)
    
    kf <- K_filter(initx, initV, t(y), A, C, R, Q)
    ks <- K_smoother(A, kf$xitt, kf$xittm, kf$Ptt, kf$Pttm, C, R, W)

    xsmooth <- ks$xitT
    Vsmooth <- ks$PtT
    Wsmooth <- ks$PtTm

    delta <- matrix(0, os, ss)
    gamma <- matrix(0, ss, ss)
    beta <- matrix(0, ss, ss)
    
    for (t in 1:T) {
      z <- y[,t]; z[is.na(z)] <- 0
# There seems to be a problem here
      delta <- delta + z %*% t(xsmooth[,t])
      gamma <- gamma + xsmooth[,t] %*% t(xsmooth[,t]) + Vsmooth[,,t]
      if (t > 1) {
        beta <- beta + xsmooth[,t] %*% t(xsmooth[,(t-1)]) + Wsmooth[,,t]
      }
    }
    
    gamma1 <- gamma - xsmooth[, T] %*% t(xsmooth[, T]) - Vsmooth[, , T]
    gamma2 <- gamma - xsmooth[, 1] %*% t(xsmooth[, 1]) - Vsmooth[, , 1]
    x1 <- xsmooth[, 1]
    V1 <- Vsmooth[, , 1]

    return(list(beta_t = beta, gamma_t = gamma, delta_t = delta, gamma1_t = gamma1,
                gamma2_t = gamma2, x1 = x1, V1 = V1, loglik_t = kf$loglik, xsmooth = xsmooth))
    
}

#' Convergence test for EM-algorithm.
#'
#' @param loglik Current value of the log-likelihood function
#' @param previous_loglik Value of the log-likelihood function at the previous
#  iteration
#' @param threshold If difference is less than threshold, then algorithm has
#' converged
#' @param check_increased TO DOCUMENT
#' @return A logical statement indicating whether EM algorithm has converged
#' according to slope convergence test
em_converged <- function(loglik, previous_loglik, threshold=1e-4, check_increased=TRUE) {
    
    converged <- FALSE
    decrease <- 0
    
    if (check_increased == TRUE) {
        if (loglik - previous_loglik < -0.001) {
#            cat("*** Likelihood decreased from ", previous_loglik, " to ", loglik, "\n")
            decrease <- 1
        }
    }
    
    delta_loglik <- abs(loglik - previous_loglik)
    avg_loglik <- (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2
    
    if ((delta_loglik/avg_loglik) < threshold) {
        converged <- TRUE
    }
    return(converged)
    # return(list('converged'=converged, 'decrease'=decrease))
}
