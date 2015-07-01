#' Implements a Kalman for dynamic factor model.
#'
#' @param initx Initial value for state space observations
#' @param initV Initial value for state covariance
#' @param x Observation matrix
#' @param A State space matrix
#' @param C System matrix
#' @param R State space covariance
#' @param Q System covariance
#' @return Filtered state space variable and its covariance matrix
#' as well as their forecast for next period for further iterations
K_filter <- function(initx, initV, x, A, C, R, Q) {
    
    T <- dim(x)[1]
    N <- dim(x)[2]
    r <- dim(A)[1]
    W <- !is.na(x)    
    y <- t(x)
    
    xittm <- matrix(0, r, (T+1))
    xitt <- matrix(0, r, T)
    
    Pttm <- array(0, c(r, r, (T+1)))
    Ptt <- array(0, c(r, r, T))

    xittm[,1] <- initx
    Pttm[,,1] <- initV
   
    logl <- c()
    Ci <- C
    Ri <- R
    for (j in 1:T) {
#      missing_data <- MissData(y[,j], C, R)
#      C <- missing_data$C
#      R <- missing_data$R
      C <- Ci[W[j,],, drop=FALSE]
      R <- Ri[W[j,], W[j,], drop=FALSE]
      if (FALSE) #(all(!W[j,])) #(all(is.na(missing_data$y) == TRUE))
      {
         xitt[,,j] <- A %*% xittm[,,j]
         Ptt[,,j] <- C %*% Pttm[,,j] %*% t(C) + R
      } else
      {
         # Innovation covariance (inverse)
         Icov <- C %*% Pttm[,,j] %*% t(C) + R
         L <- solve(Icov)
         # Innovation residual
         Ires <- as.numeric(na.omit(y[,j])) - C %*% xittm[,j]
         # Optimal Kalman gain
         G <- Pttm[,,j] %*% t(C) %*% L
         # Updated state estimate: predicted + (Kalman gain)*fitted
         xitt[,j] <- xittm[,j] + G %*% Ires
         # Updated covariance estimate
         Ptt[,,j] <- Pttm[,,j] - G %*% C %*% Pttm[,,j]
         # State space variable and covariance predictions E[f_t | t-1]
         xittm[,(j+1)] <- A %*% xitt[,j]
         Pttm[,,(j+1)] <- A %*% Ptt[,,j] %*% t(A) + Q

         # Compute log-likelihood with Mahalanobis distance
         d <- length(Ires)
         S <- C %*% Pttm[,,j] %*% t(C) + R
         Sinv <- solve(S)
         if (nrow(R) == 1)
         {
           GG <- t(C) %*% solve(R) %*% C
           detS <- prod(R) %*% det(diag(1, r) + Pttm[,,j] %*% GG)
         } else {
           GG <- t(C) %*% diag(1/diag(R)) %*% C
           detS <- prod(diag(R)) * det(diag(1, r) + Pttm[,,j] %*% GG)
         }
         denom <- (2 * pi)^(d/2) * sqrt(abs(detS))
         mahal <- sum(t(Ires) %*% Sinv %*% Ires)
         logl[j] <- -0.5 * mahal - log(denom)
        }
    }
    loglik <- sum(logl, na.rm=TRUE)
    return(list(xitt = xitt, xittm = xittm, Ptt = Ptt, Pttm = Pttm, loglik = loglik))
}

#' Implements Kalman smoothing and is used along with Kalman filter.
#' Kalman filter outputs enter Kalman smoother as inputs.
#'
#' @param A State space matrix
#' @param xitt State space variable
#' @param xittm Predicted state space variable
#' @param Ptt State space covariance
#' @param Pttm Predicted state space covariance
#' @param C System matrix
#' @param R State space covariance
#' @param W Logical matrix (T x n) indicating missing data.
#' TRUE if observation is present, FALSE if it is missing.
#' @return Smoothed state space variable and state space covariance matrix
K_smoother <- function(A, xitt, xittm, Ptt, Pttm, C, R, W) {
    T <- dim(xitt)[2]
    r <- dim(A)[1]
    
    Pttm <- Pttm[,,(1:(dim(Pttm)[3] - 1)), drop = FALSE]
    xittm <- xittm[,(1:(dim(xittm)[2] - 1)), drop = FALSE]
   
    # Whereas J is of constant dimension, L and K dimensions may vary
    # depending on existence of NAs
    J <- array(0, c(r, r, T))
    L <- list()
    K <- list()

    for (i in 1:(T-1)) {
        J[,,i] <- Ptt[,,i] %*% t(A) %*% solve(Pttm[,,(i+1)], tol = 1e-32)
    }

    Ci <- C 
    Ri <- R
    for (i in 1:T) {
      # Only keep entries for non-missing data
      C <- Ci[W[i,],, drop=FALSE]
      R <- Ri[W[i,], W[i,], drop=FALSE]
      L[[i]] <- solve(C %*% Pttm[,,i] %*% t(C) + R)
      K[[i]] <- Pttm[,,i] %*% t(C) %*% L[[i]]
    }
    
    xitT <- cbind(matrix(0, r, (T-1)), xitt[,T])
    PtT <- array(0, c(r, r, T))
    PtTm <- array(0, c(r, r, T))
    PtT[,,T] <- Ptt[,,T]
    PtTm[,,T] <- (diag(1, r) - K[[T]] %*% C) %*% A %*% Ptt[,,(T-1)]

    for (j in 1:(T-1)) {
        xitT[,(T-j)] <- xitt[,(T-j)] + J[,,(T-j)] %*% (xitT[,(T+1-j)] - xittm[,(T+1-j)])
        PtT[,,(T-j)] <- Ptt[,,(T-j)] + J[,,(T-j)] %*% (PtT[,,(T+1-j)] - Pttm[,,(T+1-j)]) %*% t(J[,,(T-j)])
    }
    
    for (j in 1:(T-2)) {
        PtTm[,,(T-j)] <- Ptt[,,(T-j)] %*% t(J[,,(T-j-1)]) + J[,,(T-j)] %*% (PtTm[,,(T-j+1)] - A %*% Ptt[,,(T-j)]) %*% t(J[,,(T-j-1)])
    }

    return(list(xitT = xitT, PtT = PtT, PtTm = PtTm))
}
