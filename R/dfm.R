#' Estimates a dynamic factor model based on Doz, Gianone & Reichlin (2011)
#'
#' @param X Data matrix (T x n)
#' @param r Number of static factors
#' @param p Lag order for factors
#' @param q Number of dynamic factors, must be equal to or less than r
#' @param max_iter Maximum number of iterations in EM-algorithm
#' @param threshold Threshold for algorithm convergence
#' @param rQ Restrictions on system state covariance
#' @param rC Restrictions on factor loading matrix
#' @return 3 types of factor estimates, namely principal component estimate, two step
#' estimate based on PCA and Kalman filtering and QML estimate based on EM-algorithm
#' @examples
#' x <- matrix(rnorm(50*10), 50, 10)
#' W <- as.logical(matrix(rbinom(50*10, 1, 0.1), 50, 10))
#' x[W] <- NA
#' dfm(x, 2, 2, 1)
dfm <- function(X, r, p, q, max_iter=100, threshold=1e-4, rQ, rC) {

  if(missing(rQ)) { rQ <- '' }
  if(missing(rC)) { rC <- '' }
  if(missing(q)) { q <- r }
  if (q > r) { stop("r must be larger than q.")}

  T <- dim(X)[1]; N <- dim(X)[2]
  x <- apply(X, 2, function(z) { (z - mean(z, na.rm=TRUE))/sd(z, na.rm=TRUE) })

  Mx <- apply(X, 2, mean, na.rm=TRUE)
  Wx <- apply(X, 2, sd, na.rm=TRUE)
  W <- !is.na(x)
  # A consists of two parts. In particular, the upper dynamic part and the
  # lower invariable identity part to ensure time lag coherence.

  A <- rbind(matrix(0, nrow=r, ncol=r*p),
             diag(1, nrow=r*(p-1), ncol=r*p))

  Q <- matrix(0, nrow=p*r, ncol=p*r)
  Q[1:r, 1:r] <- diag(1, r)

  eigen.decomp <- eigen(cov(x, use="complete.obs"))
  v <- eigen.decomp$vectors[,1:r]
  d <- eigen.decomp$values[1:r]

  chi <- x %*% v %*% t(v)
  d <- diag(1, r)
  F <- x %*% v
  F_pc <- F
  F <- na.omit(F)
  # If there are any dynamic factors, VAR(p) model is estimated
  # to initialize their parameters.

  if (p > 0)
  {
    # ML estimator for VAR(p) model when Q is restricted
    if (rQ == 'identity')
    {
      fit <- VAR(F,p)
      A[1:r, 1:(r*p)] <- t(fit$A)
      Q[1:r, 1:r] <- diag(1, r)
    }
    else
    {
      fit <- VAR(F, p)
      A[1:r, 1:(r*p)] <- t(fit$A)
      H <- cov(fit$res)

      # This only extracts the variance explained by 
      # the dynamic components

      if (r > q)
      {
        q.decomp <- eigen(H)
        P <- q.decomp$vectors[,1:q, drop=FALSE]
        M <- q.decomp$values[1:q]
        if (q == 1)
        {
          P <- P * P[1,]
          Q[1:r, 1:r] <- P %*% t(P) * M
        } else {
          P <- P %*% diag(sign(P[1,]))
          Q[1:r, 1:r] <- P %*% diag(M) %*% t(P)
        }
      } else
      {
        Q[1:r, 1:r] <- H
      }
    }
  }
  R <- diag(diag(cov(x - chi, use="complete.obs")))
  Z <- fit$X
  initx <- Z[1,]
  initV <- matrix(ginv(kronecker(A,A)) %*% as.numeric(Q),
                  ncol=r*p, nrow=r*p)
  C <- cbind(v, matrix(0, nrow=N, ncol=r*(p-1)))

  previous_loglik <- -.Machine$double.xmax
  loglik <- 0
  num_iter <- 0
  LL <- c()

  converged <- 0
  kf_res <- K_filter(initx, initV, x, A, C, R, Q)
  ks_res <- K_smoother(A, kf_res$xitt, kf_res$xittm,
                       kf_res$Ptt, kf_res$Pttm, C, R, W)

  xsmooth <- ks_res$xitT
  Vsmooth <- ks_res$PtT
  Wsmooth <- ks_res$PtTm

  F_kal <- t(xsmooth[1:r,, drop=FALSE])
  if (rC == 'upper' & (r > 1))
  {
    dimC <- dim(C[,1:r])
    rK <- rep(0, (r-1)*r/2)
    irC <- which(matrix(upper.tri(C[,1:r]) + 0) == 1)
    rH <- matrix(0, nrow=length(rK), ncol=prod(dimC))
    for (i in 1:length(rK))
    {
      rH[i,irC[i]] <- 1
    }
  }

  while ((num_iter < max_iter) & !converged)
  {

    # E-step will return a list of sufficient statistics, namely second (cross)-moments
    # for latent and observed data. This is then plugged back into M-step.
    em_res <- Estep(t(x), A, C, Q, R, initx, initV, W)
    beta <- em_res$beta_t
    gamma <- em_res$gamma_t
    delta <- em_res$delta_t
    gamma1 <- em_res$gamma1_t
    gamma2 <- em_res$gamma2_t
    P1sum <- em_res$V1 + em_res$x1 %*% t(em_res$x1)
    x1sum <- em_res$x1
    loglik <- em_res$loglik_t

    num_iter <- num_iter + 1

    # M-step computes model parameters as a function of the sufficient statistics that
    # were computed with the E-step. Iterate the procedure until convergence. Due to the
    # model specification, likelihood maximiation in the M-step is just an OLS estimation.
    # In particular, X_t = C*F_t and F_t = A*F_(t-1).

    if (rC == 'upper' & (r > 1))
    {
      fp <- matrix(delta[,1:r] %*% ginv(gamma[1:r,1:r]))
      kronCR <- kronecker(ginv(gamma[1:r,1:r]), R)
      sp <- kronCR %*% t(rH) %*% ginv(rH %*% kronCR %*% t(rH)) %*% (rK - rH %*% fp)
      C[,1:r] <- matrix(fp + sp, nrow=dimC[1], ncol=dimC[2])  
    }
    else
    {
      C[,1:r] <- delta[,1:r] %*% ginv(gamma[1:r,1:r])
    }

    if (p > 0)
    {

      A_update <- beta[1:r,1:(r*p), drop=FALSE] %*% solve(gamma1[1:(r*p),1:(r*p)])
      A[1:r,1:(r*p)] <- A_update
      if (rQ != 'identity')
      {
        H <- (gamma2[1:r, 1:r] - A_update %*% t(beta[1:r, 1:(r*p), drop=FALSE])) / (T-1)
        if (r > q)
        {
          h.decomp <- svd(H)
          P <- h.decomp$v[,1:q, drop=FALSE]
          M <- h.decomp$d[1:q]
          if (q == 1)
          {
            P <- P * P[1,]
            Q[1:r, 1:r] <- P %*% t(P) * M
          } else {
            P <- P %*% diag(sign(P[1,]))
            Q[1:r, 1:r] <- P %*% diag(M) %*% t(P)
          }
        } else {
          Q[1:r, 1:r] <- H
        }
      }
    }
    xx <- as.matrix(na.omit(x))
    R <- (t(xx) %*% xx - C %*% t(delta)) / T
    RR <- diag(R); RR[RR < 1e-7] <- 1e-7; R <- diag(RR)

    R <- diag(diag(R))
    LL <- c(LL, loglik)

    initx <- x1sum
    initV <- P1sum - initx %*% t(initx)

    converged <- em_converged(loglik, previous_loglik, threshold=threshold)
    previous_loglik <- loglik
    if (num_iter < 25) { converged <- FALSE }

  }

  if (converged == TRUE)
  {
    cat("Converged after", num_iter, "iterations.\n")
  } else
  {
    cat("Maximum number of iterations reached.\n")
  }

  kf <- K_filter(initx, initV, x, A, C, R, Q)
  ks <- K_smoother(A, kf$xitt, kf$xittm, kf$Ptt, kf$Pttm, C, R, W)

  xsmooth <- ks$xitT

  chi <- t(xsmooth) %*% t(C) %*% diag(Wx) + kronecker(matrix(1,T,1), t(Mx))
  F_hat <- t(xsmooth[1:r,, drop=FALSE])

  final_object <- list("pca"=F_pc, "qml"=F_hat, "twostep"=F_kal,
                       "A"=A[1:r,], "C"=C[,1:r], "Q"=Q[1:q,1:q], "R"=R,
                       "p"=p, "data"=x)
  class(final_object) <- c("dfm", "list")

  return(final_object)
}

#' National Bank of Belgium business and consumer surveys
#'
#' @format A time-stamped dataframe (xts object) with
#' @source \url{http://stat.nbb.be}
"NBBsurvey"
