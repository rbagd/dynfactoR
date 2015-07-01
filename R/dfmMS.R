#' Dynamic factor model with Markov-switching states
#'
#' @param X Data matrix (Txn)
#' @param s Number of states. Only \code{2} are supported for now
#' @param J Number of factors. Currently, only \code{1} is supported
#' @param p Lag number for VAR model of factors
#' @param x0 Initial value for state
#' @param P0 Initial value for state covariance
#' @return Return an estimator. Currently, it runs via a likelihood
#' maximization an so is rather slow. Do not call it with large number
#' of input variables or missing data. Some box constraints are currently
#' to speed up the process, these are based on usual rationale for economic
#' data.
dfmMS <- function(X, J=1, s=2, p=1, x0, P0)
{
  
  T <- nrow(X)
  n <- ncol(X) 
  x <- apply(X, 2, function(z) { (z-mean(z, na.rm=TRUE))/sd(z, na.rm=TRUE) })

  J <- 1
  s <- 2

  if (missing(x0)) { x0 <- rep(0,J) }
  if (missing(P0)) { P0 <- diag(1,J) }
  
  initV <- list()
  initV$A <- c(0.96, 0.6)
  initV$F <- runif(n*J*s, -1, 1)
  initV$R <- rep(1,n)
  initV$p <- c(0.98, 0.95)

  # State equation covariance is restricted to identity
  Q <- diag(1, J)

  dimA <- J*J*s; dimF <- n*J*s; dimAF <- dimA + dimF; dimAFR <- dimAF + n

  KimFilterOptim <- function(pars)
  {
    A <- array(head(pars, dimA), c(J,J,s))
    F <- array(pars[(dimA+1):(dimAF)], c(n,J,s))
    R <- diag(pars[(dimAF+1):(dimAFR)])
    dp <- tail(pars, s); r <- 1-dp
    p <- diag(dp) + r - diag(r)
    results <- KimFilter(x0, P0, x, F, A, R, Q, p)
    return(results$result)
  }

  optimP <- optimx(par=unlist(initV),
                   fn = KimFilterOptim,
                   lower = c(rep(-0.2, J*J*s), # A
                             rep(-2.5, n*J*s), # F
                             rep(0.15, n),     # R
                             rep(0.9,s)),      # p
                   upper = c(rep(0.98, J*J*s), # A
                             rep(2.5, n*J*s),  # F
                             rep(1.5, n),      # R
                             rep(0.991,2)),    # p
                   method = "L-BFGS-B",
                   control=list(maximize=TRUE,maxit=1000, trace=1, kkt=FALSE))

  pars <- as.numeric(optimP)[1:length(unlist(initV))]
  A_hat <- array(pars[1:dimA], c(J,J,s))
  F_hat <- array(pars[(dimA+1):dimAF], c(n,J,s))
  R_hat <- diag(pars[(dimAF+1):dimAFR])
  p_hat <- tail(pars, s); p_hat <- diag(p_hat) + (1-p_hat) - diag(1-p_hat)

  kf <- KimFilter(x0, P0, x, F_hat, A_hat, R_hat, Q, p_hat)
  ks <- KimSmoother(kf$xA, kf$Pa, A_hat, kf$P, kf$x, p_hat, kf$stateP, kf$stateP_fut)

return(list("A"=A_hat, "F"=F_hat, "R"=R_hat, "p"=p_hat,
            "xF"=ks$xF, "Pf"=ks$Pf, "ProbS"=ks$ProbS))
}
