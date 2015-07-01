#' Summary information on dynamic factor model estimation
#' 
#' @param x An object of \code{dfm} class
#' @param plot If \code{TRUE}, returns a bunch of summary plots
#' @return Prints out a summary information following a dynamic factor
#' model estimation. Also can return summary plots.
summary.dfm <- function(x, plot=FALSE) {

  cl <- match.call()
  
  nf <- dim(x$qml)[2]

  if (plot == TRUE)
  {
    par(mfrow=c((nf+1),1))
    for (i in 1:nf)
    {
      plot.title <- paste0("QML estimated factor ", i)
      plot(x$qml[,i], type='l', main=plot.title, ylab="Value", xlab="Time") 
    }
    boxplot(x$data - t(x$C %*% t(x$qml)), main="Residuals by input variable")
  }
  cat("Observation equation matrix: \n")
  print(x$C)

  cat("\nObservation residual covariance matrix: \n")
  print(cov(x$data - t(x$C %*% t(x$qml))))

  cat("\nSystem equation transition matrix: \n")
  print(x$A)

}

#' Predict factors and observables based on an estimated dynamic
#' factor model
#'
#' @param x An object of \code{dfm} class
#' @param h Forecasting horizon
#' @return Returns prediction matrices for factors and input variables
#' based on model estimation. Each row corresponds to a prediction
#' horizon. First row is \code{T+1}, second -\code{T+2} and so on...
predict.dfm <- function(x, h=1) {
  
  cl <- match.call()

  nf <- dim(x$qml)[2]
  ny <- dim(x$C)[1]

  factor_forecast <- matrix(NA, nrow=h, ncol=nf)
  y_forecast <- matrix(NA, nrow=h, ncol=ny)
  factor_last <- matrix(NA, nrow=x$p, ncol=nf)

  factor_last <- tail(x$qml, x$p)

  for (i in 1:h)
  {
    reg_factor <- tail(factor_last, x$p)
    factor_forecast[i,] <- x$A %*% matrix(t(reg_factor[rev(1:nrow(reg_factor)),]))
    y_forecast[i,] <- x$C %*% factor_forecast[i,]

    factor_last <- rbind(factor_last, factor_forecast[i,])
  }
  return(list("y"=y_forecast, "f"=factor_forecast))
}
