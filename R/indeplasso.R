#' @title Independent LASSO using summary statistics (a.k.a. soft-thresholding)
#' 
#' @param coef vector of regression coefficients (\eqn{r})
#' @param lambda a vector of \eqn{\lambda}s 
#' @param lambda_ct a vector of \eqn{\lambda_{ctp}}s
#' @param thr threshold to stop CD algorithm
#' @param maxiter the maximum number of iterations
#' @param trace controls the amount of output
#' 
#' @details A function to find the minimum of \eqn{\beta} in  
#' \deqn{f(\beta)=\beta'\beta - 2\beta'r + 2\lambda||\beta||_1}
#' where \eqn{r} is the vector of regression coefficients. The analytical solution
#' is given by
#' \deqn{\hat{\beta}=sign(r)(max(|r| - \lambda))}
#' @export
indeplasso_ct <- function(coef, lambda=exp(seq(log(0.001), log(0.1), length.out=20)), lambda_ct, thr=1e-4,maxiter=10000, trace=1) {
  coef <- as.matrix(coef)
  traits <- ncol(coef)
  p <- nrow(coef)
  # single trait
  #results <- outer(coef, rep(NA, length(lambda)))
  # for(i in 1:length(lambda)) {
  #   results[,i] <- sign(coef) * pmax((abs(coef) - lambda[i]),0)
  # }
  
  indeplasso_fixed_ctp <- function(lambda_ct){
    results <- matrix(0,ncol = length(lambda), nrow = length(coef))
    conv <- rep(NA,length(lambda))
    
    for(i in 1:length(lambda)) {
      for(k in 1:maxiter) {
        dlx <- 0.0
        #results[,i] <- sign(coef) * pmax((abs(coef) - lambda[i]),0)
        for (tt in 1:traits) {
          for (jj in 1:p) {
            xold <- results[(tt-1)*p+jj,i]
            results[(tt-1)*p+jj,i] <- 0.0
            ctp <- lambda_ct * sum(results[(0:(traits-1))*p+jj,i])
            if((coef[jj,tt]+ctp)>lambda[i]){
              results[(tt-1)*p+jj,i] <- (coef[jj,tt]+ctp-lambda[i])/(1+lambda_ct*(traits-1))
            } else if((coef[jj,tt]+ctp)< -lambda[i]){
              results[(tt-1)*p+jj,i] <- (coef[jj,tt]+ctp+lambda[i])/(1+lambda_ct*(traits-1))
            } else{
            }
            del <- results[(tt-1)*p+jj,i]-xold
            dlx <- max(dlx,abs(del))
          }
        }
        if(dlx<thr){
          conv[i] <- 1
          break
        }
      }
    }
    return(list(lambda=lambda, beta=results, conv=conv))
  }
  
  ls <- list()
  if(length(lambda_ct) > 0) {
    if(trace) cat("Running independent lassosum ...\n")
    ls <- lapply(lambda_ct, function(ct) {
      if(trace) cat("lambda_ct = ", ct, "\n")
      indeplasso_fixed_ctp(ct)
    })
  }
  names(ls) <- as.character(lambda_ct)
  
  #' @return A list with the length equal to the number of lambda_ct, each element of the list has teh following elements
  #' \item{lambda}{Same as \code{lambda} in input}
  #' \item{beta}{A matrix of estimates of \eqn{\beta}}

  return(ls)
}
