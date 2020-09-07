#' @title Function to convert p-values to correlation via the t-statistic
#' @param p Matrix of p-values
#' @param n Sample size
#' @param traits The number of traits
#' @param sign A matrix giving the sign of the correlations (e.g. the log odds ratios)
#' @param min.n The minimum sample size to be considered a valid p-value
#' @return A matrix of correlations
#' @export
p2cor <- function(p, n, sign=rep(1, length(p)), min.n=max(n, 30)/10) {
  
  stopifnot(length(n)==1 || length(n) == length(p))
  stopifnot(length(p) == length(sign))
  
  t <- sign(sign) * qt(p/2, df=n-2, lower.tail=F)
  invalid <- n < min.n
  if(any(invalid)) {
    warning(paste(sum(invalid), "statistics has n < ", min.n, "and are coded as NA."))
  }
  t[invalid] <- NA
  
  return(t / sqrt(n - 2 + t^2))
  
}

