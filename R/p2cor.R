#' @title Function to convert p-values to correlation via the t-statistic
#' @param p Matrix of p-values
#' @param n Sample size
#' @param traits The number of traits
#' @param sign A matrix giving the sign of the correlations (e.g. the log odds ratios)
#' @param min.n The minimum sample size to be considered a valid p-value
#' @return A matrix of correlations
#' @export
p2cor_ct <- function(p, n, traits, sign, min.n=max(n, 30)/10) {
  # if(traits==1){
  #   p <- as.matrix(p)
  #   sign <- as.matrix(sign)
  # }
  stopifnot(length(n)==1 || length(n) == nrow(p))
  stopifnot(nrow(p) == nrow(sign))
  
  t <- lapply(1:traits, function(i){
    sign(sign[,i,with=FALSE]) * qt(as.matrix(p)[,i]/2, df=n-2, lower.tail=F)
  })
  t <- matrix(unlist(t),byrow = F, ncol=traits)
  
  invalid <- n < min.n
  if(any(invalid)) {
    warning(paste(sum(invalid), "statistics has n < ", min.n, "and are coded as NA."))
  }
  t[invalid,] <- NA
  
  return(t / sqrt(n - 2 + t^2))
  
}
