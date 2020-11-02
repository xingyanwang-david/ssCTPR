#' @title Subset a ssCTPR.pipeline object by lambda, s and lambda_ct
#'
#' @param ssCTPR.pipeline An object returned by ssCTPR.pipeline
#' @param s Value(s) of s to restrict to 
#' @param lambda Value(s) of lambda to restrict to 
#' @param lambda_ct Values(s) of lambda_ct to restrict to
#' @details This function is usually used to reapply a validated pgs to a new data.set. 
#' See example below. 
#' 
#' @rdname subset.ssCTPR.pipeline
#' @export
subset.ssCTPR.pipeline <- function(ssCTPR.pipeline, s=NULL, lambda=NULL, lambda_ct=NULL) {
  
  err <- function(param, value) {
    stop(paste("There is no", param, "equalling", value, "in ssCTPR.pipeline"))
  }
  lp <- ssCTPR.pipeline
  
  
  if(!is.null(lambda_ct)) {
    w <- which(lp$lambda_ct %in% lambda_ct)
    if(length(w) == 0) err("lambda_ct", lambda_ct)
    lp$lambda_ct <- lp$lambda_ct[w]
    lp$beta <- lp$beta[w]
    lp$pgs <- lp$pgs[w]
  }
  
  if(!is.null(s)) {
    w <- which(lp$s %in% s)
    if(length(w) == 0) err("s", s)
    lp$s <- lp$s[w]
    for(i in 1:length(lp$lambda_ct)){
      lp$beta[[i]] <- lp$beta[[i]][w]
      lp$pgs[[i]] <- lp$pgs[[i]][w] 
    }
  }
  
  if(!is.null(lambda)) {
    w <- which(lp$lambda %in% lambda)
    if(length(w) == 0) err("lambda", lambda)
    lp$lambda <- lp$lambda[w]
    for(j in 1:length(lp$lambda_ct)){
      for(i in 1:length(lp$s)) {
        lp$beta[[j]][[i]] <- lp$beta[[j]][[i]][,w, drop=F]
        lp$pgs[[j]][[i]] <- lp$pgs[[j]][[i]][,w, drop=F]
      } 
    }
  }

  #' @return A ssCTPR.pipeline object
  class(lp) <- "ssCTPR.pipeline"
  return(lp)
  
  #' @examples 
  #' \dontrun{
  #'  ### Run ssCTPR using standard pipeline ### 
  #'  lp <- ssCTPR.pipeline(cor=cor, traits=ncol(cor), lambda_ct = lambda_ct,
  #'                           chr=ss$Chr, pos=ss$Position, 
  #'                           A1=ss$A1, A2=ss$A2,
  #'                           ref.bfile=ref.bfile, test.bfile=test.bfile, 
  #'                           LDblocks = ld)
  #'  v <- validate(lp)
  #'  lp2 <- subset(lp, s=v$best.s, lambda=v$best.lambda, lambda_ct=v$best.lambda_ct)
  #'  v2 <- validate(lp2)
  #' }
}
