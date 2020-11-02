merge.ssCTPR <- function(...) {
  
  #' @title Merge ssCTPR results 
  #' @description e.g. when calculated over different blocks/chromosomes
  #' @param ... \code{ssCTPR.pipeline} objects, separated by commas
  #' @method merge ssCTPR
  #' @export
  #' 
  ll <- list(...)
  stopifnot(all(sapply(ll, "class") == "ssCTPR"))
  for (ii in 1:length(ll)) {
    stopifnot(all(sapply(ll[[ii]], function(x) all(x$lambda == ll[[1]][[1]]$lambda))))
    shrink <- sapply(ll[[ii]], function(x) x$shrink)
    stopifnot(all(shrink == shrink[1]))
  }
  
  Cumsum <- function(...) {
    mat <- do.call("rbind", list(...))
    if(ncol(mat) > 0) return(as.vector(colSums(mat))) else
      return(numeric(0))
  }
  results <- list()
  len <- length(ll[[1]])
  for(ii in 1:len){
    results[[as.character(ii)]]$lambda <- ll[[1]][[1]]$lambda
    results[[as.character(ii)]]$beta <- do.call("rbind", lapply(ll, function(x) x[[ii]]$beta))
    results[[as.character(ii)]]$conv <- do.call("pmin", lapply(ll, function(x) x[[ii]]$conv))
    pred <- do.call("Cumsum", lapply(ll, function(x) as.vector(x[[ii]]$pred)))
    results[[as.character(ii)]]$pred <- matrix(pred, ncol=length(results[[ii]]$lambda), nrow=nrow(ll[[1]][[1]]$pred))
    results[[as.character(ii)]]$loss <- do.call("Cumsum", lapply(ll, function(x) x[[ii]]$loss))
    results[[as.character(ii)]]$fbeta <- do.call("Cumsum", lapply(ll, function(x) x[[ii]]$fbeta))
    results[[as.character(ii)]]$sd <- do.call("c", lapply(ll, function(x) x[[ii]]$sd))
    results[[as.character(ii)]]$shrink <- ll[[1]][[1]]$shrink
    results[[as.character(ii)]]$nparams <- do.call("Cumsum", lapply(ll, function(x) x[[ii]]$nparams))
  }
  names(results) <- names(ll[[1]])
  class(results) <- "ssCTPR"
  return(results)
}