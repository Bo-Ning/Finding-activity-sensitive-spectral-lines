#' \code{setup.initial.value} choose the optimal value for the tuning parameter of lasso
#' 
#' @param T: number of days
#' @param activity: activity index
#' @param pixels: normalized flux
#' @param lambda.candidate: a set of lambda values for the lasso tuning parameter
#' @param v0: 1/variance value for spike density
#' @param v1: 1/variance value for slab density
#' @param counts.per.day: number of observation in each night
#' @param missing.period: missing day with no observations
#' @param parallel: parallel computing, default is TRUE
#' @param fold: number of fold for cv
#' @param num.cores: number of cores for parallel computing, if input na, the maximum of cores-1 will be used
#' @param max.iter: maximum iteration for selecting intial values
#' @param epsilon: threshold decides when the algorithm stops
#' @return log posterior using different lambda for generating intial values

setup.initial.value <- 
  function(T, activity, pixels, lambda.candidate, counts.per.day, missing.period, 
           parallel = TRUE, num.cores = NA, v0 = 100, v1 = 1, max.iter = 30, 
           epsilon = 1) {

  # load data
  y <- activity
  X <- pixels
  y.unlist <- unlist(y)
  y.unlist <- y.unlist[-which(is.na(y.unlist))]
  
  if (parallel == TRUE) {
    # run parallel computing
    
    if (is.na(num.cores) == TRUE) {
      num.cores <- detectCores() - 1
    }
    
    log.post.lambda <- 
      mclapply(lambda.candidate, FUN = function(l) {
        lasso.res <- glmnet(x = X, y = y.unlist, lambda = l)
        beta.int <- c(lasso.res$beta[, 1])
        beta.int <- unname(beta.int)  # initialize beta
        res <- EMVS.SSLASSO(T = T, y = y, X = X, counts.per.day = counts.per.day, 
                            missing.period = missing.period, beta.int = beta.int, 
                            v0 = v0, v1 = v1, max.iter = max.iter, epsilon = 1)
        log.post <- res$log.post
        return(log.post)
      }, mc.cores = num.cores)
    
    log.post.lambda <- unlist(log.post.lambda)
    
  } else {
    
    # if parallel not availble
    logpost.lambda <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      lasso.res <- glmnet(x = X, y = y.unlist, lambda = lambda[i])
      beta.int <- c(lasso.res$beta[, 1])
      beta.int <- unname(beta.int)  # initialize beta
      res <- EMVS.SSLASSO(T, y, X, counts.per.day, missing.period,
                          beta.int = beta.int, v0 = v0, v1 = v1, 
                          max.iter = max.iter, epsilon = epsilon)
      logpost.lambda[i] <- res$log.post
    }
  }
  
  return(log.post.lambda)
}
