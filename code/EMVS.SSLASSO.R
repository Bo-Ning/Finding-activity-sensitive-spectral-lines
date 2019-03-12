#' \code{EMVS.SSLASSO} the EM algorithm for Bayesian variable selection using SSLASSO prior
#' 
#' @param T: number of days
#' @param y: activity index
#' @param X: normalized flux
#' @param v0: 1/variance for spike density 
#' @param v1: 1/variance for slab density
#' @param beta.int: initial value for beta, if not provided, the code will generate one
#' @param counts.per.day: number of observation in each night
#' @param missing.period: missing day with no observations
#' @param lambda: lambda for the lasso tuning parameter, used for getting an initial value
#' @param max.iter: maximum iteration for EM algorithm
#' @param epsilon: threshold decides when the algorithm stops
#' @return parameters estimates from EM algorithm

EMVS.SSLASSO <- function(T, y, X, counts.per.day, missing.period, beta.int = NULL,
                         v0 = 100, v1 = 0.5, max.iter = 1000, epsilon = 1) {
  
  # load functions
  source("code/coord.descent.R")
  
  p <- dim(X)[2]
  y.unlist <- unlist(y)
  y.unlist <- y.unlist[-which(is.na(y.unlist))]
  # initial value
  if (is.null(beta.int) == TRUE) {
    lasso.res <- glmnet(x = X, y = y.unlist, lambda = 0.1)
    beta.int <- c(lasso.res$beta[, 1])
    beta.int <- unname(beta.int)  # initialize beta
  } else {
    beta.int <- beta.int
  }
  theta.int <- 0.5 # initialize theta
  sigma.sq_t.int <- rep(1, T)
  
  # prior values
  c1 <- 0.2
  d1 <- 0.1
  v0 <- v0
  v1 <- v1
  theta.a <- 1
  theta.b <- p
  
  # collect results
  # theta.res <- rep(0, max.iter)
  # sigma.sq_t.res <- matrix(0, T, max.iter)
  log.post.res <- rep(0, max.iter)
  
  ## -------------------------- Starting EMVS --------------------------- ##
  #-------------------------------- E_step --------------------------------#
  pb  <- txtProgressBar(1, max.iter, style=3)    # report progress
  cat("\nStarting EM algorithm: \n")     # report progress
  for (iter in 1:max.iter) {
   #  print(iter)
    
    # report progress
    setTxtProgressBar(pb, iter)
    
    if (iter == 1) {
      beta.update <- beta.int
      theta.update <- theta.int
      sigma.sq_t.update <- sigma.sq_t.int
      # -------------- initialize A --------------- #
      gamma0 <- dlaplace(beta.update, m = 0, s = 1/v0) # corresp. to the spike
      gamma1 <- dlaplace(beta.update, m = 0, s = 1/v1) # corresp. to the slab 
      pstar <- (gamma1*theta.update) / (gamma0*(1-theta.update) + (gamma1*theta.update))
    }
    
    #-------------------------- M_step --------------------------#
    # update beta
    y.deMean.scale <- NULL
    for (t in (1:T)[-missing.period]) {
      y.deMean.scale <- c(y.deMean.scale, y[[t]]/sqrt(sigma.sq_t.update[t]))
    }
    
    sigma.sq_ext <- NULL
    for (t in (1:T)[-missing.period]) {
      sigma.sq_ext <- c(sigma.sq_ext, rep(1, counts.per.day[t])*sqrt(sigma.sq_t.update[t]) )
    }
    X.scale <- X/sigma.sq_ext
    
    # ------------- update beta ------------- #
    # organize dataset
    ptm <- proc.time()
    beta.update <- coord.descent(y = y.deMean.scale, X = X.scale, beta.int = beta.update,
                                 v0, v1, pstar, eps = 0.1)
    proc.time() - ptm
    
    gamma0 <- dlaplace(beta.update, m = 0, s = 1/v0) # corresp. to the spike
    gamma1 <- dlaplace(beta.update, m = 0, s = 1/v1) # corresp. to the slab
    pstar <- (gamma1*theta.update) / ((gamma0*(1-theta.update)) + ((gamma1*theta.update)))
    
    # ------------- update theta ------------- #
    theta.update <- (sum(pstar) + theta.a - 1) / (theta.a + theta.b + p - 2)
    # theta.update <- 100/p
      
    # ------ update sigma.sq_t ------- #
    X.times.beta.update <- X %*% beta.update

    for (t in (1:T)[-missing.period]) {
      vector.one <- rep(1, counts.per.day[t])
      if (t == 1) {
        y.tilde.deReg <- y[[t]] - X.times.beta.update[1:counts.per.day[t]]
      } else {
        y.tilde.deReg <- y[[t]] -
          X.times.beta.update[(sum(counts.per.day[1:(t-1)])+1):sum(counts.per.day[1:t])]
      }
    
      sigma.sq.part1 <- crossprod(y.tilde.deReg)
      sigma.sq_t.update[t] <- (sigma.sq.part1/2 + d1)/(counts.per.day[t]/2 + c1 + 1)
    }
    
    # ------ calculate log-posterior distribution -------- # 
    neg.log.like <- 0
    residual <- NULL
    for (t in (1:T)[-missing.period]) {
      vector.one <- rep(1, counts.per.day[t])
      if (t == 1) {
        y.tilde <- y[[t]] - X.times.beta.update[1:counts.per.day[t]]
      } else {
        y.tilde <- y[[t]] -
          X.times.beta.update[(sum(counts.per.day[1:(t-1)])+1):sum(counts.per.day[1:t])]
      }

      residual <- c(residual, y.tilde/sqrt(sigma.sq_t.update[t]))
      neg.log.like.part1 <- sum(y.tilde^2)/sigma.sq_t.update[t]
      neg.log.like <- neg.log.like + neg.log.like.part1/2 +
        counts.per.day[t]*log(sigma.sq_t.update[t]) / 2
    }
    if (theta.update == 0) theta.update <- 1e-8
    sslasso.density <- (1-theta.update)*dlaplace(beta.update, m = 0, s = 1/v0) +
      theta.update*dlaplace(beta.update, m = 0, s = 1/v1)
    
    log.post <- -neg.log.like -
      (c1 + 1)*sum(log(sigma.sq_t.update[-missing.period])) -
      sum(d1/sigma.sq_t.update[-missing.period]) +
      sum(log(sslasso.density)) + (theta.a - 1) * log(theta.update) +
      (theta.b - 1) * log(1 - theta.update)
    
    # print(log.post)
    if (iter > 1) {
      if (abs(log.post - log.post.res[iter-1]) < epsilon)  {
        break
      }
    }
    
  }
  # return value
  list(beta.hat = beta.update, theta = theta.update, log.post = log.post,
       residual = residual, sigma.sq = sigma.sq_t.update, pstar = pstar, 
       iter = iter)
}