#' \code{main} main function 
#' 
#' @param T: number of days
#' @param activity: activity index
#' @param pixels: normalized flux
#' @param lambda.candidate: a set of lambda values for the lasso tuning parameter
#' @param v0.candidate: a vector of 1/variance values for spike density 
#' @param v1: 1/variance value for slab density
#' @param counts.per.day: number of observation in each night
#' @param missing.period: missing day with no observations
#' @param parallel: parallel computing, default is TRUE
#' @param fold: number of fold for cv
#' @param num.cores: number of cores for parallel computing, if input na, the maximum of cores-1 will be used
#' @param max.iter.step1: maximum iteration for selecting intial values
#' @param max.iter.step2: maximum iteration for cross validation
#' @param max.iter.step3: maximum iteration model fitting
#' @param epsilon: threshold decides when the algorithm stops
#' @param seed: set seed for replication, default is 2019
#' @param return: model fitting results

model.fitting <- 
  function(T, activity, pixels, lambda.candidate, v0.candidate, counts.per.day, 
           missing.period, parallel = TRUE, num.cores = NA, v1 = 1, 
           fold = 10, seed = 2019, max.iter.step1 = 30, max.iter.step2 = 20, 
           max.iter.step3 = 50, epsilon = 1) {
  
  y <- activity
  X <- pixels
  
  # check if a package has been installed
  pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    } 
  }
  pkgTest("parallel")
  pkgTest("glmnet")
  pkgTest("psych")
  pkgTest("rmutil")
  pkgTest("Matrix")
  
  # load packages
  require(parallel)
  require(glmnet)
  require(psych)
  require(rmutil)
  require(Matrix)
  
  # load functions
  source("code/EMVS.SSLASSO.R")
  source("code/setup.intial.value.R")
  source("code/cv.R")
  
  ######################################################
  ##################### Step 1 #########################
  ######################################################
  cat("\nStep 1: choose initial value \n")     # report progress
  
  if (length(lambda.candidate) == 1) {
    lambda.select <- lambda.candidate
    
  } else {
    # get intial values
    
    parallel <- TRUE
    num.cores <- 7
    v1 <- 1
    max.iter.step1 <- 2
    epsilon <- 1
    
    logpost.lambda <- 
      setup.initial.value(T, activity = y, pixels = X, lambda.candidate, 
                          counts.per.day, missing.period,
                          parallel = parallel, num.cores = num.cores, v0 = 100,
                          v1 = v1, max.iter = max.iter.step1, epsilon = epsilon)
    
    # choose which lambda to use
    lambda.select <- lambda.candidate[which(logpost.lambda == max(logpost.lambda))]
    print(logpost.lambda)
  }
  print(lambda.select)
  
  cat("\nStep 1 finished \n")
  
  ######################################################
  ##################### Step 2 #########################
  ######################################################
  
  cat("\nStep 2: choose v0 \n")     # report progress
  
  if (length(v0.candidate) == 1) {
    v0.select <- v0.candidate
    
  } else {
    
    # 10-fold cross-validation
    pred.error <- cv(T, activity = y, pixels = X, v0.candidate, counts.per.day, 
                     missing.period, lambda.select, parallel = parallel,
                     fold = fold, num.cores = num.cores, v1 =v1, 
                     max.iter = max.iter.step2, epsilon = epsilon, seed = seed)
    
    # select the v0 which has the smallest pred.error
    pred.error.sum <- rowSums(pred.error)
    print(pred.error.sum)
    v0.select <- v0.candidate[pred.error.sum == min(pred.error.sum)]
    if (length(v0.select) > 1) {
      v0.select <- v0.select[1]
    }
  }
  
  print(v0.select)
  
  cat("\nStep 2 finished \n")
  
  ######################################################
  ##################### Step 3 #########################
  ######################################################
  
  cat("\nStep 3: model fitting \n") 
  
  # Model fitting
  y.unlist <- unlist(y)
  y.unlist <- y.unlist[-which(is.na(y.unlist))]
  lasso.res <- glmnet(x = X, y = y.unlist, lambda = lambda.select)
  beta.int <- c(lasso.res$beta[, 1])
  beta.int <- unname(beta.int) 
  model.fitting <- EMVS.SSLASSO(T, y = y, X = X, v0 = v0.select, v1 = v1,
                                max.iter = max.iter.step3, beta.int = beta.int,
                                counts.per.day = counts.per.day,
                                missing.period = missing.period)
  beta.hat <- model.fitting$beta.hat
  pstar.hat <- model.fitting$pstar
  logpost <- model.fitting$log.post
  sigma.sq_t.res <- model.fitting$sigma.sq
  residual <- model.fitting$residual
  
  return(list(beta.hat = beta.hat, pstar.hat = pstar.hat, logpost = logpost,
              sigma.sq_t.res = sigma.sq_t.res, residual = residual, 
              logpost.lambda = logpost.lambda, pred.error = pred.error.sum))
  
}
  
  