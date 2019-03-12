#' \code{cv} cross-validation function
#' 
#' @param T: number of days
#' @param activity: activity index
#' @param pixels: normalized flux
#' @param v0.candidate: a vector of 1/variance values for spike density 
#' #' @param v1: 1/variance value for slab density
#' @param counts.per.day: number of observation in each night
#' @param missing.period: missing day with no observations
#' @param lambda: lambda for the lasso tuning parameter, used for getting an initial value
#' @param parallel: parallel computing, default is TRUE
#' @param fold: number of fold for cv
#' @param num.cores: number of cores for parallel computing, if input na, the maximum of cores-1 will be used
#' @param max.iter: maximum iteration for EM algorithm
#' @param epsilon: threshold decides when the algorithm stops
#' @param seed: set seed for replication, default is 2019
#' @return prediction errors for each v0.candidate

cv <- function(T, activity, pixels, v0.candidate, counts.per.day, 
               missing.period, lambda, parallel = TRUE, fold = 10, 
               num.cores = NA, v1 = 1, max.iter = 20, epsilon = 1, seed = 2019) {
  
  set.seed(seed)
  
  y <- activity
  X <- pixels
  num.date <- which(counts.per.day > 0)
  cv.fold <- sample(num.date, replace = F)
  len.num.date <- length(num.date)
  size <- ceiling(len.num.date/fold)
  pred.error <- matrix(0, length(v0.candidate), fold)
  
  pb  <- txtProgressBar(1, length(v0.candidate), style=3)    # report progress
  cat("\nStart cross-validation: \n")     # report progress
  
  for (j in 1:length(v0.candidate)) {
    
    setTxtProgressBar(pb, j)  # report progress
    
    if (parallel == TRUE) {
      
      # run parallel computing
      
      if (is.na(num.cores) == TRUE) {
        num.cores <- detectCores() - 1
      }
      
      pred.err <- mclapply(1:fold, FUN = function(l) {
        
        if (l == 10) {
          day.index <- cv.fold[((l-1)*size+1):len.num.date]
        } else {
          day.index <- cv.fold[((l-1)*size+1):(l*size)]
        }
        # print(index)
        missing.period.cv <- c(missing.period, day.index)
        counts.per.day.cv <- counts.per.day
        counts.per.day.cv[day.index] <- 0
        
        # test.dataset
        y.test <- y[day.index]
        y.test.unlist <- unlist(y.test)
        # fitting.dataset
        y.fitting <- y
        y.fitting[day.index] <- NA
        y.fitting.unlist <- unlist(y.fitting)
        y.fitting.unlist <- y.fitting.unlist[-which(is.na(y.fitting.unlist))]
        
        day.dataset.index <- unlist(dataset.index[day.index])
        X.test <- X[day.dataset.index, ]
        X.fitting <- X[-day.dataset.index, ]
        
        # start model fitting
        lasso.res.cv <- glmnet(x = X.fitting, y = y.fitting.unlist, lambda = lambda)
        beta.int.cv <- c(lasso.res.cv$beta[, 1])
        beta.int.cv <- unname(beta.int.cv)  # initialize beta
        model.fitting <-
          EMVS.SSLASSO(T = T, y = y.fitting, X = X.fitting,
                       counts.per.day = counts.per.day.cv, time = time,
                       missing.period = missing.period.cv, beta.int = beta.int.cv,
                       v1 = v1, v0 = v0.candidate[j], max.iter = max.iter,
                       epsilon = epsilon)
        beta.cv <- model.fitting$beta.hat
        error <- sum((y.test.unlist - X.test %*% beta.cv)^2)
        return(error)
        
      }, mc.cores = num.cores)
      pred.error[j, ] <- unlist(pred.err)
      
    } else {
      
      for (num.fold in 1:10) {
        
        if (num.fold == 10) {
          day.index <- cv.fold[((num.fold-1)*size+1):len.num.date]
        } else {
          day.index <- cv.fold[((num.fold-1)*size+1):(num.fold*size)]
        }
        # print(index)
        missing.period.cv <- c(missing.period, day.index)
        counts.per.day.cv <- counts.per.day
        counts.per.day.cv[day.index] <- 0
        
        # test.dataset
        y.test <- y[day.index]
        y.test.unlist <- unlist(y.test)
        # fitting
        y.fitting <- y
        y.fitting[day.index] <- NA
        y.fitting.unlist <- unlist(y.fitting)
        y.fitting.unlist <- y.fitting.unlist[-which(is.na(y.fitting.unlist))]
        
        day.dataset.index <- unlist(dataset.index[day.index])
        X.test <- X[day.dataset.index, ]
        X.fitting <- X[-day.dataset.index, ]
        
        # start model fitting
        lasso.res.cv <- glmnet(x = X.fitting, y = y.fitting.unlist, lambda = lambda)
        beta.int.cv <- c(lasso.res.cv$beta[, 1])
        beta.int.cv <- unname(beta.int.cv)  # initialize beta
        model.fitting <-
          EMVS.SSLASSO(T = T, y = y.fitting, X = X.fitting,
                       counts.per.day = counts.per.day.cv,
                       missing.period = missing.period.cv, beta.int = beta.int.cv,
                       v1 = v1, v0 = v0.candidate[j], max.iter = max.iter,
                       epsilon = epsilon)
        beta.cv <- model.fitting$beta.hat
        pred.error[j, num.fold] <- sum((y.test.unlist - X.test %*% beta.cv)^2)
      }
    }
  }
  
  return(pred.error)
}