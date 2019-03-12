#' \code{coord.descent} the coordinate descent algorithm
#' 
#' @param y: response variable
#' @param X: regressor
#' @param beta.int: initial value for beta
#' @param v0: 1/variance for spike density
#' @param v1: 1/variance for slab density
#' @param pstar: starting value for E(gamma)
#' @param max.iter: maximum iteration for the algorithm
#' @param eps: threshold
#' @return estimated beta

coord.descent <- function(y, X, beta.int, v0, v1, pstar, max.iter = 1000, eps = 0.5) {
  
  n <- length(y)
  p <- dim(X)[2]
  beta <- beta.int
  lambda.star <- v0 * (1-pstar) + v1 * pstar
  r <- c(y - X %*% beta)
  
  for (iter in 1:max.iter) {
    # print(iter)
    
    beta.old <- beta
    
    for (j in 1:p) {
      # if (j %% 50000 == 0) print(j)
      
      # partial residuals
      r <- r + X[,j]*beta[j]
      
      # update beta
      xr <- sum(X[,j]*r)
      xx <- sum(X[,j]^2)   
      beta[j] <- sign(xr) * max(0, (abs(xr)-lambda.star[j])/xx)
      
      r <- r - X[,j]*beta[j]
      
    }
    
    diff.beta.update <- sum((beta - beta.old)^2)
    
    # print(diff.beta.update)
    if (iter > 1) {
      if (abs(diff.beta.update) < eps) {
        break
      } 
    }
    
  }
  return(beta)
}