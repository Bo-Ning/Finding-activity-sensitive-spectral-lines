rm(list = ls())
setwd("PATH")
##################### EMVS ########################
# load activtiy indices
S_index <- read.table("data/activity-indices/S-index.txt")$V1
BIS <- read.table("data/activity-indices/BIS.txt")$V1
FWHM <- read.table("data/activity-indices/FWHM.txt")$V1
NaD <- read.table("data/activity-indices/NaD.txt")$V1
Halpha <- read.table("data/activity-indices/Halpha.txt")$V1

# load outliers
outlier_index <- read.csv("data/outlier_index.csv")$x
# choose a activity indices
activity.index <- scale(BIS[-outlier_index])

# load wavelength
Highest.order <- 22
pb  <- txtProgressBar(1, Highest.order-20, style=3)    # report progress
cat("\nStart loading wavelength: \n")     # report progress
wvlth <- NULL
p.wvlth.vec <- rep(0, Highest.order-20)
for (order in 21:Highest.order) {
  # report progress
  setTxtProgressBar(pb, order-20)
  
  wvlth_per_order <- read.table(paste("data/wavelengths/wavelength",
                                      order,".txt", sep = ""))$x
  wvlth_per_order <- c(wvlth_per_order)
  names(wvlth_per_order) <- NULL
  
  wvlth <- c(wvlth, wvlth_per_order)
  p.wvlth.vec[order-20] <- length(wvlth_per_order)
}

# load data
pb  <- txtProgressBar(1, Highest.order-20, style=3)    # report progress
cat("\nStart loading pixels: \n")     # report progress
data <- NULL
p.vec <- rep(0, Highest.order-20)
for (order in 21:Highest.order) {
  # report progress
  setTxtProgressBar(pb, order-20)
  
  data_per_order <- read.table(paste("data/normalized-spectrum/normRVInterp",
                                     order,".txt", sep = ""))
  data_per_order <- as.matrix(data_per_order)
  rownames(data_per_order) <- colnames(data_per_order) <- NULL
  
  data <- cbind(data, data_per_order)
  p.vec[order-20] <- dim(data_per_order)[2]
}
data <- data[-outlier_index, ] # delete outliers

# load Julian.date
time <- read.table(file = "data/JDs.txt")$V1
time <- time[-outlier_index]
counts.per.day <- NULL
# insert missing time points
for (t in 1:82) {
  counts.per.day <- 
    c(counts.per.day, length(which(time >= 2455278+t-1 & time < 2455278+t)))
}

missing.period <- which(counts.per.day == 0)
for (t in 1:82) {
  if (t == 1) {
    start <- 1
    end <- counts.per.day[1:t]
    y <- list(activity.index[start:end])
    dataset.index <- list(start:end)
  } else {
    if (counts.per.day[t] == 0) {
      y <- c(y, NA)
      dataset.index <- c(dataset.index, NA)
    } else {
      start <- sum(counts.per.day[1:(t-1)]) + 1 
      end <- sum(counts.per.day[1:t])
      y <- c(y, list(activity.index[start:end])) 
      dataset.index <- c(dataset.index, list(start:end))
    }
  }
}

y <- y
X <- t(scale(t(data))) # scale and center data
T <- 82

source("code/main.R")
lambda.candidate <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
v0.candidate <- 10^c(1:10)
res <- model.fitting(T, activity = y, pixels = X, lambda.candidate, v0.candidate, 
                     counts.per.day, time, missing.period, 
                     max.iter.step1 = 30, max.iter.step2 = 20, max.iter.step3 = 100)
# collect results
beta.hat <- res$beta.hat
pstar.hat <- res$pstar.hat
logpost <- res$logpost
sigma.sq_t.res <- res$sigma.sq_t.res
residual <- res$residual
logpost.lambda <- res$logpost.lambda
pred.error <- res$pred.error.sum

# save results
write.csv(sigma.sq_t.res, file = "sigma.sq.csv")
write.csv(beta.hat, file = "beta.hat.csv")
write.csv(pstar.hat, file = "pstar.csv")
write.csv(logpost, file = "logpost.csv")
write.csv(residual, file = "residual.csv")
write.csv(logpost.lambda, file = "logpost.lambda")
write.csv(pred.error, file = "pred.error")
