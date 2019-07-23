rm(list=ls()); gc()
if (Sys.info()[1]=="Linux")
{
  setwd("~/Documents/GitHub/MLE_Thurstone")
} else {
  setwd("~/GitHub/MLE_Thurstone")
}
library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(foreach)
sourceCpp('./lib/test.cpp')

timesTwo(1)
numCores = 4
registerDoParallel(numCores)  # use multicore, set to the number of our cores
a = c(2,1,3)
foreach (i=a, .noexport = "timesTwo") %dopar% {
 timesTwo(i)
}
?foreach
