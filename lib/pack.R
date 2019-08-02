# package
library(MASS)
if (!require(truncnorm)) {install.packages("truncnorm") ; library(truncnorm)}
if (!require(doParallel)) {install.packages("doParallel") ; library(doParallel)}
if (!require(foreach)) {install.packages("foreach") ; library(foreach)}
source("./lib/conv.R")
source("./lib/em_lib.R")
source("./lib/em_lib_C.R")  
