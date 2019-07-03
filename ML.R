rm(list = ls()); gc()
library(MASS)
setwd("~/GitHub/MLE_Thurstone")
if (!require(truncnorm)) {install.packages("truncnorm") ; library(truncnorm)}
source("./lib/conv.R")
source("./lib/em_lib.R")
p = 4
mu = seq(5,0, length = p)
Sigma = diag(1,p)
n = 10000
burn_num = 1e+2
restore_num = 5e+3
verbose = T
#set.seed(1)
z = mvrnorm(n, mu, Sigma)
# pi_mat denotes the index of items
# pi_mat[i,k]: the index of item with rank k in the ith experi
pi_mat =t(apply(z, 1, order, decreasing = T))
rankIndex_list = rankIndex_fun(pi_mat)

# initialization:: pairwise comparisons (package)
fit_mat = convToMat_fun(pi_mat)
d_mat = genDesignR_fun(p)
wvec = convToW_fun(fit_mat)
fit_model = glm.fit(x = d_mat$x, y=d_mat$y, 
        weight = wvec, family = binomial(link='probit'),
        intercept = FALSE)
# Note: ignore warning:In eval(family$initialize)
mu_e = c(fit_model$coefficients*sqrt(2),0) # note: probit reg z~N(0, 1/sqrt(2))
Sig_e = diag(1,p)
Omg_e = solve(Sig_e)
# sampling function

Ez = E_fun.prob(rankIndex_list, mu_e, Omg_e,
           burn_num, restore_num, verbose)
## 
mu_e = colMeans(Ez)
Sig_e = cov(Ez)
Omg_e = solve(Sig_e)

# debugging
pi_hat =t(apply(zmat_mean, 1, order, decreasing = T))
pi_mat
for (i in 1:length(rank_id))
{
  # iteration 
  cat("#iter: ", i, '\n')
  idx = which(rank_index == i)
  x = pi_mat[idx[1],]
  tx = pi_hat[i,]
  if (sum(x == tx)!= p ) stop()
}



