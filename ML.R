rm(list = ls()); gc()
library(MASS)
if (Sys.info()[1]=="Linux")
{
  setwd("~/Documents/GitHub/MLE_Thurstone")
} else {
  setwd("~/GitHub/MLE_Thurstone")
}

if (!require(truncnorm)) {install.packages("truncnorm") ; library(truncnorm)}
source("./lib/conv.R")
source("./lib/em_lib.R")
p = 5
mu = seq(2,0, length = p)
Sigma = matrix(0.5, p, p)
Sigma[,p] =  Sigma[p,] = 0
diag(Sigma) = 1

n = 10000
burn_num = 1e+2
restore_num = 1e+3
verbose = F
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
for (iter in 1:100)
{
  cat("outer iter:", iter ,'\n')
  E_fit = E_fun.prob(rankIndex_list, mu_e, Sig_e, Omg_e,
                     burn_num, restore_num, verbose)
  ##
  mu_e = colMeans(E_fit$Ez)
  for (j in 1:n)
  {
    if (j == 1) Sig_e = E_fit$Covz[[j]] else Sig_e = Sig_e + E_fit$Covz[[j]]
  }
  Sig_e = Sig_e/n
  Sig_e = Sig_e - mu_e%*%t(mu_e)
  Sig_e[1,-1] = Sig_e[1,-1]/sqrt(Sig_e[1,1])
  Sig_e[-1,1] = Sig_e[-1,1]/sqrt(Sig_e[1,1])
  Sig_e[1,1] = 1
  Omg_e = solve(Sig_e)
  cat("me:", mu_e, '\n')
  cat("Frobenius norm:", sum((Sigma-Sig_e)^2), '\n')
  print(Sig_e)
}
