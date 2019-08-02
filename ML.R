rm(list = ls()); gc()
if (Sys.info()[1]=="Linux")
{
  setwd("~/Documents/GitHub/MLE_Thurstone")
} else {
  setwd("~/GitHub/MLE_Thurstone")
}
source("./lib/pack.R")
library(mlThurstone)
p = 5
mu = seq(2,0, length = p)
Sigma = matrix(0.5, p, p)
Sigma[,p] =  Sigma[p,] = 0
diag(Sigma) = 1

n = 1000
burn_num = 1e+2
restore_num = 5e+3
parallel = T
verbose = F
numCores = 4

if (parallel == T) registerDoParallel(numCores) 
set.seed(1)
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
mu_e = c(fit_model$coefficients*sqrt(2),0) # note: probit reg z~N(0, 1/sqrt(2))
Sig_e = diag(1,p)
Omg_e = solve(Sig_e)
# sampling function
system.time({
  for (iter in 1:1000)
  {
    cat("outer iter:", iter ,'\n')
    E_fit = E_fun.prob_C(pi_mat, rankIndex_list, mu_e, Sig_e, Omg_e,
                       burn_num, restore_num, parallel, verbose)
    mu_e = E_fit$Ez
    Sig_e = E_fit$Covz
    Sig_e = Sig_e - mu_e%*%t(mu_e)
    Sig_e[1,-1] = Sig_e[1,-1]/sqrt(Sig_e[1,1])
    Sig_e[-1,1] = Sig_e[-1,1]/sqrt(Sig_e[1,1])
    Sig_e[1,1] = 1
    ### start:sparse coding
    
    
    ### end:sparse coding
    Omg_e = solve(Sig_e)
    cat("me:", mu_e, '\n')
    cat("Frobenius norm:", sum((Sigma-Sig_e)^2), '\n')
    print(Sig_e)
  }
}
)

