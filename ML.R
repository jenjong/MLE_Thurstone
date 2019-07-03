rm(list = ls()); gc()
library(MASS)
setwd("~/GitHub/MLE_Thurstone")
if (!require(truncnorm)) {install.packages("truncnorm") ; library(truncnorm)}
source("./lib/conv.R")
p = 5
mu = seq(5,0, length = p)
Sigma = diag(1,p)
n = 200
verbose = T
#set.seed(1)
z = mvrnorm(n, mu, Sigma)
# pi_mat denotes the ranks of items
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
burn_num = 1e+4
restore_num = 5e+3
p = ncol(pi_mat)
rank_id = rankIndex_list$rank_id
rank_index = rankIndex_list$rank_index
zmat_mean = matrix(0,length(rank_id),p)
for (i in 1:length(rank_id))
{
  # 
  cat("#iter: ", i, '\n')
  idx = which(rank_index == i)
  x = pi_mat[idx[1],]
  # set an initial of z
  while(T)
  {
    z = mvrnorm(1, mu_e, Sig_e)  
    zo = order(z, decreasing = T)
    if (!any(x != zo)) 
    {
      if (verbose) cat("Initialization is complete! ")
      break
    }
  }
  # burning
  # ix[j_c] : the index of item with rank j_c
  # iterates according to ranks (j_c from 1 to p)
  j_c = 1
  ix = sort(x, index.return = T)$ix
  for (j in rep(ix, burn_num))
  {
    m = mu_e[j]- sum(Omg_e[j,][-j]*(z[-j] - mu_e[-j]))/Omg_e[j,j]
    v = sqrt(1/Omg_e[j,j])
    if (j_c == 1) b = Inf else b = z[ix[j_c-1]]
    if (j_c == p) a = -Inf else a = z[ix[j_c+1]]
    if (a>=b) stop("error 1")
    z[j] = rtruncnorm(1, a=a, b=b, mean = m, sd = v)
    if (is.na(z[j])) stop("NA occurs")
    j_c = j_c + 1
    if (j_c>p) j_c = 1
  }
  
  # restore
  j_k = 1
  j_c = 1
  zmat = matrix(0, restore_num, p)
  for (j in rep(ix, restore_num))
  {
    m = mu_e[j]- sum(Omg_e[j,][-j]*(z[-j] - mu_e[-j]))/Omg_e[j,j]
    v = sqrt(1/Omg_e[j,j])
    if (j_c == 1) b = Inf else b = z[ix[j_c-1]]
    if (j_c == p) a = -Inf else a = z[ix[j_c+1]]
    tz = rtruncnorm(1, a=a, b=b, mean = m, sd = v)
    z[j] = tz
    j_c = j_c + 1
    if (j_c>p)
    {
      zmat[j_k,] = z
      j_c = 1
      j_k = j_k + 1
    }
  }
  
  zmat_mean[i,] = colMeans(zmat)
}





boxplot(zmat)

#while(T)
{}



