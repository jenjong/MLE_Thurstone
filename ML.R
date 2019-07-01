rm(list = ls()); gc()
library(MASS)
p = 5
mu = seq(5,0, length = p)
Sigma = diag(1,p)
n = 100
set.seed(1)
z = mvrnorm(n, mu, Sigma)
pi_mat =t(apply(z, 1, order, decreasing = T))

# ranking index
rank_index_t = c("", n)
for (i in 1:n)
{
  x = pi_mat[i,]
  if (p<10)
  {
    xt = as.character(x)
    rank_index_t[i] = paste(xt,sep='',collapse = "")
  }
  if (p>=10)
  {
    xt = x
    idx = (nchar(xt)==1)
    xt[idx] = paste0("0",xt[idx])
    rank_index_t[i] = paste(xt,sep='',collapse = "")
  }
}
rank_index_t = factor(rank_index_t)
rank_id = levels(rank_index_t)
rank_index = as.integer(rank_index_t)

# initialization:: pairwise comparisons (package)
mu_e = seq(5,0, length = p)
Sig_e = diag(1,p)
Omg_e = solve(Sig_e)


# sampling
which(rank_index==1)



