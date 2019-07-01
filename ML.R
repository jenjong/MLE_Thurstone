rm(list = ls()); gc()
library(MASS)
library(truncnorm)
p = 5
mu = seq(5,0, length = p)
Sigma = diag(1,p)
n = 100
set.seed(1)
z = mvrnorm(n, mu, Sigma)
# pi_mat denotes the ranks of items
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
i = 9 ; length(rank_id)
idx = which(rank_index == i)
x = pi_mat[idx[1],]
# set an initial of z
while(T)
{
  z = mvrnorm(1, mu_e, Sig_e)  
  zo = order(z, decreasing = T)
  if (!any(x != zo)) break
}


burn_num = 1e+4
restore_num = 5e+3
# burning
ix = sort(x, index.return = T)$ix
j_c = 1
# ix[j_c] : the index of item with rank j_c
# iterates according to ranks
for (j in rep(ix, burn_num))
{
  m = mu_e[j]- sum(Omg_e[j,][-j]*(z[-j] - mu_e[-j]))/Omg_e[j,j]
  v = sqrt(1/Omg_e[j,j])
  if (j_c == 1) b = Inf else b = z[ix[j_c-1]]
  if (j == p) a = -Inf else a = z[ix[j_c+1]]
  z[j] = rtruncnorm(1, a=a, b=b, mean = m, sd = v)
  j_c = j_c + 1
  if (j_c>p) j_c = 1
}
# restore
z_tmp = rep(0, p*restore_num)
j_k = 1
for (j in rep(1:p, restore_num))
{
  m = mu_e[j]- sum(Omg_e[j,][-j]*(z[-j] - mu_e[-j]))/Omg_e[j,j]
  v = sqrt(1/Omg_e[j,j])
  if (j == 1) b = Inf else b = z[j-1]
  if (j == p) a = -Inf else a = z[j+1]
  tz = rtruncnorm(1, a=a, b=b, mean = m, sd = v)
  z[j] = tz
  z_tmp[j_k]  = tz
  j_k = j_k + 1
}
zmat = matrix(z_tmp, restore_num, p, byrow = T)
colMeans(zmat)
boxplot(zmat)


#while(T)
{}



