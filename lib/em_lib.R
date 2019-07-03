# return expectation of latent variables 
# input: rankIndex_list(in conv.R)
# output: expectation of latent variables 
E_fun.prob= function(rankIndex_list, mu_e, Omega_e,
                     burn_num, restore_num,
                     verbose = T)
{
  p = ncol(pi_mat)
  n = nrow(pi_mat)
  rank_id = rankIndex_list$rank_id
  rank_index = rankIndex_list$rank_index
  zmat_mean = matrix(0,length(rank_id),p)
  for (i in 1:length(rank_id))
  {
    # 
    if (verbose) cat("#iter:", i, "  ")
    idx = which(rank_index == i)
    x = pi_mat[idx[1],]
    # set an initial of z
    while(T)
    {
      z = mvrnorm(1, mu_e, diag(1,p))  
      zo = order(z, decreasing = T)
      if (!any(x != zo)) 
      {
        if (verbose) cat("Initialization is complete!\n")
        break
      }
    }
    # burning
    # ix[j_c] : the index of item with rank j_c
    # iterates according to ranks (j_c from 1 to p)
    j_c = 1
    for (j in rep(x, burn_num))
    {
      m = mu_e[j]- sum(Omg_e[j,][-j]*(z[-j] - mu_e[-j]))/Omg_e[j,j]
      v = sqrt(1/Omg_e[j,j])
      if (j_c == 1) b = Inf else b = z[x[j_c-1]]
      if (j_c == p) a = -Inf else a = z[x[j_c+1]]
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
    for (j in rep(x, restore_num))
    {
      m = mu_e[j]- sum(Omg_e[j,][-j]*(z[-j] - mu_e[-j]))/Omg_e[j,j]
      v = sqrt(1/Omg_e[j,j])
      if (j_c == 1) b = Inf else b = z[x[j_c-1]]
      if (j_c == p) a = -Inf else a = z[x[j_c+1]]
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
  
  Ez = matrix(0, n, p)
  
  for (i in 1:length(rank_id))
  {
    idx = which(rank_index == i)
    for (j in idx)
    {
      Ez[j,] = zmat_mean[i,]  
    }
  }
  return(Ez)
}

