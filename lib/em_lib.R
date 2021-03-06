# return expectation of latent variables 
# input: rankIndex_list(in conv.R)
# output: expectation of latent variables 
E_fun.prob= function(pi_mat, rankIndex_list, mu_e, Sig_e, Omg_e,
                     burn_num, restore_num,
                     parallel = T, verbose = T)
{
  p = ncol(pi_mat)
  n = nrow(pi_mat)
  rank_id = rankIndex_list$rank_id
  rank_index = rankIndex_list$rank_index
  freq_v = rankIndex_list$freq_v
  Ez = rep(0, p)
  Covz = matrix(0,p,p)
  if (parallel == TRUE)
  {
    # i has a value on the range of (1, .., length(rank_id))
    r = foreach(i=1:length(rank_id), .packages = c("MASS", "truncnorm"),
                .export = 'E_fun.prob_inner') %dopar% {
      E_fun.prob_inner(i, pi_mat, rankIndex_list, mu_e, Sig_e, Omg_e,
                       burn_num, restore_num,
                       verbose = F)  
                }
    
    for (i in 1:length(rank_id))
    {
      Ez = Ez + r[[i]]$Ez
      Covz = Covz + r[[i]]$Covz
    }
  }

  if (parallel == FALSE)
  {
    for (i in 1:length(rank_id))
    {
      r = E_fun.prob_inner(i, pi_mat, rankIndex_list, mu_e, Sig_e, Omg_e,
                       burn_num, restore_num,
                       verbose = F) 
      Ez = Ez + r$Ez
      Covz = Covz + r$Covz
    }
  }
  Ez[p] = 0
  Covz[,p] = 0 ; Covz[p,] = 0 ; Covz[p,p] = 1
  return(list(Ez = Ez, Covz = Covz))
}


###

E_fun.prob_inner = function(i, pi_mat, rankIndex_list, mu_e, Sig_e, Omg_e,
                            burn_num, restore_num,
                            verbose = F)
{
  p = ncol(pi_mat)
  n = nrow(pi_mat)
  rank_index = rankIndex_list$rank_index
  freq_v = rankIndex_list$freq_v
  Ez = rep(0, p)
  Covz = matrix(0,p,p)
  idx = which(rank_index == i)
  x = pi_mat[idx[1],]
  # set an initial of z
  init_iter = 1
  while(T)
  {
    if (init_iter > 1000) 
    {
      z = mvrnorm(1, mu_e, diag(5,p))
    } else {
      z = mvrnorm(1, mu_e, Sig_e)    
    }
    zo = order(z, decreasing = T)
    init_iter = init_iter + 1
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
  zm = rep(0,p)
  zc = matrix(0,p,p)
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
      zm = zm + z/restore_num
      zc = zc + z%*%t(z)/restore_num
      j_c = 1
      j_k = j_k + 1
    }
  }
  Ez = zm*freq_v[i]
  Covz = zc*freq_v[i] 
  return(list(Ez = Ez, Covz = Covz))
}

