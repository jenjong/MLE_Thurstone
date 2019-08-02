# return expectation of latent variables 
# input: rankIndex_list(in conv.R)
# output: expectation of latent variables 
E_fun.prob_C= function(pi_mat, rankIndex_list, mu_e, Sig_e, Omg_e,
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
    r = foreach(i=1:length(rank_id), .packages = c("MASS", "truncnorm", "test2"),
                .export = 'E_fun_C') %dopar% {
      E_fun_C(i, pi_mat, rankIndex_list, mu_e, Sig_e, Omg_e,
                       burn_num, restore_num)
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
      r = E_fun_C(i, pi_mat, rankIndex_list, mu_e, Sig_e, Omg_e,
                       burn_num, restore_num)
      Ez = Ez + r$Ez
      Covz = Covz + r$Covz
    }
  }
  Ez[p] = 0
  Covz[,p] = 0 ; Covz[p,] = 0 ; Covz[p,p] = 1
  return(list(Ez = Ez, Covz = Covz))
}


