# rankIndex_fun
# input: obs of rankings, pi_mat[i,j] denotes the rank of the jth item in the ith experiment
# output: rank_id(set of the unique ranks), 
#         rank_index[i]; type of ranks in rank_id in the ith experiment

rankIndex_fun = function(pi_mat)
{ 
  if (class(pi_mat)!='matrix' & class(pi_mat)!='list')
  {
    stop("pi_mat should be matrix or list!")
  }
  
  if (class(pi_mat)=='matrix') 
  {
    n = nrow(pi_mat)
    p = ncol(pi_mat)
  }
  if (class(pi_mat)=='list') n = length(pi_mat)
  
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
  return(list(rank_id = rank_id, rank_index = rank_index))
}