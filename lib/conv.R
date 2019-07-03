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

# genDesignR_fun
# input: number of items to be ranked
# output: design matrix for glmnet or glm
genDesignR_fun <- function(p)
{
  x = matrix(0, p*(p-1), p)
  y = rep(0, p*(p-1) )
  ix = 1
  for (i in 1:p)
  {
    for (j in 1:p)
    {
      if (i == j) next
      jx1 = min(i,j)
      jx2 = max(i,j)
      x[ix,jx1] = 1; x[ix,jx2] = -1
      if (i<j) y[ix] = 1
      ix = ix + 1
    }
  }
  x = x[,-p]
  return(list(x = x, y = y))  
}


# convTpMat_fun
# inpout: obs of rankings
# output: Gmat_hat is win-loss matrix(prob); 
#         Qpmat is weight matrix of design.
convToMat_fun = function(pi_mat)
{
  p = ncol(pi_mat)
  n = nrow(pi_mat)
  Nmat = matrix(0,p,p)
  Wmat = matrix(0,p,p)
  for (i in 1:nrow(pi_mat))
  {
    z = pi_mat[i,]
    for (j in 1:(length(z)-1))
    {
      for (k in (j+1):length(z))
      {
        # ix1: the index of item with rank j
        ix1 = z[j] ; ix2 = z[k]
        Nmat[ix1,ix2] = Nmat[ix2,ix1] = Nmat[ix1,ix2] + 1
        Wmat[ix1,ix2] = Wmat[ix1,ix2] + 1
      }
    }
  }
  Gmat_hat = Wmat/Nmat
  Gmat_hat[Nmat==0] = 0
  Qpmat = Nmat/sum(Nmat)
  return(list(Gmat_hat=Gmat_hat, Qpmat=Qpmat))
}

# convToW_fun
# input: list of Qpmat and Gmat_hat
# output: weight vector used in glm
convToW_fun = function(fit_mat)
{
  Qpmat = fit_mat$Qpmat
  Gmat_hat = fit_mat$Gmat_hat
  p = ncol(Qpmat)
  wmat = Qpmat*Gmat_hat
  wmat = t(wmat)
  wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))] 
  return(wvec)
}


