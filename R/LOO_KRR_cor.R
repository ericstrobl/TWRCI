LOO_KRR_cor <- function(tr_y, tr_e, dH){
  err = tr_e/(1-dH)
  pre = tr_y - err
  # 
  # pre = pre/rowSums(pre)
  # d = ncol(pre);
  # pre[is.nan(pre)] = 1/d
  # 
  # err = tr_y - pre
  err = c()
  for (i in 1:ncol(tr_y)){
    # err = c(err,-abs(cor(tr_y[,i],pre[,i])))
    err = c(err,-cor(tr_y[,i],pre[,i])) # maximize positive correlation
  }
  # print(err)
  
  list( pre=pre, err=err, resid=err )
}