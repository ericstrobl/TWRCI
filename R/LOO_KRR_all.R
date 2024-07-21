LOO_KRR_all <- function(tr_y, tr_e, dH){
  err = tr_e/(1-dH)
  pre = tr_y - err
  # 
  # pre = pre/rowSums(pre)
  # d = ncol(pre);
  # pre[is.nan(pre)] = 1/d
  # 
  # err = tr_y - pre
  
  list( pre=pre, err=mean(err^2), resid=err )
}