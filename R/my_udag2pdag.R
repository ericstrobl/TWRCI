my_udag2pdag <- function (gInput, sepset) 
{
  g <- gInput
  if (sum(g)>0){
    p <- as.numeric(dim(g)[1])
    pdag <- g
    ind <- which(g == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(g[y, ] == 1), x)
      for (z in allZ) {
        if (g[x, z] == 0 && !(y %in% sepset[[x]][[z]] || 
                              y %in% sepset[[z]][[x]])) {

          pdag[x, y] <- pdag[z, y] <- 1
          pdag[y, x] <- pdag[y, z] <- 0
        }
      }
    }
    res2 <- pdag2dag(as(pdag, "graphNEL"))
    if (res2$success) {
      old_pdag <- matrix(0, p, p)
      while (!all(old_pdag == pdag)) {
        old_pdag <- pdag
        ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[b, ] == 1 & pdag[, b] == 
                           1) & (pdag[a, ] == 0 & pdag[, a] == 0))
          if (length(indC) > 0) {
            pdag[b, indC] <- 1
            pdag[indC, b] <- 0
           
          }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[a, ] == 1 & pdag[, a] == 
                           0) & (pdag[, b] == 1 & pdag[b, ] == 0))
          if (length(indC) > 0) {
            pdag[a, b] <- 1
            pdag[b, a] <- 0
            
          }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[a, ] == 1 & pdag[, a] == 
                           1) & (pdag[, b] == 1 & pdag[b, ] == 0))
          if (length(indC) >= 2) {
            g2 <- pdag[indC, indC]
            if (length(g2) <= 1) {
              g2 <- 0
            }
            else {
              diag(g2) <- rep(1, length(indC))
            }
            if (any(g2 == 0)) {
              pdag[a, b] <- 1
              pdag[b, a] <- 0
              
            }
          }
        }
      }
     
    }
  }
  return(g)
}