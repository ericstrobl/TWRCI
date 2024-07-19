my_skeleton <- function (suffStat, p, SNPs=NULL, alpha=0.01, Gest = NULL, order=NULL, indepTest=earth_wrap, m.max=Inf) 
{
  
  if (!is.null(Gest)){
    G = (t(Gest) + Gest) > 0
  } else{
    G = matrix(TRUE,p,p)
  }
  diag(G) = FALSE 

  seq_p <- seq_len(p)
  sepset <- lapply(seq_p, function(.) vector("list", p))
  done <- FALSE
  ord <- 0L
  n.edgetests <- numeric(1)
  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord + 1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE) ####
    ind <- ind[order(ind[, 1]), ,drop=FALSE]
    remEdges <- nrow(ind)
    
    G.l <- split(G, gl(p, p)) # potential parents
    # print(ncol(G))
    for (i in 1:remEdges) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      if (G[x,y]) {
        nbrsBool <- G.l[[x]] #### potential parents of x
        nbrsBool[y] <- FALSE
        if (!is.null(Gest)){
          nbrs <- intersect(seq_p[nbrsBool],which(Gest[,x])) #### potential parents of x
        } else{
          nbrs <- seq_p[nbrsBool]
        }

        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord)
            done <- FALSE
          S <- seq_len(ord)
          repeat {
            n.edgetests[ord1] <- n.edgetests[ord1] + 1
            if (!is.null(SNPs)){
              pval <- indepTest(x, y, nbrs[S], SNPs[[x]], suffStat)
            } else{
              pval <- indepTest(x, y, nbrs[S], NULL, suffStat)
            }
            if (is.na(pval))
              pval <- as.numeric(NAdelete)
            if (pval >= alpha) {
              G[x, y] <- G[y, x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if (nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          }
        }
      }
    }
    ord <- ord + 1L
  }


  
  return(list(G=G,sepset=sepset))
}

# pcor_test <- function(x,y,z,suffStat){
#   
#   require(Rfast)
#   
#   n <- suffStat$n
#   gp <- length(z)
#   
#   is = c(x,y,z)
#   if (length(z)>0){
#     icvx <- spdinv(suffStat$C[is,is])
#     pcor <- -cov2cor(icvx)[1, 2]
#   } else{
#     pcor <- cov2cor(suffStat$C[is,is])[1, 2]
#   }
#   
#   statistic <- pcor * sqrt((n - 2 - gp)/(1 - pcor^2)) # assumes normality, but rapid convergence with non-normal samples
#   p.value <- 2 * pt(-abs(statistic), (n - 2 - gp))
#   
#   return( p.value )
# }
