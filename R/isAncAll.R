isAncAll <- function(graph,as,b){
  
  anc = c() # ancestors of b (not including b)
  
  for (a in as){
    if(isAnc(graph,a,b)){
      anc = c(anc,a)
    }
  }
  
  return(anc)
}

isAnc <- function(graph,a,b,visited=rep(FALSE,nrow(graph)))
{
  
  if (a %in% b){
    return(TRUE)
  }
  
  visited[a] = TRUE;
  
  adj = which(graph[a,] & !visited);
  
  out=FALSE;
  for (j in adj){
    out=isAnc(graph, j, b, visited);
    if(out==TRUE){
      break;
    }
  }
  
  return(out)
  
}