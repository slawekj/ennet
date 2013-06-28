v2l = function(V,max=100000) {
  V            = data.frame(V)
  N            = nrow(V)
  M            = ncol(V)
  targets      = colnames(V)
  tfs          = rownames(V)
  connList     = data.frame(matrix(0,N*M,3))
  
  connList[,1] = rep(tfs,M)
  connList[,2] = rep(targets,1,each=N)
  connList[,3] = as.vector(unlist(V))
  
  self = seq(from=1,to=(N*M),by=(N+1))
  
  connList     = connList[-self,]
  connList     = arrange(connList,X3,decreasing=TRUE)
  
  connList[,3] = connList[,3]/connList[1,3]
  if (nrow(connList)>max) {
    connList     = connList[1:max,]
  }
  return(connList)
}