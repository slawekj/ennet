ennet = function (E      = matrix(rnorm(10000),100,100),
                  K      = matrix(0,nrow(E),ncol(E)),
                  Tf     = rep(0,ncol(E)),
                  s_s    = 1,
                  s_f    = 0.3,
                  M      = 5000,
                  nu     = 0.001,
                  scale  = TRUE,
                  center = TRUE) {
  N = ncol(E)
  S = nrow(E)
  E = scale(E,scale=scale,center=center)
  V = matrix(0, N, N)
  V = foreach(i = 1:N, .inorder = TRUE, .combine = "cbind") %dopar% {
    predictedI  = i
    predictorsI = setdiff(Tf, predictedI)
    validE      = K[,predictedI] == 0
    
    model = ennet.train(X.train = E[validE,predictorsI],
                        Y.train = E[validE,predictedI],
                        M.train = M,
                        nu = nu,
                        s_s = s_s,
                        s_f = s_f)
    result = rep(0, N)
    result[predictorsI] = model$importance
    result
  }
  
  # first stage of re-evaluation
  s  = apply(V,1,var)
  S1 = matrix(rep(s,N),N,N)
  V  = V * S1
  
  # second stage of re-evaluation
  S2 = matrix(1,N,N)
  for (i in which(colSums(K)>0 && colSums(K)<S)) {
    for (j in 1:N) {
      avg.ko   = mean(K[K[,i]==1,j])
      avg.n.ko = mean(K[K[,i]==0,j])
      std.dev  = sd(K[,j])
      S2[i,j]  = abs(avg.ko-avg.n.ko)/std.dev
    }
  }
  V = V * S2
  
  # prepare the final adjacency matrix
  colnames(V) = colnames(E)
  rownames(V) = colnames(E)
  return(V)
}
