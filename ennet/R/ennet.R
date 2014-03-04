#    This is an implementation of ENNET algorithm for Gene Regulatory Network
#    inference from mRNA expression data, in form of an R package.
#    Copyright (C) 2013  Janusz Slawek

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program, see LICENSE.

ennet = function (E      = matrix(rnorm(10000),100,100),
                  K      = matrix(0,nrow(E),ncol(E)),
                  Tf     = 1:ncol(E),
                  s_s    = 1,
                  s_f    = 0.3,
                  M      = 5000,
                  nu     = 0.001,
                  scale  = TRUE,
                  center = TRUE,
                  optimization.stage  = 2) {
  N = ncol(E)
  P = length(Tf)
  S = nrow(E)
  E = scale(E,scale=scale,center=center)
  V = matrix(0, P, N)
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
    result = rep(0, P)
    result[which(!is.na(match(Tf,predictorsI)))]] = model$importance
    result
  }
  
  if (optimization.stage > 0) {
    # first stage of re-evaluation
    s  = apply(V,1,var)
    S  = matrix(rep(s,N),P,N)
    V  = V * S
  }
  
  # prepare the final adjacency matrix
  colnames(V) = colnames(E)
  rownames(V) = colnames(E)[Tf]
  return(V)
}
