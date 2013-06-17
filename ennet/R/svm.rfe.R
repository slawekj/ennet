svm.rfe = function(Ex,Pt,Tm,Tf,delay=1,kernel="linear",degree=3,cost=1,gamma=1,scale=TRUE,center=TRUE){
  
  # N - number of regulated genes
  N <- ncol(Ex[[1]]);
  
  # scale and combine expression data
  E.all        <- NULL; # all expression data
  P.all        <- NULL; # all perturbated genes
  T.all        <- NULL; # all timestamps
  Tm.rep       <- NULL; # all time series replicates
  Tm.rep.count <- 0;    # time series replicates counter
  
  for (i in 1:length(Ex)) {
    Ex[[i]] <- scale(Ex[[i]],center=center,scale=scale);
    E.all   <- rbind(E.all,Ex[[i]]);
    P.all   <- rbind(P.all,Pt[[i]]);
    T.all   <- rbind(T.all,matrix(Tm[[i]],length(Tm[[i]]),1));
    if(sum(Tm[[i]])>0) {
      # time series replicate
      Tm.rep.count <- Tm.rep.count + 1;
      Tm.rep <- rbind(Tm.rep,matrix(Tm.rep.count,nrow(Ex[[i]]),1));
    } else {
      # non time series replicate
      Tm.rep <- rbind(Tm.rep,matrix(0,nrow(Ex[[i]]),1));
    }
  }
  S <- nrow(E.all); # all the samples
  
  V <- matrix(0,N,N);
  # solve columns of V independently
  V <- foreach(i=1:N, .inorder=TRUE, .combine="cbind") %dopar% {
    predictedI  <- i;
    predictorsI <- setdiff(Tf,predictedI);
    
    T_sample    <- rep(TRUE,S);
    F_sample    <- rep(TRUE,S);
    T_sample[P.all[,predictedI]>0] <- FALSE;
    F_sample[P.all[,predictedI]>0] <- FALSE;
    if (Tm.rep.count > 0) {
      for (j in 1:Tm.rep.count) {
        I <- which(Tm.rep==j);
        if (delay>0 && length(I)>delay) {
          T_sample[I[1:delay]] <- FALSE;
          F_sample[rev(I)[1:delay]] <- FALSE;
        } else {
          T_sample[I] <- FALSE;
          F_sample[I] <- FALSE;
        }
      }
    }
    
    survivingFeaturesIndexes = seq(1:length(predictorsI))
    result                   = rep(0,N)
    ranking                  = rep(0,length(predictorsI))
    
    for (j in 1:(length(predictorsI)-1)) {
      model = svm(
        E.all[F_sample,predictorsI[survivingFeaturesIndexes]],
        E.all[T_sample,predictedI],
        cost=cost,
        gamma=gamma,
        kernel="linear",
        degree=degree,
        scale=FALSE
      )
      w       = t(model$coefs)%*%model$SV
      temp    = w*w
      worst   = which.min(temp)
      survivingFeaturesIndexes = survivingFeaturesIndexes[-worst]
      ranking[survivingFeaturesIndexes] = ranking[survivingFeaturesIndexes] + temp[-worst]
    }
    
    result[predictorsI] <- ranking
    result
  }
  
  colnames(V)<-colnames(Ex[[1]]);
  rownames(V)<-colnames(Ex[[1]]);
  return(V);
}
