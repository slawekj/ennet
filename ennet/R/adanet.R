adanet <- function(Ex,Pt,Tm,Tf,delay=1,
                   c=30,t=10,psi=0.25,xi=0.75,delta=0.05,
                   interaction.depth=1,n.minobsinnode=5,shrinkage=0.01,bag.fraction=0.67,
                   scale=TRUE,center=TRUE) {
  
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
  #E.all <- data.frame(E.all);
  
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
    
    sortI  <- sort(E.all[F_sample,predictedI],index.return=TRUE)[[2]];
    Nsamp  <- sum(T_sample);
    Y      <- rep(0,Nsamp);
    result <- rep(0,N);
    
    for (m in seq(psi,xi,by=(xi-psi)/c)) {
      L    <- round((m-delta)*Nsamp);
      U    <- round((m+delta)*Nsamp);
      negI <- sortI[1:L];
      posI <- sortI[U:Nsamp];
      Y[posI] <-  1;
      Y[negI] <-  0;
      
      model       <- gbm.fit(
        E.all[F_sample,predictorsI],
        Y,
        distribution="adaboost",
        n.trees=t,
        interaction.depth=interaction.depth,
        n.minobsinnode=n.minobsinnode,
        shrinkage=shrinkage,
        bag.fraction=bag.fraction,
        verbose=FALSE
      );
      
      result[predictorsI] <- result[predictorsI] + summary.gbm(model,n.trees=t,plotit=FALSE,order=FALSE)[,2];
    }
    result <- result / c;
  }
  
  colnames(V)<-colnames(Ex[[1]]);
  rownames(V)<-colnames(Ex[[1]]);
  return(V);
}
