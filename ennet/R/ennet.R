ennet <- function (Ex, Pt, Tm, Tf, delay = 1, ntree = 5000, shrinkage = 0.001, 
          bag.fraction = 1, col.fraction = 0.3, scale = TRUE, center = TRUE) {
  N <- ncol(Ex[[1]])
  E.all <- NULL
  P.all <- NULL
  T.all <- NULL
  Tm.rep <- NULL
  Tm.rep.count <- 0
  for (i in 1:length(Ex)) {
    Ex[[i]] <- scale(Ex[[i]], center = center, scale = scale)
    E.all <- rbind(E.all, Ex[[i]])
    P.all <- rbind(P.all, Pt[[i]])
    T.all <- rbind(T.all, matrix(Tm[[i]], length(Tm[[i]]), 
                                 1))
    if (sum(Tm[[i]]) > 0) {
      Tm.rep.count <- Tm.rep.count + 1
      Tm.rep <- rbind(Tm.rep, matrix(Tm.rep.count, nrow(Ex[[i]]), 
                                     1))
    }
    else {
      Tm.rep <- rbind(Tm.rep, matrix(0, nrow(Ex[[i]]), 
                                     1))
    }
  }
  S <- nrow(E.all)
  V <- matrix(0, N, N)
  V <- foreach(i = 1:N, .inorder = TRUE, .combine = "cbind") %dopar% {
    predictedI <- i
    predictorsI <- setdiff(Tf, predictedI)
    T_sample <- rep(TRUE, S)
    F_sample <- rep(TRUE, S)
    T_sample[P.all[, predictedI] > 0] <- FALSE
    F_sample[P.all[, predictedI] > 0] <- FALSE
    if (Tm.rep.count > 0) {
      for (j in 1:Tm.rep.count) {
        I <- which(Tm.rep == j)
        if (delay > 0 && length(I) > delay) {
          T_sample[I[1:delay]] <- FALSE
          F_sample[rev(I)[1:delay]] <- FALSE
        }
        else {
          T_sample[I] <- FALSE
          F_sample[I] <- FALSE
        }
      }
    }
    model <- ennet.train(X.train = E.all[F_sample, predictorsI], 
                         Y.train = E.all[T_sample, predictedI], M.train = ntree, 
                         nu = shrinkage, s_s = bag.fraction, s_f = col.fraction)
    result <- rep(0, N)
    result[predictorsI] <- model$importance
    result
  }
  colnames(V) <- colnames(Ex[[1]])
  rownames(V) <- colnames(Ex[[1]])
  return(V)
}
