thetafit <- function(x, parallel = TRUE, ncores = getOption("mc.cores", 2L), intercept = FALSE, K = 20, l = 5, seed = NULL, verbose = FALSE){
  require("parallel")
  require("doParallel")
  require("foreach")
  require("doSNOW")
  require("doRNG")
  
  
  if (parallel) {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    if (verbose)
      cat("running nodewise regressions using cv lambda choice.")
    
    output <- foreach(i = seq(p), .packages = c("midasml")) %dopar% {
      fit <- tscv.sglfit(x[,-i], x[,i], gamma = 1.0, K = K, l = l, seed = seed, intercept = intercept, standardize = FALSE)
      coeffs <- as.vector(fit$cv.fit$lam.min$beta) 
      # get lambda.min
      lambda <- fit$lambda.min
      # C[-i,i] elements of ThetaHat
      c_ii <- -as.vector(coeffs)
      # tau^2_i = ||X_i - X_{-i}gammahat_i||2_T + lambda_j | gammahat_i|_1 
      t2 <- as.numeric((x[,i] %*%(x[,i] - predict(fit,x[,-i], s = lambda)))/n)
      out <- list(c_ii = c_ii, t2 = t2, lambda = lambda)
      out
    }
    stopImplicitCluster()
    stopCluster(cl)
  } else {
    stop("only parallel computation is implemented")
  } 
  sort.output <- function(output,p){
    # Purpose: sort the output that comes from parallel nodewise LASSO regression computations.
    k <- length(output)
    if (k!=p) {
      stop("output from parallel nodewise regressions != p. check.")
    }
    C <- diag(rep(1,p))
    T2 <- lambda.seq <- numeric(p)
    
    for (i in 1:p){
      tmp <- output[[i]]
      C[-i,i] <- tmp$c_ii
      T2[i] <- tmp$t2
      if (!is.null(tmp$lambda)){
        lambda.seq[i] <- tmp$lambda
      }
    }
    return(list(C=C,T2=T2,lambda.seq=lambda.seq))
  }
  output.sorted <- sort.output(output,p)
  C <- output.sorted$C
  T2 <- output.sorted$T2
  lambda.seq <- output.sorted$lambda.seq
  thetahat <- C %*% solve(diag(T2))
  ##this is thetahat ^ T!!
  thetahat <- t(thetahat)
  
  if(all(thetahat[lower.tri(thetahat)] == 0,
         thetahat[upper.tri(thetahat)] == 0) && verbose)
    cat("Thetahat is a diagonal matrix!\n")
  return(thetahat)
}