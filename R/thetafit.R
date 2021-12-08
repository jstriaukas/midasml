#' Nodewise LASSO regressions to fit the precision matrix \ifelse{html}{\out{&Theta;}}{\eqn{\Theta}}.
#' 
#' @description 
#' Fits the precision matrix \ifelse{html}{\out{&Theta;}}{\eqn{\Theta}} by running nodewise LASSO regressions. 
#' 
#' @details
#' The function runs \link{tscv.sglfit} \code{p} times by regressing \code{j}-th covariate on all other covariates excluding \code{j}-th covariate. The precision matrix is then constructed based on LASSO estimates. Each nodewise LASSO regression tuning parameter \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}} is optimized using time series cross-validation. See \link{tscv.sglfit} for more details on cross-validation implementation.
#' @usage 
#' thetafit(x, parallel = FALSE, ncores = getOption("mc.cores", NULL), 
#'          intercept = FALSE, K = 20, l = 5, seed = NULL, verbose = FALSE, 
#'          registerpar = TRUE, ...)
#' @param x T by p data matrix, where T and p respectively denote the sample size and the number of regressors.
#' @param parallel if \code{TRUE}, use parallel foreach to fit nodewise LASSO regressions. Parallel registered within the function.
#' @param ncores number of cores used in parallelization
#' @param intercept whether intercept be fitted (\code{TRUE}) or set to zero (\code{FALSE}). Default is \code{FALSE}.
#' @param K number of folds of the cv loop. Default set to \code{20}.
#' @param l the gap used to drop observations round test set data. See \link{tscv.sglfit} for more details.
#' @param seed set a value for seed to control results replication, i.e. \code{set.seed(seed)} is used. \code{seed} is stored in the output list. Default set to \code{as.numeric(Sys.Date())}.
#' @param verbose if \code{TRUE}, prints progress bar. Default set to \code{FALSE}.
#' @param registerpar if \code{TRUE}, register parallelization using \code{registerDoParallel}. Default set to \code{TRUE}.
#' @param ... Other arguments that can be passed to \link{tscv.sglfit}.
#' @return thetafit object.
#' @author Jonas Striaukas
#' @examples
#' \donttest{ 
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' thetafit(x = x, parallel = FALSE)
#' }
#' @export thetafit
thetafit <- function(x, parallel = FALSE, ncores = getOption("mc.cores", NULL), intercept = FALSE, K = 20, l = 5, seed = NULL, verbose = FALSE, registerpar = TRUE, ...){
  np <- dim(x)
  n <- np[1]
  p <- np[2]
  if (is.null(seed)){
    seed <- as.numeric(Sys.Date())
  }
  if (parallel) {
    if (is.null(ncores)){
      ncores <- parallel::detectCores(all.tests = TRUE, logical = TRUE)
    }
    cl <- NULL
    if (registerpar){
      cl <- parallel::makeCluster(ncores)
      registerDoParallel(cl)
    }
    opts <- NULL
    if (verbose){
      pb <- utils::txtProgressBar(max=p, style=3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
    }
    i <- 1
    output <- foreach(i = 1:p, .packages = c("midasml"), .options.RNG=seed, .options.snow=opts) %dorng% {
      fit <- tscv.sglfit(x[,-i], x[,i], gamma = 1.0, K = K, l = l, seed = NULL, intercept = FALSE, standardize = FALSE, ...)
      coeffs <- as.vector(fit$cv.fit$lam.min$beta) 
      # get lambda.min
      lambda <- fit$lambda.min
      # C[-i,i] elements of ThetaHat
      c_ii <- -as.vector(coeffs)
      # tau^2_i = ||X_i - X_{-i}gammahat_i||2_T + lambda_j | gammahat_i|_1 
      t2 <- as.numeric((x[,i] %*%(x[,i] - predict.cv.sglfit(fit, x[,-i], s = lambda)))/n)
      out <- list(c_ii = c_ii, t2 = t2, lambda = lambda)
      out
    }
    if (!is.null(cl)){
      parallel::stopCluster(cl)
    }
  } else {
    i <- 1
    output <- foreach(i = 1:p, .packages = c("midasml")) %do% {
      fit <- tscv.sglfit(x[,-i], x[,i], gamma = 1.0, K = K, l = l, seed = seed, intercept = intercept, standardize = FALSE)
      coeffs <- as.vector(fit$cv.fit$lam.min$beta) 
      # get lambda.min
      lambda <- fit$lambda.min
      # C[-i,i] elements of ThetaHat
      c_ii <- -as.vector(coeffs)
      # tau^2_i = ||X_i - X_{-i}gammahat_i||2_T + lambda_j | gammahat_i|_1 
      t2 <- as.numeric((x[,i] %*%(x[,i] - predict.cv.sglfit(fit, x[,-i], s = lambda)))/n)
      out <- list(c_ii = c_ii, t2 = t2, lambda = lambda)
      out
    }
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
  return(list(thetahat = thetahat, seed = seed))
}