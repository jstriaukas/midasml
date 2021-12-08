#' Time series cross-validation fit for sg-LASSO
#' 
#' @description 
#' Does k-fold time series cross-validation for sg-LASSO regression model.
#' 
#' The function runs \link{sglfit} \code{nfolds+1} times; the first to get the path solution in \ifelse{html}{\out{&lambda;}}{\code{lambda}} sequence, the rest to compute the fit with each of the folds omitted. 
#' The average error and standard deviation over the folds is computed, and the optimal regression coefficients are returned for \code{lam.min} and \code{lam.1se}. Solutions are computed for a fixed \ifelse{html}{\out{&gamma;}}{\eqn{\gamma}.}
#'
#' @details
#' \ifelse{html}{\out{The cross-validation is run for sg-LASSO linear model. The sequence of linear regression models implied by &lambda; vector is fit by block coordinate-descent. The objective function is  <br><br> <center> ||y - &iota;&alpha; - x&beta;||<sup>2</sup><sub>T</sub> + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where &iota;&#8712;R<sup>T</sup>enter> and ||u||<sup>2</sup><sub>T</sub>=&#60;u,u&#62;/T is the empirical inner product. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on  &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)|&beta;|<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{The cross-validation is run for sg-LASSO linear model. The sequence of linear regression models implied by \eqn{\lambda} vector is fit by block coordinate-descent. The objective function is \deqn{\|y-\iota\alpha - x\beta\|^2_{T} + 2\lambda \Omega_\gamma(\beta),} where \eqn{\iota\in R^T} and \eqn{\|u\|^2_T = \langle u,u \rangle / T} is the empirical inner product. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)|\beta|_{2,1},} a convex combination of LASSO and group LASSO penalty functions.}     
#' @usage 
#' tscv.sglfit(x, y, lambda = NULL, gamma = 1.0, gindex = 1:p, 
#'   K = 20, l = 5, parallel = FALSE, seed = NULL, ...)
#' @param x T by p data matrix, where T and p respectively denote the sample size and the number of regressors.
#' @param y T by 1 response variable.
#' @param lambda a user-supplied lambda sequence. By leaving this option unspecified (recommended), users can have the program compute its own \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}} sequence based on \ifelse{html}{\out{<code>nlambda</code>}}{\code{nlambda}} and \ifelse{html}{\out{&gamma;}}{\eqn{\gamma}} \code{lambda.factor.} It is better to supply, if necessary, a decreasing sequence of lambda values than a single (small) value, as warm-starts are used in the optimization algorithm. The program will ensure that the user-supplied \code{lambda} sequence is sorted in decreasing order before fitting the model.
#' @param gamma sg-LASSO mixing parameter. \ifelse{html}{\out{&gamma;}}{\eqn{\gamma}} = 1 gives LASSO solution and \ifelse{html}{\out{&gamma;}}{\eqn{\gamma}} = 0 gives group LASSO solution.
#' @param gindex p by 1 vector indicating group membership of each covariate.
#' @param K number of folds of the cv loop. Default set to \code{20}.
#' @param l the gap used to drop observations round test set data. For each test observation (in total K), we drop 2l observations such that the test observation at time t is separated by \code{l} observations. Default set to \code{5}.
#' @param parallel if \code{TRUE}, use parallel foreach to fit each fold. Must register parallel before hand, such as doMC or others. See the example below.
#' @param seed set a value for seed to control results replication, i.e. \code{set.seed(seed)} is used. \code{seed} is stored in the output list. Default set to \code{as.numeric(Sys.Date())}.
#' @param ... Other arguments that can be passed to \link{sglfit}.
#' @return tscv.sglfit object.
#' @author Jonas Striaukas
#' @examples
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' beta = c(5,4,3,2,1,rep(0, times = 15))
#' y = x%*%beta + rnorm(100)
#' gindex = sort(rep(1:4,times=5))
#' tscv.sglfit(x = x, y = y, gindex = gindex, gamma = 0.5, 
#'   standardize = FALSE, intercept = FALSE)
#' \dontrun{ 
#' # Parallel
#' require(doMC)
#' registerDoMC(cores = 2)
#' x = matrix(rnorm(1000 * 20), 1000, 20)
#' beta = c(5,4,3,2,1,rep(0, times = 15))
#' y = x%*%beta + rnorm(1000)
#' gindex = sort(rep(1:4,times=5))
#' system.time(tscv.sglfit(x = x, y = y, gindex = gindex, gamma = 0.5, 
#'   standardize = FALSE, intercept = FALSE))
#' system.time(tscv.sglfit(x = x, y = y, gindex = gindex, gamma = 0.5, 
#'   standardize = FALSE, intercept = FALSE, parallel = TRUE))
#' }
#' @export tscv.sglfit
tscv.sglfit <- function(x, y, lambda = NULL, gamma = 1.0, gindex = 1:p, K = 20, l = 5, parallel = FALSE, seed = NULL, ...) {
  N <- nrow(x)
  p <- ncol(x)
  y <- drop(y)
  sglfit.object <- sglfit(x, y, lambda = lambda, gamma = gamma, gindex = gindex, method = "single",...)
  lambda <- sglfit.object$lambda
  # predict -> coef
  nz <- sapply(coef(sglfit.object, type = "nonzero"), length)
  if (l<=1) {
    stop("l must be at least 2; l=5 recommended")
  } 
  if (K < 1) 
    stop("K must be at least 1; K=5 recommended")
  if (K >= N) 
    stop("K must be at most T-1; K=5 recommended")
  if (is.null(seed)){
    seed <- as.numeric(Sys.Date())
  }
  set.seed(seed)
  nk <- floor(N/K)
  foldid <- matrix(sample(1:N)[1:(nk*K)], nrow = nk)
  outlist <- vector("list", length = K)
  if (parallel){
    outlist <- foreach(i = seq(K), .packages = c("midasml")) %dopar%{
      whichfoldnot <- foldid[,i]
      whichgaptrain <- computegapobs(whichfoldnot, N, l)
      y_sub <- y[whichgaptrain]
      outlist[[i]] <- sglfit(x = x[whichgaptrain, , drop = FALSE], 
                             y = y_sub, lambda = lambda, gamma = gamma, gindex = gindex, method = "single", ...)
    }
  } else {
    for (i in seq(K)) {
      whichfoldnot <- foldid[,i]
      whichgaptrain <- computegapobs(whichfoldnot, N, l)
      y_sub <- y[whichgaptrain]
      outlist[[i]] <- sglfit(x = x[whichgaptrain, , drop = FALSE], 
                             y = y_sub, lambda = lambda, gamma = gamma, gindex = gindex, method = "single", ...)
    }
  }
  cvstuff <- tscv.sglpath(outlist, lambda = lambda, x = x, y = y, foldid = foldid, ...)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  lamin <- getmin(lambda, cvm, cvsd)
  idxmin <- which(lamin$lambda.min == lambda)
  idx1se <- which(lamin$lambda.1se == lambda)
  cv.fit <- list(lam.min = list(b0 = sglfit.object$b0[idxmin], beta = sglfit.object$beta[,idxmin]), 
                 lam.1se = list(b0 = sglfit.object$b0[idx1se], beta = sglfit.object$beta[,idx1se]))
  
  obj <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
                cvsd, cvlower = cvm - cvsd, nzero = nz, name = cvname, lamin = lamin, set.seed = set.seed,
              sgl.fit = sglfit.object, cv.fit = cv.fit)
  class(obj) <- "tscv.sglfit"
  obj
} 

#' Sorts time series cross-validation output
#' 
#' @description 
#' Computes \code{cvm} and \code{cvsd} based on times series cross-validation fit. 
#' @usage 
#' tscv.sglpath(outlist, lambda, x, y, foldid, ...)
#' @param outlist list of cross-validation fits.
#' @param lambda a sequence of \eqn{\lambda} parameter.
#' @param x regressors
#' @param y response
#' @param foldid the fold assignment.
#' @param ... other arguments passed to \code{predict.sgl}
#' @return \code{cvm} and \code{cvsd}.
#' @export tscv.sglpath
#' @keywords internal
tscv.sglpath <- function(outlist, lambda, x, y, foldid, ...) {
  typenames <- "Single outcome sg-LASSO"
  y <- as.double(y)
  K <- dim(foldid)[2]
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(K)
  for (i in seq(K)) {
    whichfold <- foldid[,i]
    fitobj <- outlist[[i]]
    preds <- predict.sglpath(fitobj, x[whichfold, , drop = FALSE], ...)
    nlami <- length(outlist[[i]]$lambda)
    predmat[whichfold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- (y-predmat)^2
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames)
} 

computegapobs <- function(whichfoldnot, N, l){
  band <- cbind(whichfoldnot-l,whichfoldnot+l)
  vec <- 1:N
  bandvec <- function(x){
    return(x[1]:x[2])
  }
  exclude <- as.vector(apply(band,1, bandvec))
  exclude <- unique(sort(exclude))
  out <- vec[!vec%in%exclude]
  out
}
