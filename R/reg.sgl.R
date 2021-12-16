#'  Fit for sg-LASSO regression
#' 
#' @description 
#' Fits sg-LASSO regression model.
#' 
#' The function fits sg-LASSO regression based on chosen tuning parameter selection \code{method_choice}. Options include cross-validation and information criteria. 
#' @details
#' \ifelse{html}{\out{The sequence of linear regression models implied by &lambda; vector is fit by block coordinate-descent. The objective function is  <br><br> <center> ||y - &iota;&alpha; - x&beta;||<sup>2</sup><sub>T</sub> + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where &iota;&#8712;R<sup>T</sup>enter> and ||u||<sup>2</sup><sub>T</sub>=&#60;u,u&#62;/T is the empirical inner product. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on  &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)|&beta;|<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{The sequence of linear regression models implied by \eqn{\lambda} vector is fit by block coordinate-descent. The objective function is \deqn{\|y-\iota\alpha - x\beta\|^2_{T} + 2\lambda \Omega_\gamma(\beta),} where \eqn{\iota\in R^T} and \eqn{\|u\|^2_T = \langle u,u \rangle / T} is the empirical inner product. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)|\beta|_{2,1},} a convex combination of LASSO and group LASSO penalty functions.}     
#' @usage 
#' reg.sgl(x, y, gamma = NULL, gindex, intercept = TRUE, 
#'         method_choice = c("tscv","ic","cv"), verbose = FALSE, ...)
#' @param x T by p data matrix, where T and p respectively denote the sample size and the number of regressors.
#' @param y T by 1 response variable.
#' @param gamma sg-LASSO mixing parameter. \eqn{\gamma} = 1 gives LASSO solution and \eqn{\gamma} = 0 gives group LASSO solution.
#' @param gindex p by 1 vector indicating group membership of each covariate.
#' @param intercept whether intercept be fitted (\code{TRUE}) or set to zero (\code{FALSE}). Default is \code{TRUE}.
#' @param method_choice choose between \code{tscv} \code{ic} and \code{cv}. \code{tscv} fits sg-LASSO based on time series cross-validation (see \link{tscv.sglfit}), \code{ic} fits sg-LASSO based on information criteria (BIC, AIC or AICc, see \link{ic.sglfit}), \code{cv} fits sg-LASSO based on cross-validation (see \link{cv.sglfit}). Additional arguments for each method choice are passed on to the relevant functions. 
#' @param verbose flag to print information.
#' @param ... Other arguments that can be passed to \code{sglfit}.
#' @return reg.sgl object.
#' @author Jonas Striaukas
#' @examples
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' beta = c(5,4,3,2,1,rep(0, times = 15))
#' y = x%*%beta + rnorm(100)
#' gindex = sort(rep(1:4,times=5))
#' reg.sgl(x = x, y = y, gamma = 0.5, gindex = gindex)
#' @export reg.sgl
reg.sgl <- function(x, y, gamma = NULL, gindex, intercept = TRUE, method_choice = c("tscv","ic","cv"), verbose = FALSE, ...){
  if(any(is.na(y)))
    stop("y has NA entries, check and rerun")
  if(any(is.na(x)))
    stop("X has NA entries, check and rerun")
  # get settings
  method_choice <- match.arg(arg = method_choice, choices = c("tscv","ic","cv"))
  
  tp <- dim(x)
  t <- tp[1]
  p <- tp[2]
  
  if (is.null(gamma)){
    gamma <- 0.8 
    warning("gamma parameter is set to 0.8 as it is un-supplied.")
  }
  if (method_choice=="ic"){
    if(verbose)
      message(paste0("computing solution using information criteria"))
    if (length(gamma) == 1){
      fit <- ic.sglfit(x, y, gamma = gamma, gindex = gindex, intercept = intercept, ...)
    } else {
      cvms <- matrix(0, nrow = numeric(length(gamma)), ncol = 3)
      for (j in seq(length(gamma))){
        cvms[j, ] <- apply(ic.sglfit(x, y, gamma = gamma[j], gindex = gindex, intercept = intercept, ...)$cvm, 2, min)
      }
      mins <- apply(cvms, 2, min)
      fit <- NULL
      for (i in seq(3)){
        idx <- which(cvms == mins[i], arr.ind = T)[1]
        fit[[i]] <- ic.sglfit(x, y, gamma = gamma[idx], gindex = gindex, intercept = intercept, ...)
      }
    }
  }
  if (method_choice=="cv") {
    if(verbose)
      message(paste0("computing cross-validation fit"))
    
    if (length(gamma) == 1){
      fit <- cv.sglfit(x, y, gamma = gamma, gindex = gindex, intercept = intercept, ...)
    } else {
      cvms <- numeric(length(gamma))
      for (j in seq(length(gamma))){
        cvms[j] <- min(cv.sglfit(x, y, gamma = gamma[j], gindex = gindex, intercept = intercept, ...)$cvm)
      }
      idx <- which(cvms == min(cvms))
      fit <- cv.sglfit(x, y, gamma = gamma[idx], gindex = gindex, intercept = intercept, ...)
    }
  }
  if (method_choice=="tscv") {
    if(verbose)
      message(paste0("computing time series cross-validation fit"))
    
    if (length(gamma) == 1){
      fit <- cv.sglfit(x, y, gamma = gamma, gindex = gindex, intercept = intercept, ...)
    } else {
      cvms <- numeric(length(gamma))
      for (j in seq(length(gamma))){
        cvms[j] <- min(tscv.sglfit(x, y, gamma = gamma[j], gindex = gindex, intercept = intercept, ...)$cvm)
      }
      idx <- which(cvms == min(cvms))
      fit <- tscv.sglfit(x, y, gamma = gamma[idx], gindex = gindex, intercept = intercept, ...)
    }
  }
  class(fit) <- c("reg.sgl", class(fit))
  fit
}

