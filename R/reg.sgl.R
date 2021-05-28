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
#'         method_choice = c("ic","cv"), nfolds = 10, 
#'         verbose = FALSE, ...)
#' @param x T by p data matrix, where T and p respectively denote the sample size and the number of regressors.
#' @param y T by 1 response variable.
#' @param gamma sg-LASSO mixing parameter. \eqn{\gamma} = 1 gives LASSO solution and \eqn{\gamma} = 0 gives group LASSO solution.
#' @param gindex p by 1 vector indicating group membership of each covariate.
#' @param intercept whether intercept be fitted (\code{TRUE}) or set to zero (\code{FALSE}). Default is \code{TRUE}.
#' @param method_choice choose between \code{ic} and \code{cv}. \code{ic} gives fit based on information criteria (BIC, AIC or AICc) by running \code{ic.fit}, while \code{cv} gives fit based on cross-validation by running \code{cv.sglfit}. If \code{cv} is chosen, optional number of folds \code{nfolds} can be supplied. 
#' @param nfolds number of folds of the cv loop. Default set to \code{10}.
#' @param verbose flag to print information.
#' @param ... Other arguments that can be passed to \code{sglfit}.
#' @return reg.sgl object.
#' @author Jonas Striaukas
#' @examples
#' \donttest{ 
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' beta = c(5,4,3,2,1,rep(0, times = 15))
#' y = x%*%beta + rnorm(100)
#' gindex = sort(rep(1:4,times=5))
#' reg.sgl(x = x, y = y, gindex = gindex, gamma = 0.5,
#'   standardize = FALSE, intercept = FALSE)
#' }
#' @export reg.sgl
reg.sgl <- function(x, y, gamma = NULL, gindex, intercept = TRUE, method_choice = c("ic","cv"), nfolds = 10, verbose = FALSE, ...){
  if(any(is.na(y)))
    stop("y has NA entries, check and rerun")
  if(any(is.na(x)))
    stop("X has NA entries, check and rerun")
  # get settings
  method_choice <- match.arg(method_choice)
  
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
    
    fit <- ic.sglfit(x, y, gamma = gamma, gindex = gindex, intercept = intercept, method = "single", ...)
  }
  if (method_choice=="cv") {
    if(verbose)
      message(paste0("computing ", nfolds, "-fold cross-validation"))
    
    fit <- cv.sglfit(x, y, gamma = gamma, gindex = gindex, nfolds = nfolds, intercept = intercept, method = "single", ...)
  }
  class(fit) <- c("reg.sgl", class(fit))
  fit
}

