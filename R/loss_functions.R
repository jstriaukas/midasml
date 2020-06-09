#' Mean squared error loss for ARDL-MIDAS or DL-MIDAS regression
#' 
#' @description 
#'  For a set of parameter in \code{param} vector, data inputs y, x and z, and weight function specification, computes the mse loss.
#' @param param parameter vector.
#' @param y response variable. 
#' @param x high-frequency covariate lags.
#' @param z low-frequency variable lags.
#' @param weight non-linear weight function specification.
#' @param ... currently ignored optional parameters. 
#' @return objective function value evaluated at \code{param}.
#' @export mse_loss
#' @keywords internal
mse_loss <- function(param, y, x, z, weight, ...){
  d <- dim(x)[2]
  n <- length(y)
  iota <- rep(1,times=n)
  c <- param[1]
  param <- param[-1]
  if(!is.null(z)){
    z <- matrix(z, nrow=n)
    p_ar <- dim(z)[2]
    param_ar <- param[1:p_ar]
    param_midas <- param[-c(1:p_ar)]
    yhat <- iota*param[1]+z%*%param_ar+param_midas[1]*x%*%weight(param_midas[-1],d)
  } else {
    param_midas <- param
    yhat <- iota*param[1]+param_midas[1]*x%*%weight(param_midas[-1],d)
  }
  r <- y-yhat
  obj <- 1/n*sum(r^2)
  obj
}

#' Asymmetric least squares loss for ARDL-MIDAS or DL-MIDAS regression
#' 
#' @description 
#'  For a set of parameter in \code{param} vector, data inputs y, x and z, and weight function specification, computes the mse loss.
#' @param param parameter vector.
#' @param y response variable. 
#' @param x high-frequency covariate lags.
#' @param z low-frequency variable lags.
#' @param weight non-linear weight function specification.
#' @param tau quantile level.
#' @param ... currently ignored optional parameters. 
#' @return objective function value evaluated at \code{param}.
#' @export mse_loss
#' @keywords internal
als_loss <- function(param, y, x, z, weight, tau, ...){
  d <- dim(x)[2]
  n <- length(y)
  iota <- rep(1,times=n)
  c <- param[1]
  param <- param[-1]
  if(!is.null(z)){
    z <- matrix(z, nrow=n)
    p_ar <- dim(z)[2]
    param_ar <- param[1:p_ar]
    param_midas <- param[-c(1:p_ar)]
    yhat <- iota*param[1]+z%*%param_ar+param_midas[1]*x%*%weight(param_midas[-1],d)
  } else {
    param_midas <- param
    yhat <- iota*param[1]+param_midas[1]*x%*%weight(param_midas[-1],d)
  }
  r <- y-yhat
  obj <- 1/n*sum(r^2*abs( tau - as.numeric((r < 0)) ))
  obj
}

#' Quantile regression loss for ARDL-MIDAS or DL-MIDAS regression
#' 
#' @description 
#'  For a set of parameter in \code{param} vector, data inputs y, x and z, and weight function specification, computes the mse loss.
#' @param param parameter vector.
#' @param y response variable. 
#' @param x high-frequency covariate lags.
#' @param z low-frequency variable lags.
#' @param weight non-linear weight function specification.
#' @param tau quantile level.
#' @param ... currently ignored optional parameters. 
#' @return objective function value evaluated at \code{param}.
#' @export mse_loss
#' @keywords internal
rq_loss <- function(param, y, x, z, weight, tau, ...){
  d <- dim(x)[2]
  n <- length(y)
  iota <- rep(1,times=n)
  c <- param[1]
  param <- param[-1]
  if(!is.null(z)){
    z <- matrix(z, nrow=n)
    p_ar <- dim(z)[2]
    param_ar <- param[1:p_ar]
    param_midas <- param[-c(1:p_ar)]
    yhat <- iota*param[1]+z%*%param_ar+param_midas[1]*x%*%weight(param_midas[-1],d)
  } else {
    param_midas <- param
    yhat <- iota*param[1]+param_midas[1]*x%*%weight(param_midas[-1],d)
  }
  r <- y-yhat
  obj <- 1/n*sum(r*( tau - as.numeric((r < 0)) ))
  obj
}
