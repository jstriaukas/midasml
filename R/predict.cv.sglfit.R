#' Computes prediction 
#' @description 
#' Similar to other predict methods, this functions predicts fitted values from a fitted sglfit object.
#' @details 
#' \code{s} is the new vector at which predictions are to be made. If s is not in the lambda sequence used for fitting the model, the predict function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \eqn{lambda} indices.
#' @param object fitted \code{\link{cv.sglfit}} model object.
#' @param newx matrix of new values for x at which predictions are to be made. NOTE: \code{newx} must be a matrix, predict function does not accept a vector or other formats of newx.
#' @param s choose between 'lam.min' and 'lam.1se'.
#' @param type type of prediction required. Only response is available. Gives predicted response for regression problems.
#' @param method choose between 'single', 'pooled', and 'fe'. 
#' @param ... Not used. Other arguments to predict.
#' @return The object returned depends on type.
#' @export predict.cv.sglfit
predict.cv.sglfit <- function(object, newx, s = c("lam.min","lam.1se"), type = c("response"), method = c("single","pooled","fe"),...) {
  type <- match.arg(type)
  method <- match.arg(method)
  s <- match.arg(s)
  if (s == "lam.min") {
    if (method == "single"){
      object <- object$cv.fit$lam.min
    } else{
      object <- object$cv.panel.fit$lam.min
    }
  }
  if (s == "lam.1se") {
    if (method == "single"){
      object <- object$cv.fit$lam.1se
    } else{
      object <- object$cv.panel.fit$lam.1se
    }
  }
  if (method == "single" || method == "pooled"){
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- c(b0, object$beta)
    nfit <- c(1, newx) %*% nbeta
  } 
  if (method == "fe"){
    a0 <- object$a0
    N <- object$nf
    T <- dim(newx)[1]/N
    nbeta <- object$beta
    nfit <- newx %*% nbeta + a0
  }
  nfit
} 