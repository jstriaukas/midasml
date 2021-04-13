#' Computes prediction 
#' @description 
#' Similar to other predict methods, this functions predicts fitted values from a fitted sglfit object.
#' @details 
#' \code{s} is the new vector at which predictions are to be made. If s is not in the lambda sequence used for fitting the model, the predict function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \eqn{lambda} indices.
#' @param object fitted \code{\link{sglfit}} model object.
#' @param newx matrix of new values for x at which predictions are to be made. NOTE: \code{newx} must be a matrix, predict function does not accept a vector or other formats of newx.
#' @param s value(s) of the penalty parameter \eqn{lambda} at which predictions are to be made. Default is the entire sequence used to create the model.
#' @param type type of prediction required. Only response is available. Gives predicted response for regression problems.
#' @param method choose between 'single', 'pooled', and 'fe'. 
#' @param ... Not used. Other arguments to predict.
#' @return The object returned depends on type.
#' @export predict.sglpath
predict.sglpath <- function(object, newx, s = NULL, type = c("response"), method = c("single","pooled","fe"), ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  if (method == "single" || method == "pooled"){
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- methods::rbind2(b0, object$beta)
    if (!is.null(s)) {
      vnames <- dimnames(nbeta)[[1]]
      dimnames(nbeta) <- list(NULL, NULL)
      lambda <- object$lambda
      lamlist <- lambda.interp(lambda, s)
      nbeta <- nbeta[,lamlist$left,drop=FALSE]%*%Matrix::Diagonal(x=lamlist$frac) +
        nbeta[,lamlist$right,drop=FALSE]%*%Matrix::Diagonal(x=1-lamlist$frac)
      dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    nfit <- as.matrix(as.matrix(methods::cbind2(1, newx)) %*% nbeta)
  } 
  if (method == "fe"){
    a0 <- object$a0
    N <- object$nf
    T <- dim(newx)[1]/N
    nbeta <- methods::rbind2(a0, object$beta)
    if (!is.null(s)) {
      vnames <- dimnames(nbeta)[[1]]
      dimnames(nbeta) <- list(NULL, NULL)
      lambda <- object$lambda
      lamlist <- lambda.interp(lambda, s)
      nbeta <- nbeta[,lamlist$left,drop=FALSE]%*%Matrix::Diagonal(x=lamlist$frac) +
        nbeta[,lamlist$right,drop=FALSE]%*%Matrix::Diagonal(x=1-lamlist$frac)
      dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    nfit <- as.matrix(as.matrix(methods::cbind2(kronecker(diag(N), rep(1, times=T)),newx)) %*% nbeta)
  }
  nfit
} 