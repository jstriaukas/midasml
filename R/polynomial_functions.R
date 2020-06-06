#' Legendre polynomials shifted to [a,b]
#' 
#' @description For a given set of points in X, computes the orthonormal Legendre polynomials basis of L2 [a,b] for a given degree.
#' @param degree polynomial degree.
#' @param a lower shift value (default - 0).
#' @param b upper shift value (default - 1).
#' @param jmax number of high-frequency lags.
#' @param X optional evaluation grid vector.
#' @return Psi weight matrix with Legendre functions upto \code{degree}.
#' @export lb
lb <- function(degree,a=0,b=1,jmax=NULL,X=NULL){
  # Notes:
  #   References:
  #   H. J. Weber, G. B. Arfken, Essential Mathematical Methods for Physicists,
  #   Elsevier Academic Press, San Diego, USA, 2004.
  #   Translated from matlab function lb to R. 2019-03-01, Jonas Striaukas (jonas.striaukas@gmail.com)
  if (!is.null(jmax)){
    X <- seq(0,1,length.out=jmax)
  }
  if (is.null(X)){
    stop("X is not provided. Either set X or set jmax.")
  }
  n <- length(X)
  P <- matrix(1,nrow=n,ncol=degree+2)
  Psi <- matrix(1,nrow=n,ncol=degree+1)/ sqrt(b-a)
  P[, 2] <-  2*X/(b-a) - ((b+a)/(b-a))
  
  for (i in 1:degree){
    P[, i+2]   <- ((2*i+1)/(i+1)) * P[, 2]*P[, i+1] - i/(i+1) * P[, i]
    Psi[, i+1] <- sqrt((2*i + 1) / (b-a)) %*% P[, i+1]
  }
  
  return(Psi)
}

#' Gegenbauer polynomials shifted to [a,b]
#' 
#' @description For a given set of points in X, computes the orthonormal Gegenbauer polynomials basis of L2 [a,b] for a given degree and \eqn{\alpha} parameter. The Gegenbauer polynomials are special cases of Jacobi polynomials. In turn, you may get Legendre polynomials from Gegenbauer by setting \eqn{\alpha} = 0, or Chebychev's polynomials 
#'  by setting \eqn{\alpha} = 1/2 or -1/2.
#' @param degree polynomial degree.
#' @param alpha Gegenbauer polynomials parameter.
#' @param a lower shift value (default - 0).
#' @param b upper shift value (default - 1).
#' @param jmax number of high-frequency lags.
#' @param X optional evaluation grid vector.
#' @return Psi weight matrix with Gegenbauer functions upto \code{degree}.
#' @export gb
gb <- function(degree,alpha,a=0,b=1,jmax=NULL,X=NULL){
  if (!is.null(jmax)){
    X <- seq(0,1,length.out=jmax)
  }
  if (is.null(X)){
    stop("X is not provided. Either set X or set jmax.")
  }
  n <- length(X)
  P <- matrix(1,nrow=n,ncol=degree+2)
  Psi <- matrix(1,nrow=n,ncol=degree+1)/ sqrt(b-a)
  P[, 2] <-  2*X/(b-a) - ((b+a)/(b-a))
  
  for (i in 1:degree){
    a <- (2*i + 2*alpha + 1)*(2*i + 2*alpha + 2)/(2*(i+1)*(i + 2*alpha + 1))
    c <- (alpha + i)^2*(2*i + 2*alpha + 2)/( (i+1)*(i + 2*alpha + 1)*(2*i + 2*alpha)) 
    P[, i+2]   <- a * P[, 2]*P[, i+1] - c * P[, i]
    Psi[, i+1] <- sqrt((2*i + 1) / (b-a)) %*% P[, i+1]
  }
  return(Psi)
}

#' Exponential Almon polynomial weights 
#' 
#' @param param two-dimensional parameter vector \eqn{\theta}.
#' @param dayLag number of high-frequency lags.
#' @return (normalized) weights vector.
#' @description For a given set of parameters \eqn{\theta} and the number of high-frequency lags, returns weights implied by exponential Almon functional form.
#' @export expalmon_w
expalmon_w <- function(param,dayLag){
  grid <- seq(1,dayLag,length.out = dayLag)
  w <- exp(param[1]*grid+param[2]*grid^2)/sum(exp(param[1]*grid+param[2]*grid^2))
  w <- w/sum(w)
  w
}

#' Beta density polynomial weights 
#' 
#' @param param two-dimensional parameter vector \eqn{\theta}.
#' @param dayLag number of high-frequency lags.
#' @return (normalized) weights vector.
#' @description For a given set of parameters \eqn{\theta} and the number of high-frequency lags, returns weights implied by Beta density functional form.
#' @export expalmon_w
beta_w <- function(param,dayLag){
  eps <- .Machine$double.eps
  u <- seq(eps,1-eps,length.out = dayLag)
  w <- u^(param[1]-1)*(1-u)^(param[2]-1)
  w <- w/sum(w)
  w
}

#' Restricted Beta density polynomial weights 
#' 
#' @param param one-dimensional parameter vector \eqn{\theta}.
#' @param dayLag number of high-frequency lags.
#' @return (normalized) weights vector.
#' @description For a given set of parameters \eqn{\theta} and the number of high-frequency lags, returns weights implied by Restricted Beta density functional form.
#' @export expalmon_w
rbeta_w <- function(param,dayLag){
  eps <- .Machine$double.eps
  u <- seq(eps,1-eps,length.out = dayLag)
  w <- (1-u)^(param[1]-1)
  w <- w/sum(w)
  w
}