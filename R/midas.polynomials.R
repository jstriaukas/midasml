#' Legendre polynomials shifted to [a,b]
#' 
#' @description For a given set of points in X, computes the orthonormal Legendre polynomials basis of L2 [a,b] for a given degree.
#' @param degree polynomial degree.
#' @param a lower shift value (default - 0).
#' @param b upper shift value (default - 1).
#' @param jmax number of high-frequency lags.
#' @param X optional evaluation grid vector.
#' @return Psi weight matrix with Legendre functions upto \code{degree}.
#' @author Jonas Striaukas
#' @examples
#' degree <- 3
#' jmax <- 66
#' lb(degree = degree, a = 0, b = 1, jmax = jmax)
#' @export lb
lb <- function(degree,a=0,b=1,jmax=NULL,X=NULL){
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
  if(degree>0){
    for (i in 1:degree){
      P[, i+2]   <- ((2*i+1)/(i+1)) * P[, 2]*P[, i+1] - i/(i+1) * P[, i]
      Psi[, i+1] <- sqrt((2*i + 1) / (b-a)) %*% P[, i+1]
    }
  }
  col.name <- NULL
  for (i in 0:degree)
    col.name <- c(col.name, paste0("poly-degree-",i))
    
  row.name <- NULL
  for (i in 1:n)
    row.name <- c(row.name, paste0("lag-",n))
  
  colnames(Psi) <- col.name
  rownames(Psi) <- row.name
  class(Psi) <- "midasml"
  Psi
  return(Psi)
}

#' Gegenbauer polynomials shifted to [a,b]
#' 
#' @description For a given set of points in X, computes the orthonormal Gegenbauer polynomials basis of L2 [a,b] for a given degree and \eqn{\alpha} parameter. The Gegenbauer polynomials are a special case of more general Jacobi polynomials. In turn, you may get Legendre polynomials from Gegenbauer by setting \eqn{\alpha} = 0, or Chebychev's polynomials 
#'  by setting \eqn{\alpha} = 1/2 or -1/2.
#' @param degree polynomial degree.
#' @param alpha Gegenbauer polynomials parameter.
#' @param a lower shift value (default - 0).
#' @param b upper shift value (default - 1).
#' @param jmax number of high-frequency lags.
#' @param X optional evaluation grid vector.
#' @return Psi weight matrix with Gegenbauer functions upto \code{degree}.
#' @author Jonas Striaukas
#' @examples
#' degree <- 3
#' alpha <- 1
#' jmax <- 66
#' gb(degree = degree, alpha = alpha, a = 0, b = 1, jmax = jmax)
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
  if(degree>0){
    for (i in 1:degree){
      d <- (2*i + 2*alpha + 1)*(2*i + 2*alpha + 2)/(2*(i+1)*(i + 2*alpha + 1))
      c <- (alpha + i)^2*(2*i + 2*alpha + 2)/( (i+1)*(i + 2*alpha + 1)*(2*i + 2*alpha)) 
      P[, i+2]   <- d * P[, 2]*P[, i+1] - c * P[, i]
      Psi[, i+1] <- sqrt((2*i + 1) / (b-a)) %*% P[, i+1]
    }
  }
  col.name <- NULL
  for (i in 0:degree)
    col.name <- c(col.name, paste0("poly-degree-",i))
  
  row.name <- NULL
  for (i in 1:n)
    row.name <- c(row.name, paste0("lag-",n))
  
  colnames(Psi) <- col.name
  rownames(Psi) <- row.name
  class(Psi) <- "midasml"
  Psi
  return(Psi)
}