#' Information criteria fit for panel sg-LASSO
#' 
#' @description 
#' Does information criteria for panel data sg-LASSO regression model.
#' 
#' The function runs \link{sglfit} 1 time; computes the path solution in \ifelse{html}{\out{<code>lambda</code>}}{\code{lambda}} sequence.
#' Solutions for \code{BIC}, \code{AIC} and \code{AICc} information criteria are returned. 
#' @details
#' \ifelse{html}{\out{The sequence of linear regression models implied by &lambda; vector is fit by block coordinate-descent. The objective function is (case <code>method='pooled'</code>) <br><br> <center> ||y - &iota;&alpha; - x&beta;||<sup>2</sup><sub>NT</sub> + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where &iota;&#8712;R<sup>NT</sup> and &alpha; is common intercept to all N items or (case <code>method='fe'</code>) <br><br> <center> ||y - B&alpha; - x&beta;||<sup>2</sup><sub>NT</sub> + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where B = I<sub>N</sub>&#8855;&iota; and ||u||<sup>2</sup><sub>NT</sub>=&#60;u,u&#62;/NT is the empirical inner product. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on  &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)|&beta;|<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{The sequence of linear regression models implied by \eqn{\lambda} vector is fit by block coordinate-descent. The objective function is  either (case \code{method='pooled'})  \deqn{\|y-\iota\alpha - x\beta\|^2_{NT} + 2\lambda \Omega_\gamma(\beta),} where \eqn{\iota\in R^{NT}} and \eqn{\alpha} is common intercept to all N items or (case \code{method='fe'}) \deqn{\|y-B\alpha - x\beta\|^2_{NT} + 2\lambda \Omega_\gamma(\beta),} where \eqn{B=I_N \times \iota} and \eqn{\|u\|^2_{NT} = \langle u,u \rangle / NT} is the empirical inner product. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)|\beta|_{2,1},} a convex combination of LASSO and group LASSO penalty functions.}     
#' @usage 
#' ic.panel.sglfit(x, y, lambda = NULL, gamma = 1.0, gindex = 1:p, 
#'                 method = c("pooled","fe"), nf = NULL, ...)
#' @param x NT by p data matrix, where NT and p respectively denote the sample size of pooled data and the number of regressors.
#' @param y NT by 1 response variable.
#' @param lambda a user-supplied lambda sequence. By leaving this option unspecified (recommended), users can have the program compute its own \eqn{\lambda} sequence based on \code{nlambda} and \code{lambda.factor.} It is better to supply, if necessary, a decreasing sequence of lambda values than a single (small) value, as warm-starts are used in the optimization algorithm. The program will ensure that the user-supplied \eqn{\lambda} sequence is sorted in decreasing order before fitting the model.
#' @param gamma sg-LASSO mixing parameter. \eqn{\gamma} = 1 gives LASSO solution and \eqn{\gamma} = 0 gives group LASSO solution.
#' @param gindex p by 1 vector indicating group membership of each covariate.
#' @param method choose between 'pooled' and 'fe'; 'pooled' forces the intercept to be fitted in \link{sglfit}, 'fe' computes the fixed effects. User must input the number of fixed effects \code{nf} for \code{method = 'fe'}. Default is set to \code{method = 'pooled'}.
#' @param nf number of fixed effects. Used only if \code{method = 'fe'}.
#' @param ... Other arguments that can be passed to \link{sglfit}.
#' @return ic.panel.sglfit object.
#' @author Jonas Striaukas
#' @examples
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' beta = c(5,4,3,2,1,rep(0, times = 15))
#' y = x%*%beta + rnorm(100)
#' gindex = sort(rep(1:4,times=5))
#' ic.panel.sglfit(x = x, y = y, gindex = gindex, gamma = 0.5, 
#'   standardize = FALSE, intercept = FALSE)
#' @export ic.panel.sglfit
ic.panel.sglfit <- function(x, y, lambda = NULL, gamma = 1.0, gindex = 1:p, method = c("pooled","fe"), nf = NULL, ...){
  method <- match.arg(method)
  if (method == "fe" && is.null(nf))
    stop("for fe method nf must supplied.")
  
  p <- ncol(x)
  NT <- nrow(x)
  if(!is.null(nf)) {
    T <- NT/nf
  } else {
    T <- NT
    nf <- 1
  }
  sglfit.object <- sglfit(x, y, gamma = gamma, gindex = gindex, method = method, nf = nf, ...)
  lambda <- sglfit.object$lambda
  nlam <- length(lambda)
  cvm <- matrix(NA, nrow = nlam, ncol = 3)
  yhats <- predict.sglpath(sglfit.object, newx = x, method = method)
  df <- sglfit.object$df+sglfit.object$nf
  sigsqhat <- sum((y-mean(y))^2)/NT
  mse <- colSums((replicate(nlam, as.numeric(y))-yhats)^2)/NT
  cvm[,1] <- mse/sigsqhat + ic.pen("bic", df, NT) 
  cvm[,2] <- mse/sigsqhat + ic.pen("aic", df, NT) 
  cvm[,3] <- mse/sigsqhat + ic.pen("aicc", df, NT) 
  idx <- numeric(3)
  for(i in seq(3)){
    lamin <- lambda[which(cvm[,i]==min(cvm[,i]))]
    if (length(lamin)>1)
      lamin <- min(lamin)
    
    min.crit <- min(cvm[,i])
    idx[i] <- which(lamin==lambda)
  }
  if (method == "pooled"){
    ic.panel.fit <- list(bic.fit = list(b0 = sglfit.object$b0[idx[1]], beta = sglfit.object$beta[,idx[1]]),
                      aic.fit = list(b0 = sglfit.object$b0[idx[2]], beta = sglfit.object$beta[,idx[2]]),
                      aicc.fit = list(b0 = sglfit.object$b0[idx[3]], beta = sglfit.object$beta[,idx[3]]))
  }
  if (method == "fe"){
    ic.panel.fit <- list(bic.fit = list(a0 = sglfit.object$a0[,idx[1]], beta = sglfit.object$beta[,idx[1]]),
                      aic.fit = list(a0 = sglfit.object$a0[,idx[2]], beta = sglfit.object$beta[,idx[2]]),
                      aicc.fit = list(a0 = sglfit.object$a0[,idx[3]], beta = sglfit.object$beta[,idx[3]]))
  }
  
  obj <- list(lambda = lambda, cvm = cvm, lamin = lamin, 
              fit = sglfit.object, ic.panel.fit = ic.panel.fit)
  class(obj) <- "ic.panel.sglfit"
  obj
}