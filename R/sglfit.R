#' Fits sg-LASSO regression
#' 
#' @description 
#' Fits sg-LASSO regression model.
#' The function fits sg-LASSO regression model for a sequence of \eqn{\lambda} tuning parameter and fixed \eqn{\gamma} tuning parameter. The optimization is based on block coordinate-descent. Optionally, fixed effects are fitted. 
#' @details
#' \ifelse{html}{\out{The sequence of linear regression models implied by &lambda; vector is fit by block coordinate-descent. The objective function is <br><br> <center> ||y - &iota;&alpha; - x&beta;||<sup>2</sup><sub>T</sub> + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where ||u||<sup>2</sup><sub>T</sub>=&#60;u,u&#62;/T is the empirical inner product. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on  &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)|&beta;|<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{The sequence of linear regression models implied by \eqn{\lambda} vector is fit by block coordinate-descent. The objective function is \deqn{\|y-\iota\alpha - x\beta\|^2_T + 2\lambda * \Omega_\gamma(\beta),} where \eqn{\|u\|^2_T = \langle u,u \rangle / T} is the empirical inner product. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)|\beta|_{2,1},} a convex combination of LASSO and group LASSO penalty functions.}     
#' @usage 
#' sglfit(x, y, gamma = 1.0, nlambda = 100L, method = c("single", "pooled", "fe"), 
#'        nf = NULL, lambda.factor = ifelse(nobs < nvars, 1e-02, 1e-04), 
#'        lambda = NULL, pf = rep(1, nvars), gindex = 1:nvars, 
#'        dfmax = nvars + 1, pmax = min(dfmax * 1.2, nvars), standardize = TRUE, 
#'        intercept = TRUE, eps = 1e-08, maxit = 1000000L, peps = 1e-08)
#' @param x T by p data matrix, where t and p respectively denote the sample size and the number of regressors.
#' @param y T by 1 response variable.
#' @param gamma sg-LASSO mixing parameter. \eqn{\gamma} = 1 gives LASSO solution and \eqn{\gamma} = 0 gives group LASSO solution.
#' @param nlambda number of \eqn{\lambda}'s to use in the regularization path; used if \code{lambda = NULL}.
#' @param method choose between 'single', 'pooled' and 'fe'; 'single' implies standard sg-LASSO regression, 'pooled' forces the intercept to be fitted, 'fe' computes the fixed effects. User needs to input the number of fixed effects \code{nf}. Default is set to 'single'.
#' @param nf number of fixed effects. Used only if \code{method = 'fe'}.
#' @param lambda.factor The factor for getting the minimal \eqn{\lambda} in the \eqn{\lambda} sequence, where \code{min(lambda) = lambda.factor * max(lambda)}. max(lambda) is the smallest value of lambda for which all coefficients are zero. \eqn{\lambda_{max}} is determined for each \eqn{\gamma} tuning parameter separately. The default depends on the relationship between \code{T} (the sample size) and \code{p} (the number of predictors). If \code{T < p}, the default is \code{0.01}. If \code{T > p}, the default is \code{0.0001}, closer to zero. The smaller the value of \code{lambda.factor} is, the denser is the fit for \ifelse{html}{\out{&lambda;<sub>min</sub>}}{\eqn{\lambda_{min}}}. Used only if \code{lambda = NULL}.
#' @param lambda a user-supplied lambda sequence. By leaving this option unspecified (recommended), users can have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor.} It is better to supply, if necessary, a decreasing sequence of lambda values than a single (small) value, as warm-starts are used in the optimization algorithm. The program will ensure that the user-supplied \eqn{\lambda} sequence is sorted in decreasing order before fitting the model.
#' @param pf \ifelse{html}{\out{&#8467;<sub>1</sub>}}{\eqn{\ell_1}}penalty factor of length \code{p} used for the adaptive sg-LASSO. Separate L1 penalty weights can be applied to each coefficient to allow different \ifelse{html}{\out{&#8467;<sub>1</sub>}}{\eqn{\ell_1}} + \ifelse{html}{\out{&#8467;<sub>2,1</sub>}}{\eqn{\ell_{2,1}}} shrinkage. Can be 0 for some variables, which imposes no shrinkage, and results in that variable always be included in the model. Default is 1 for all variables.
#' @param gindex p by 1 vector indicating group membership of each covariate.
#' @param dfmax the maximum number of variables allowed in the model. Useful for very large \code{p} when a partial path is desired. Default is \code{p+1}. In case \code{method='fe'}, \code{dfmax} is ignored.
#' @param pmax the maximum number of coefficients allowed ever to be nonzero. For example, once \ifelse{html}{\out{&beta;<sub>i</sub> &#8800; 0}}{\eqn{\beta_i \neq 0}}  for some \ifelse{html}{\out{i &#8712; [p]}}{\eqn{i\in[p]}}, no matter how many times it exits or re-enters the model through the path, it will be counted only once. Default is \code{min(dfmax*1.2, p)}.
#' @param standardize logical flag for variable standardization, prior to fitting the model sequence. The coefficients are always returned to the original scale. It is recommended to keep \code{standardize=TRUE}. Default is \code{TRUE}.
#' @param intercept whether intercept be fitted (\code{TRUE}) or set to zero (\code{FALSE}). Default is \code{TRUE}. In case \code{method='pooled'}, \code{intercept=TRUE} is forced. In case \code{method='fe'}, \code{intercept=FALSE} is forced and \code{entity} specific intercepts are fitted in a separate output variable \code{a0}.
#' @param eps convergence threshold for block coordinate descent. Each inner block coordinate-descent loop continues until the maximum change in the objective after any coefficient update is less than thresh times the null deviance. Defaults value is \code{1e-8}.
#' @param maxit maximum number of outer-loop iterations allowed at fixed lambda values. Default is \code{1e6}. If the algorithm does not converge, consider increasing \code{maxit}.
#' @param peps convergence threshold for proximal map of sg-LASSO penalty. Each loop continues until G group difference sup-norm, \ifelse{html}{\out{|| &beta;<sup>k</sup><sub>G</sub> - &beta;<sup>k-1</sup><sub>G</sub> ||<sub>&#8734;</sub>}}{\eqn{\|\beta^{k}_{G} - \beta^{k-1}_{G} \|_\infty}}, is less than \code{peps}. Defaults value is \code{1e-8}.
#' @return sglfit object.
#' @author Jonas Striaukas
#' @examples
#' \donttest{ 
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' beta = c(5,4,3,2,1,rep(0, times = 15))
#' y = x%*%beta + rnorm(100)
#' gindex = sort(rep(1:4,times=5))
#' sglfit(x = x, y = y, gindex = gindex, gamma = 0.5)
#' }
#' @export sglfit 
sglfit <- function(x, y, gamma = 1.0, nlambda = 100L, method = c("single", "pooled", "fe"), nf = NULL,
                  lambda.factor = ifelse(nobs < nvars, 1e-02, 1e-04), 
                  lambda = NULL, pf = rep(1, nvars),
                  gindex = 1:nvars, dfmax = nvars + 1, 
                  pmax = min(dfmax * 1.2, nvars), standardize = TRUE, 
                  intercept = TRUE, eps = 1e-08, maxit = 1000000L, peps = 1e-08) {
    #################################################################################
    ## data setup
    method <- match.arg(method)
    this.call <- match.call()
    y <- drop(y)
    x <- as.matrix(x)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    vnames <- colnames(x)
    ngroups <- as.integer(max(gindex))
    if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")
    if (NROW(y) != nobs) stop("x and y have different number of observations")
    if (NCOL(y) > 1L) stop("Multivariate response is not supported now")
    #################################################################################
    ## parameter setup
    if (length(pf) != nvars) 
      stop("Size of L1 penalty factors does not match the number of input variables")
    if (length(gindex) != nvars) 
      stop("Group index does not match the number of input variables")
    maxit <- as.integer(maxit)
    pf <- as.double(pf)
    gindex <- as.integer(gindex)
    isd <- as.integer(standardize)
    intr <- as.integer(intercept)
    eps <- as.double(eps)
    peps <- as.double(peps)
    dfmax <- as.integer(dfmax)
    pmax <- as.integer(pmax)
    jd <- as.integer(0)
    #################################################################################
    # panel regression setup
    if (method == "single"){
      nf <- as.integer(0)
    }
    if (method != "single"){
      isd <- as.integer(standardize)
      if (method == "pooled"){
          intr <- as.integer(1)
          nf <- as.integer(0)
      }
      if (method == "fe"){
        intr <- as.integer(0)
        if (is.null(nf)) 
          stop("Mehtod set as fixed effects without specifying the number of fixed effects (nf).")
          
        T <- nobs/nf
        if (round(T) != T)
          stop("Number of fixed effects (nf) is not a multiple of the time series dimension, ie nf * T != nobs.")
      }
    }
    # lambda setup
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
      if (lambda.factor >= 1) stop("lambda.factor should be less than 1")
      flmin <- as.double(lambda.factor)
      ulam <- double(nlambda)
      ulam[1] <- -1
      ulam <- as.double(ulam)
    } else {
        flmin <- as.double(1)
        if (any(lambda < 0)) stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    #################################################################################
    fit <- sglfitpath(x, y, nlam, flmin, ulam, isd, intr, nf, eps, peps, dfmax, pmax, jd, 
                pf, gindex, ngroups, maxit, gamma, nobs, nvars, vnames)
    fit$call <- this.call
    #################################################################################
    class(fit) <- c("sglfit", class(fit))
    fit
} 
