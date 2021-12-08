#' MIDAS regression 
#' 
#' @description 
#' Fits MIDAS regression model with single high-frequency covariate. Options include linear-in-parameters polynomials (e.g. Legendre) or non-linear polynomials (e.g. exponential Almon).
#' Nonlinear polynomial optimization routines are equipped with analytical gradients, which allows fast and accurate optimization. 
#' 
#' @details
#' \ifelse{html}{\out{Several polynomial functional forms are available (<code>poly_choice</code>): <br><br> - <code>beta</code>: Beta polynomial</center> <br> - <code>expalmon</code>: exponential Almon polynomial <br> - <code>legendre</code>: Legendre polynomials. <br> <br> The ARDL-MIDAS model is: <br> <center> y<sub>t</sub> =  &mu; + &Sigma;<sub>p</sub> &rho;<sub>p</sub> y<sub>t-p</sub> + &beta;  &Sigma;<sub>j</sub> &omega;<sub>j</sub>(&theta;)x<sub>t-1</sub> </center> <br> where &mu;, &beta;, &theta; and  &rho;<sub>p</sub> are model parameters, p is the number of low-frequency lags and &omega; is the weight function.}}{Several polynomial functional forms are available (\code{poly_choice}): \cr  - \code{beta}: Beta polynomial  \cr - \code{expalmon}: Exp Almon polynomial \cr - \code{legendre}: Legendre polynomials. \cr\cr The ARDL-MIDAS model is: \cr \deqn{y_t =  \mu + \sum_p \rho_p y_{t-p} + \beta \sum_j \omega_j(\theta)x_{t-1}} \cr where \eqn{\mu}, \eqn{\beta}, \eqn{\theta}, \eqn{\rho_p}  are model parameters, p is number of low-frequency and \eqn{\omega} is the weight function.}     
#' @usage 
#' midas.ardl(y, x, z = NULL, loss_choice = c("mse","logit"), 
#'            poly_choice = c("legendre","expalmon","beta"), 
#'            poly_spec = 0, legendre_degree = 3, nbtrials = 500)
#' @param y response variable. Continuous for \code{loss_choice = "mse"}, binary for \code{loss_choice = "logit"}.
#' @param x high-frequency covariate lags. 
#' @param z other lower-frequency covariate(s) or AR lags (both can be supplied in an appended matrix). Either must be supplied. 
#' @param loss_choice which loss function to fit: \code{loss_choice="mse"} fits least squares MIDAS regression, \code{loss_choice="logit"} fits logit MIDAS regression. 
#' @param poly_choice which MIDAS lag polynomial function to use: \code{poly_choice="expalmon"} - exponential Almon polynomials, \code{poly_choice="beta"} - Beta density function (need to set \code{poly_spec}), \code{poly_choice="legendre"} - legendre polynomials (need to set \code{legendre_degree}). Default is set to \code{poly_choice="expalmon"}.
#' @param poly_spec which Beta density function specification to apply (applicable only for \code{poly_choice="beta"}). \code{poly_spec = 0} - all three parameters are fitted,  \code{poly_spec = 1} (\eqn{\theta_2,\theta_3}) are fitted, \code{poly_spec = 2} (\eqn{\theta_1,\theta_2}) are fitted, \code{poly_spec = 3} (\eqn{\theta_2}) is fitted. Default is set to \code{poly_spec = 0}.
#' @param legendre_degree the degree of legendre polynomials (applicable only for \code{legendre="beta"}). Default is set to 3.
#' @param nbtrials number of initial values tried in multistart optimization. Default is set to \code{poly_spec = 500}.  
#' @return midas.ardl object.
#' @author Jonas Striaukas
#' @examples
#' set.seed(1)
#' x = matrix(rnorm(100 * 20), 100, 20)
#' z = rnorm(100)
#' y = rnorm(100)
#' midas.ardl(y = y, x = x, z = z)
#' @export midas.ardl
midas.ardl <- function(y, x, z = NULL, loss_choice = c("mse","logit"), poly_choice = c("legendre","expalmon","beta"), poly_spec = 0, legendre_degree = 3, nbtrials = 500){
  loss <- match.arg(loss_choice)
  poly <- match.arg(poly_choice)
  if (is.null(z))
    stop("z must be supplied. Either AR lags or other variables.")

  if (loss == "mse"){
    if (poly == "expalmon"){
      fit <- optim_ardl_expalmon(y = y, z = z, x = x, nbtrials = nbtrials)
    }
    if (poly == "beta"){
      fit <- optim_ardl_beta(y = y, z = z, x = x, poly_spec = poly_spec, nbtrials = nbtrials)
    }
    if (poly == "legendre"){
      fit <- optim_ardl_legendre(y = y, z = z, x = x, legendre_degree = legendre_degree)
    }
  } else if (loss == "logit") {
    if (poly == "expalmon"){
      fit <- optim_logit_expalmon(y = y, z = z, x = x)
    }
    if (poly == "legendre"){
      fit <- optim_logit_legendre(y = y, z = z, x = x, legendre_degree = legendre_degree)
    }
  }
  class(fit) <- c("midas.ardl",class(fit))
  fit
}








