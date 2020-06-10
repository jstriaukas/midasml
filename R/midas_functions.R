#' DL-MIDAS regression  
#' 
#' @description 
#'  Estimates and predicts using a single variate DL-MIDAS model. 
#' @details 
#'  \ifelse{html}{\out{Several polynomial functional forms are available (input variable <code>polynomial</code>): <br><br> - <code>beta_w</code>: Beta polynomial</center> <br> - <code>rbeta_w</code>: restricted Beta polynomial <br> - <code>expalmon_w</code>: exponential Almon polynomial <br>  - <code>umidas_w</code>: unrestricted lags (U-MIDAS) <br> - <code>step_fun</code>: polynomial with step functions <br> - <code>legendre_w</code>: Legendre polynomials <br> <br> different forecasting schemes (input variable <code>scheme</code>): <br><br>  - <code>fixed</code>: fixed scheme <br> - <code>rolling</code>: rolling window <br> - <code>expand</code>: expanding window <br><br> and different loss functions (input variable <code>loss</code>): <br><br> - <code>mse</code>: least squares <br> - <code>als</code>: asymmetric least squares <br> - <code>rq</code>: quantile. <br> <br> The DL-MIDAS model is: <br> <center> y<sub>t</sub> =  &mu; +  &beta;  &Sigma;<sub>j</sub> &omega;<sub>j</sub>(&theta;)x<sub>t-1</sub> </center> <br> where &mu;, &beta; and &theta; are model parameters and &omega; is the weight function.}}{Several polynomial functional forms are available (input variable \code{polynomial}): \cr  - \code{beta_w}: Beta polynomial \cr - \code{rbeta_w}: restricted Beta polynomial \cr - \code{expalmon_w}: Exp Almon polynomial \cr - \code{umidas_w}: unrestricted lags (U-MIDAS) \cr - \code{step_fun}: polynomial with step functions  \cr - \code{legendre_w}: Legendre polynomials \cr different forecasting schemes (input variable \code{scheme}): \cr - \code{fixed}: fixed scheme \cr - \code{rolling}: rolling window \cr - \code{expand}: expanding window \cr  and different loss functions (input variable \code{loss}) \cr - \code{mse}: least squares \cr - \code{als}: asymmetric least squares \cr - \code{rq}: quantile.  \cr\cr The DL-MIDAS model is: \cr \deqn{y_t =  \mu +  \beta  \sum_j \omega_j(\theta)x_{t-1}} \cr where \eqn{\mu}, \eqn{\beta} and \eqn{\theta} are model parameters and \eqn{\omega} is the weight function.}     
#' @usage 
#' midas_dl(data.x, data.xdate, data.y, data.ydate, x.lag, est.start, est.end, horizon = 1,
#'   polynomial = c("legendre_w", "beta_w", "rbeta_w", "expalmon_w", "umidas_w", "step_fun"),
#'   scheme = c("fixed", "rolling", "expand"), loss = c("mse", "rq", "als"), ...)
#' @param data.x predictor variable series.
#' @param data.xdate predictor variable dates.
#' @param data.y dependent variable series (can leave unspecified, see \code{midas_gen} option).
#' @param data.ydate dependent variable dates (can leave unspecified, see \code{midas_gen} option).
#' @param x.lag number of high-frequency data lags.
#' @param est.start start date of the estimation sample (referenced with data.xdate).
#' @param est.end end date of the estimation sample (referenced with data.xdate).
#' @param horizon forecast horizon measured in predictor variable sampling frequency (default set 1 unit ahead).
#' @param polynomial MIDAS lag polynomial specification. Options are: Legendre (legendre_w), Beta desity (beta_w), restricted Beta density (rbeta_w), exponential Almon (expalmon_w), unrestricted MIDAS (umidas_w), step functions (step_fun).
#' @param scheme forecasting scheme. Options are: fixed scheme (fixed), rolling window scheme (rolling), expanding window scheme (expand).
#' @param loss loss function. Options are: mean squared error (mse), quantile (rq), asymmetric least squares (als).
#' @param ... optional parameters to feed into other functions. forecast.flag - TRUE/FALSE to compute out-of-sample predictions (default TRUE), disp.flag - TRUE/FALSE  to display MIDAS data structures (default FALSE), 
#'    num.evals - number of objective function evaluations using random starting parameter values in the case of non-linear MIDAS polynomial (default 1e4), 
#'    num.coef - number of best coefficients to use as starting values in nonlinear optimization (default 10),
#'    seed - value used in set.seed for randomly drawing initial starting values around OLS optimal solution,
#'    profiling - TRUE/FALSE to use MIDAS paramater profiling, coded only for rbeta_w polynomial, (default FALSE),
#'    step_idx - index of step function lags. If step_fun is used as a polynomial, it is best to specify this option too, otherwise, the program figures out the sampling frequency ratio and computes \code{step_idx} accordingly (message is displayed in this case),
#'    legendre_degree - a degree of legendre polynomials. If legendre_w is used as a polynomial, it is best to specify this option too, otherwise, the value is set to 3 (message is displayed in this case),
#'    tau - quantile level for als and rq regressions. If eithr als or rq loss is used, this option must be specified, program stops if not, 
#'    midas_gen - option on how to generate the low-frequency variable. 'from_hf' - computes from high-frequency variable (see \code{mixed_freq_data_mhorizon}, \code{aggregation} method could be specified as an additional input) or 'as_ref' - computes MIDAS data structures using low-frequency variable (default 'from_hf').
#' @return returns midas_dl list which contains parameter estimates, in- and out-of-sample statistics and predictions, and some information about the specification of the method used.
#' @export midas_dl
midas_dl <- function(data.x, data.xdate, data.y, data.ydate, x.lag, est.start, est.end, horizon = 1,
                     polynomial = c("legendre_w","beta_w","rbeta_w","expalmon_w","umidas_w","step_fun"), 
                     scheme = c("fixed","rolling","expand"),
                     loss = c("mse","rq","als"),...){
  polynomial <- match.arg(polynomial)
  scheme <- match.arg(scheme)
  loss <- match.arg(loss)
  # deal with options
  options <- list(...)
  if(is.null(options$forecast.flag)){
    forecast.flag <- TRUE
  } else {
    forecast.flag <- options$forecast.flag
  }
  if(is.null(options$disp.flag)){
    disp.flag <- FALSE
  } else {
    disp.flag <- options$disp.flag
  }
  info <- tau <- legendre_degree <- step_idx <- num.evals <- num.coef <- seed <- profiling <- NULL
  if(polynomial%in%c("beta_w","rbeta_w","expalmon_w")){
    if(is.null(options$num.evals)) {
      num.evals <- 1e4
    } else {
      num.evals <- options$num.evals
    }
    info$num.evals <- num.evals
    if(is.null(options$num.coef)){
      num.coef <- 10
    } else {
      num.coef <- options$num.coef
    }
    info$num.coef <- num.coef
    if(is.null(options$seed)){
      seed <- 100
    } else {
      seed <- options$seed
    }
    info$seed <- seed
  }
  if(polynomial%in%c("rbeta_w")){
    if(is.null(options$profiling)){
      profiling <- FALSE
    } else {
      profiling <- options$profiling
    }
  }
  if(polynomial%in%"step_fun"){
    step_idx <- options$step_idx
    info$step_idx <- step_idx   
  }
  if(polynomial%in%"legendre_w"){
    legendre_degree <- options$legendre_degree
    info$legendre_degree <- legendre_degree   
  }
  if(loss=="als"){
    if(is.null(options$tau))
      stop("set quantile level tau in options")
    tau <- options$tau
    info$tau <- tau
  }
  if(loss=="rq"){
    if(is.null(options$tau))
      stop("set quantile level tau in options")
    tau <- options$tau
    info$tau <- tau
  }
  if(is.null(options$midas_gen)){
    midas_gen <- "from_hf"
  } else {
    midas_gen <- options$midas_gen
  }
  if(midas_gen=="as_ref"){
    data.ydate <- options$data.ydate
    data.y <- options$data.ydate
    if(is.null(data.ydate))
      stop("option to generate MIDAS data with low-frequency reference date was chosen without specifying reference date. Please set 'data.ydate' in options and re-run.")
    if(is.null(data.y))
      stop("option to generate MIDAS data with low-frequency reference date was chosen without low frequency data vector. Please set 'data.y' in options and re-run.")
    
    mf.data <- mixed_freq_data_single(data.refdate = data.ydate, data.x, data.xdate, x.lag, horizon, est.start, est.end, disp.flag) 
    est.y <- data.y[data.ydate%in%mf.data$est.refdate]
    est.ydate <- data.ydate[data.ydate%in%mf.data$out.refdate]
    out.y <- data.y[data.ydate%in%mf.data$est.refdate]
    out.ydate <- data.ydate[data.ydate%in%mf.data$est.refdate]
  } 
  if(midas_gen=="from_hf"){
    mf.data <- mixed_freq_data_mhorizon(data.x, data.xdate, x.lag, est.start, est.end, horizon, disp.flag, aggregation = options$aggregation)     
    est.y <- mf.data$est.y
    est.ydate <- mf.data$est.ydate
    out.y <- mf.data$out.y
    out.ydate <- mf.data$out.ydate
  }
  # in-sample data
  est.x <- mf.data$est.x
  est.xdate <- mf.data$est.xdate
  # out-of-sample data
  out.x <- mf.data$out.x
  out.xdate <- mf.data$out.xdate
  
  info$est.y <- est.y
  info$est.x <- est.x
  info$x.lag <- x.lag
  info$out.y <- out.y
  info$out.x <- out.x
  
  nobs <- length(est.y)
  nforecast <- length(out.y)
  if(nforecast==0){
    if (forecast.flag)
      message("'forecast.flag' was set TRUE but the end of of the estimation sample coincides with the end of overall sample therefore there are no sample points to evaluate predictions. 'forecast.flag' is set to FALSE - please change the sample sizes in case out-of-sample predictions are desired.")
    forecast.flag <- FALSE
  }
  
  # estimate: 
  
  params <- midas_estimate(est.y, est.x, NULL, est.xdate, polynomial, loss = loss, num.evals = num.evals, num.coef = num.coef, seed = seed, tau = tau, step_idx = step_idx, legendre_degree = legendre_degree, profiling = profiling)
  pred_in <- midas_forecast(params, est.x, NULL, polynomial, step_idx = step_idx, legendre_degree = legendre_degree)
  if(loss=="mse"){
    fit <- sqrt(mean((est.y - pred_in)^2))
  } 
  if (loss=="als"){
    r <- est.y - pred_in
    fit <- 1/length(r)*sum(r^2*abs( tau - as.numeric((r < 0)) ))
  }
  if (loss=="rq"){
    r <- est.y - pred_in
    fit <- 1/length(r)*sum(r*( tau - as.numeric((r < 0)) ))
  }
  pred <- NULL
  # schemes:
  if (forecast.flag){
    if (scheme=="fixed"){
      #params <- midas_estimate(est.y,est.x,est.lag.y,est.xdate,polynomial,loss = loss,num.evals = num.evals,num.coef = num.coef, seed = seed, tau = tau)
      pred <- midas_forecast(params,out.x,NULL,polynomial)
    } else {
      nroll <- nforecast
      if (nroll == 0)
        stop('Rolling window does not apply because there are no rolling periods. Decrease "EstEnd".')
      y.big <- c(est.y,out.y)
      x.big <- rbind(est.x,out.x)
      x.date.big <- rbind(est.xdate,out.xdate)
      y.date.big <- c(est.ydate,out.ydate)
      pred <- matrix(0,nrow=nroll,ncol=1)
      for (t in 1:nroll){
        if (scheme=="rolling"){
          est.y.roll <- y.big[t:(nobs-1+t)]
          est.x.roll <- x.big[t:(nobs-1+t),]
          est.date.roll <- x.date.big[t:(nobs-1+t),]             
        } else { 
          if (scheme=="expand"){
            est.y.roll <- y.big[1:(nobs-1+t)]
            est.x.roll <- x.big[1:(nobs-1+t),]
            est.x.date.roll <- x.date.big[1:(nobs-1+t),]    
          } else {
            stop('scheme should be set to either: fixed, rolling, expand. Check!')
          }
        }
        out.y.roll <- y.big[nobs+t]
        out.x.roll <- x.big[nobs+t,]
        out.y.dateroll <- y.date.big[nobs+t]
        if (t == 1){
          tmp_params <- midas_estimate(est.y.roll,est.x.roll,NULL,est.x.date.roll,polynomial,loss = loss,num.evals = num.evals,num.coef = num.coef, seed = seed, tau = tau, step_idx = step_idx, legendre_degree = legendre_degree)
        } else  {
          startx_all <- as.numeric(tmp_params)
          tmp_params <- midas_estimate(est.y.roll, est.x.roll, NULL, est.x.date.roll, polynomial, loss = loss, num.evals = num.evals, num.coef = num.coef, startx_all = startx_all, seed = seed, tau = tau, step_idx = step_idx, legendre_degree = legendre_degree)
        }
        tmp <- midas_forecast(tmp_params,t(as.matrix(out.x.roll)),NULL,polynomial)
        pred[t] <- tmp
        params <- rbind(params,tmp_params)
      }
    }
  }
  pred.obj <- NULL
  pred.obj$pred <- data.frame(dates=out.ydate,pred=pred)
  if(loss=="mse"){
    pred.obj$rmse <- sqrt(mean((out.y-pred)^2))
  }
  if(loss=="als"){
    r <- out.y-pred
    pred.obj$loss <- 1/length(r)*sum(r^2*abs( tau - as.numeric((r < 0)) ))
  }
  if(loss=="rq"){
    r <- out.y-pred
    pred.obj$loss <- 1/length(r)*sum(r*( tau - as.numeric((r < 0)) ))
  }
  est.obj <- NULL
  est.obj$params <- params
  est.obj$fit <- fit
  est.obj$pred <- pred_in
  est.obj$polynomial <- polynomial
  est.obj$info <- info
  return(list(est.obj=est.obj,pred.obj=pred.obj))
}


#' ARDL-MIDAS regression  
#' 
#' @description 
#'  Estimates and predicts using a single variate ARDL-MIDAS model.
#' @details 
#'  \ifelse{html}{\out{Several polynomial functional forms are available (input variable <code>polynomial</code>): <br><br> - <code>beta_w</code>: Beta polynomial</center> <br> - <code>rbeta_w</code>: restricted Beta polynomial <br> - <code>expalmon_w</code>: exponential Almon polynomial <br>  - <code>umidas_w</code>: unrestricted lags (U-MIDAS) <br> - <code>step_fun</code>: polynomial with step functions <br> - <code>legendre_w</code>: Legendre polynomials <br> <br> different forecasting schemes (input variable <code>scheme</code>): <br><br>  - <code>fixed</code>: fixed scheme <br> - <code>rolling</code>: rolling window <br> - <code>expand</code>: expanding window <br><br> and different loss functions (input variable <code>loss</code>): <br><br> - <code>mse</code>: least squares <br> - <code>als</code>: asymmetric least squares <br> - <code>rq</code>: quantile. <br> <br> The ARDL-MIDAS model is: <br> <center> y<sub>t</sub> =  &mu; + &Sigma;<sub>p</sub> &rho;<sub>p</sub> y<sub>t-p</sub> + &beta;  &Sigma;<sub>j</sub> &omega;<sub>j</sub>(&theta;)x<sub>t-1</sub> </center> <br> where &mu;, &beta;, &theta; and  &rho;<sub>p</sub> are model parameters, p is the number of low-frequency lags and &omega; is the weight function.}}{Several polynomial functional forms are available (input variable \code{polynomial}): \cr  - \code{beta_w}: Beta polynomial \cr - \code{rbeta_w}: restricted Beta polynomial \cr - \code{expalmon_w}: Exp Almon polynomial \cr - \code{umidas_w}: unrestricted lags (U-MIDAS) \cr - \code{step_fun}: polynomial with step functions  \cr - \code{legendre_w}: Legendre polynomials \cr different forecasting schemes (input variable \code{scheme}): \cr - \code{fixed}: fixed scheme \cr - \code{rolling}: rolling window \cr - \code{expand}: expanding window \cr  and different loss functions (input variable \code{loss}) \cr - \code{mse}: least squares \cr - \code{als}: asymmetric least squares \cr - \code{rq}: quantile.  \cr\cr The ARDL-MIDAS model is: \cr \deqn{y_t =  \mu + \sum_p \rho_p y_{t-p} + \beta \sum_j \omega_j(\theta)x_{t-1}} \cr where \eqn{\mu}, \eqn{\beta}, \eqn{\theta}, \eqn{\rho_p}  are model parameters, p is number of low-frequency and \eqn{\omega} is the weight function.}     
#' @usage 
#' midas_ardl(data.y, data.ydate, data.x, data.xdate, x.lag, 
#'   y.lag, est.start, est.end, horizon = 1, 
#'   polynomial = c("legendre_w", "beta_w", "rbeta_w", "expalmon_w", "umidas_w","step_fun"),
#'   scheme = c("fixed", "rolling", "expand"), loss = c("mse", "rq", "als"), ...)
#' @param data.y response variable series.
#' @param data.ydate response variable dates.
#' @param data.x predictor variable series.
#' @param data.xdate predictor variable dates.
#' @param x.lag number of high-frequency data lags.
#' @param y.lag number of low-frequency data lags.
#' @param est.start start date of the estimation sample (referenced with data.xdate).
#' @param est.end end date of the estimation sample (referenced with data.xdate).
#' @param horizon forecast horizon measured in predictor variable sampling frequency (default set 1 unit ahead).
#' @param polynomial MIDAS lag polynomial specification. Options are: Legendre (legendre_w), Beta desity (beta_w), restricted Beta density (rbeta_w), exponential Almon (expalmon_w), unrestricted MIDAS (umidas_w), step functions (step_fun).
#' @param scheme forecasting scheme. Options are: fixed scheme (fixed), rolling window scheme (rolling), expanding window scheme (expand).
#' @param loss loss function. Options are: mean squared error (mse), quantile (rq), asymmetric least squares (als).
#' @param ... optional parameters to feed into other functions. forecast.flag - TRUE/FALSE to compute out-of-sample predictions (default TRUE), disp.flag - TRUE/FALSE  to display MIDAS data structures (default FALSE), 
#'    num.evals - number of objective function evaluations using random starting parameter values in the case of non-linear MIDAS polynomial (default 1e4), 
#'    num.coef - number of best coefficients to use as starting values in nonlinear optimization (default 10),
#'    seed - value used in set.seed for randomly drawing initial starting values around OLS optimal solution,
#'    profiling - TRUE/FALSE to use MIDAS paramater profiling, coded only for rbeta_w polynomial, (default FALSE),
#'    step_idx - index of step function lags. If step_fun is used as a polynomial, it is best to specify this option too, otherwise, the program figures out the sampling frequency ratio and computes \code{step_idx} accordingly (message is displayed in this case),
#'    legendre_degree - degree of legendre polynomials. If legendre_w is used as a polynomial, it is best to specify this option too, otherwise, the value is set to 3 (message is displayed in this case),
#'    tau - quantile level for als and rq regressions. If eithr als or rq loss is used, this option must be specified, program stops if not, 
#' @return returns  midas_ardl list which contains parameter estimates, in- and out-of-sample statistics and predictions, and some information about the specification of the method used.
#' @export midas_ardl
midas_ardl <- function(data.y, data.ydate, data.x, data.xdate, x.lag, y.lag, est.start, est.end, horizon = 1,
                       polynomial = c("legendre_w","beta_w","rbeta_w","expalmon_w","umidas_w","step_fun"), 
                       scheme = c("fixed","rolling","expand"),
                       loss = c("mse","rq","als"),...){
  polynomial <- match.arg(polynomial)
  scheme <- match.arg(scheme)
  loss <- match.arg(loss)
  # deal with options
  options <- list(...)
  if(is.null(options$forecast.flag)){
    forecast.flag <- TRUE
  } else {
    forecast.flag <- options$forecast.flag
  }
  if(is.null(options$disp.flag)){
    disp.flag <- FALSE
  } else {
    disp.flag <- options$disp.flag
  }
  info <- tau <- legendre_degree <- step_idx <- num.evals <- num.coef <- seed <- profiling <- NULL
  if(polynomial%in%c("beta_w","rbeta_w","expalmon_w")){
    if(is.null(options$num.evals)) {
      num.evals <- 1e4
    } else {
      num.evals <- options$num.evals
    }
    info$num.evals <- num.evals
    if(is.null(options$num.coef)){
      num.coef <- 10
    } else {
      num.coef <- options$num.coef
    }
    info$num.coef <- num.coef
    if(is.null(options$seed)){
      seed <- 100
    } else {
      seed <- options$seed
    }
    info$seed <- seed
  }
  if(polynomial%in%c("rbeta_w")){
    if(is.null(options$profiling)){
      profiling <- FALSE
    } else {
      profiling <- options$profiling
    }
  }
  if(polynomial%in%"step_fun"){
    step_idx <- options$step_idx
    info$step_idx <- step_idx   
  }
  if(polynomial%in%"legendre_w"){
    legendre_degree <- options$legendre_degree
    info$legendre_degree <- legendre_degree   
  }
  if(loss=="als"){
    if(is.null(options$tau))
      stop("set quantile level tau in options")
    tau <- options$tau
    info$tau <- tau
  }
  if(loss=="rq"){
    if(is.null(options$tau))
      stop("set quantile level tau in options")
    tau <- options$tau
    info$tau <- tau
  }

  mf.data <- mixed_freq_data(data.y,data.ydate,data.x,data.xdate,x.lag,y.lag,horizon,est.start,est.end,disp.flag)
  # in-sample data
  est.y <- mf.data$est.y
  est.x <- mf.data$est.x
  est.lag.y <- mf.data$est.lag.y
  est.ydate <- mf.data$est.ydate
  est.xdate <- mf.data$est.xdate
  # out-of-sample data
  out.y <- mf.data$out.y
  out.x <- mf.data$out.x
  out.ydate <- mf.data$out.ydate
  out.xdate <- mf.data$out.xdate
  out.lag.y <- mf.data$out.lag.y
  
  info$est.y <- est.y
  info$est.x <- est.x
  info$out.y <- out.y
  info$x.lag <- x.lag
  info$y.lag <- y.lag
  info$out.x <- out.x
  info$out.lag.y <- out.lag.y
  
  nobs <- length(est.y)
  nforecast <- length(out.y)
  if(nforecast==0){
    if (forecast.flag)
      message("'forecast.flag' was set TRUE but the end of of the estimation sample coincides with the end of overall sample therefore there are no sample points to evaluate predictions. 'forecast.flag' is set to FALSE - please change the sample sizes in case out-of-sample predictions are desired.")
    forecast.flag <- FALSE
  }
  
  # estimate: 
  
  params <- midas_estimate(est.y, est.x, est.lag.y, est.xdate, polynomial, loss = loss, num.evals = num.evals, num.coef = num.coef, seed = seed, tau = tau, step_idx = step_idx, legendre_degree = legendre_degree, profiling = profiling)
  pred_in <- midas_forecast(params, est.x, est.lag.y, polynomial, step_idx = step_idx, legendre_degree = legendre_degree)
  if (loss=="mse"){
    fit <- sqrt(mean((est.y - pred_in)^2))
  } 
  if (loss=="als"){
    r <- est.y - pred_in
    fit <- 1/length(r)*sum(r^2*abs( tau - as.numeric((r < 0)) ))
  }
  if (loss=="rq"){
    r <- est.y - pred_in
    fit <- 1/length(r)*sum(r*( tau - as.numeric((r < 0)) ))
  }
  pred <- NULL
  # schemes:
  if (forecast.flag){
    if (scheme=="fixed"){
      #params <- midas_estimate(est.y,est.x,est.lag.y,est.xdate,polynomial,loss = loss,num.evals = num.evals,num.coef = num.coef, seed = seed, tau = tau)
      pred <- midas_forecast(params,out.x,out.lag.y,polynomial)
    } else {
      nroll <- nforecast
      if (nroll == 0)
        stop('Rolling window does not apply because there are no rolling periods. Decrease "EstEnd".')
      
      y.big <- c(est.y,out.y)
      x.big <- rbind(est.x,out.x)
      lag.y.big <-  rbind(est.lag.y,out.lag.y)
      x.date.big <- rbind(est.xdate,out.xdate)
      y.date.big <- c(est.ydate,out.ydate)
      pred <- matrix(0,nrow=nroll,ncol=1)
      for (t in 1:nroll){
        if (scheme=="rolling"){
          est.y.roll <- y.big[t:(nobs-1+t)]
          est.x.roll <- x.big[t:(nobs-1+t),]
          est.lag.y.roll <- lag.y.big[t:(nobs-1+t),]
          est.date.roll <- x.date.big[t:(nobs-1+t),]             
        } else { 
          if (scheme=="expand"){
            est.y.roll <- y.big[1:(nobs-1+t)]
            est.x.roll <- x.big[1:(nobs-1+t),]
            est.lag.y.roll <- lag.y.big[1:(nobs-1+t),]
            est.x.date.roll <- x.date.big[1:(nobs-1+t),]    
          } else {
            stop('scheme should be set to either: fixed, rolling, expand. Check!')
          }
        }
        out.y.roll <- y.big[nobs+t]
        out.x.roll <- x.big[nobs+t,]
        out.lag.y.roll <- lag.y.big[nobs+t,]
        out.y.dateroll <- y.date.big[nobs+t]
        if (t == 1){
          tmp_params <- midas_estimate(est.y.roll,est.x.roll,est.lag.y.roll,est.x.date.roll,polynomial,loss = loss,num.evals = num.evals,num.coef = num.coef, seed = seed, tau = tau, step_idx = step_idx, legendre_degree = legendre_degree)
        } else  {
          startx_all <- as.numeric(tmp_params)
          tmp_params <- midas_estimate(est.y.roll, est.x.roll, est.lag.y.roll, est.x.date.roll, polynomial, loss = loss, num.evals = num.evals, num.coef = num.coef, startx_all = startx_all, seed = seed, tau = tau, step_idx = step_idx, legendre_degree = legendre_degree)
        }
        tmp <- midas_forecast(tmp_params,t(as.matrix(out.x.roll)),t(as.matrix(out.lag.y.roll)),polynomial)
        pred[t] <- tmp
        params <- rbind(params,tmp_params)
      }
    }
  }
  pred.obj <- NULL
  pred.obj$pred <- data.frame(dates=out.ydate,pred=pred)
  if(loss=="mse"){
    pred.obj$rmse <- sqrt(mean((out.y-pred)^2))
  }
  if(loss=="als"){
    r <- out.y-pred
    pred.obj$loss <- 1/length(r)*sum(r^2*abs( tau - as.numeric((r < 0)) ))
  }
  if(loss=="rq"){
    r <- out.y-pred
    pred.obj$loss <- 1/length(r)*sum(r*( tau - as.numeric((r < 0)) ))
  }
  est.obj <- NULL
  est.obj$params <- params
  est.obj$pred <- pred_in
  est.obj$fit <- fit
  est.obj$polynomial <- polynomial
  est.obj$info <- info
  return(list(est.obj=est.obj,pred.obj=pred.obj))
}

#' MIDAS regression estimation function
#' 
#' @description 
#'  Estimates a single variate MIDAS model. 
#' @details 
#'  For specficiation details, \code{midas_dl} or \code{midas_ardl} function descriptions for more details.
#' @param est.y response variable. 
#' @param est.x predictor variable lags in MIDAS data format.
#' @param est.lag.y autoregressive lags of response variable (if NULL DL-MIDAS model is estimated).
#' @param est.xdate predictor variable lag dates in MIDAS data format.
#' @param polynomial MIDAS lag polynomial specification. 
#' @param loss loss function.
#' @param num.evals number of objective function evaluations using random starting parameter values.
#' @param num.coef number of best coefficients to use as starting values in nonlinear optimization.
#' @param startx_all starting values to feed into optimization algorithm.
#' @param seed value used in set.seed for randomly drawing initial starting values.
#' @param ... optional parameters to feed into other functions. 
#' @return returns estimates of coefficient vector for a desired model specification.
#' @export midas_estimate
#' @keywords internal
midas_estimate <- function(est.y,est.x,est.lag.y,est.xdate,polynomial,loss,num.evals,num.coef,startx_all = NULL,seed = NULL,...){
  options <- list(...)
  if(!is.null(est.lag.y))
    est.lag.y <- matrix(est.lag.y, nrow = length(est.y))
  if (loss == "mse")
    obj_fun <- mse_loss
  if (loss == "rq"){
    #stop("check step functions, legendre, umidas...")
    tau <- options$tau
    if(is.null(tau)){
      tau <- 0.95
      message("MIDAS quantile regression is used without specifying the quantile level. 'tau' was set set to default value (0.95). set 'tau' in options if this is not a good choice")
    }
    obj_fun <- rq_loss
  }
  if (loss == "als"){
    #stop("check step functions, legendre, umidas...")
    tau <- options$tau
    if(is.null(tau)){
      tau <- 0.95
      message("MIDAS quantile regression is used without specifying the quantile level. 'tau' was set set to default value (0.95). set 'tau' in options if this is not a good choice")
    }
    obj_fun <- als_loss
  }
  if (is.null(options$profiling)){
    profiling <- FALSE
  } else {
    profiling <- options$profiling
  }
  
  # nls-midas:
  if (polynomial=="beta_w" || polynomial=="expalmon_w") {
    if (polynomial=="beta_w")
      weight <- beta_w
    if (polynomial=="expalmon_w")
      weight <-  expalmon_w
    if(is.null(startx_all)){# generate initial param guess:
      if(!profiling)  
        startx_all <- get_start_midas(y = est.y, X = est.x, z = est.lag.y, loss = loss, weight, polynomial, num.evals=num.evals, num.coef=num.coef, seed=seed, tau=tau)
    }
    if(!profiling){
      p <- dim(est.lag.y)[2]
      if(is.null(p))
        p <- 0
      lower <- min(startx_all[,1])-.Machine$double.eps
      upper <- max(startx_all[,1])+.Machine$double.eps
      if(p !=0 ){
        for(n_lag in seq(p)){
          lower <- c(lower,min(startx_all[,n_lag+1]))-.Machine$double.eps
          upper <- c(upper,max(startx_all[,n_lag+1]))+.Machine$double.eps
        }
      }
      if (polynomial=="beta_w"){
        lower <- c(lower, min(startx_all[,p+2]),  1, 1)-.Machine$double.eps
        upper <- c(upper, max(startx_all[,p+2]), 70, 70)+.Machine$double.eps
      } else if (polynomial=="expalmon_w"){ 
        lower <- c(lower, min(startx_all[,p+2]), -0.1, -0.3)-.Machine$double.eps
        upper <- c(upper, max(startx_all[,p+2]), 0.3, -0.01)+.Machine$double.eps
      }
      est <- NULL
      for (j in seq(num.coef))
        est[[j]] <-  suppressWarnings(optimx::optimx(startx_all[j,], obj_fun, y = est.y, x = est.x, z = est.lag.y, weight = weight, tau = tau,lower = lower, upper = upper, scheme=c("L-BFGS-B")))
      
      estim <- NULL
      for (j in seq(num.coef))
        estim <- rbind(estim,est[[j]]$value)
      
      est <- est[[which.min(estim[,dim(estim)[2]])]]
      coef <- est[-((length(est)-7):(length(est)))]
      rownames(coef) <- ""
      ar_lags <- NULL
      if(!is.null(dim(est.lag.y))){
        for (i in 1:dim(est.lag.y)[2]) 
          ar_lags <- c(ar_lags, paste0("AR-",i))
      }
      colnames(coef) <- c("(Intercept)",ar_lags,"beta","k1","k2")
      
    } else {
      stop("MIDAS regression estimation with parameter profiling has only been implemented for restricted Beta polynomial. Please reset 'polynomial' to 'rbeta_w' and re-run.")
    }
  } else if (polynomial=="rbeta_w"){
    weight <- rbeta_w
    if (!profiling){
      if(is.null(startx_all)){# generate initial param guess:
        startx_all <- get_start_midas(y = est.y, X = est.x, z = est.lag.y, loss = loss, weight, polynomial, num.evals=num.evals, num.coef=num.coef,seed=seed, tau=tau)
      }
      p <- dim(est.lag.y)[2]
      if(is.null(p))
        p <- 0
      lower <- min(startx_all[,1])-.Machine$double.eps
      upper <- max(startx_all[,1])+.Machine$double.eps
      if( p!=0 ){
        for(n_lag in seq(p)){
          lower <- c(lower,min(startx_all[,n_lag+1]))-.Machine$double.eps
          upper <- c(upper,max(startx_all[,n_lag+1]))+.Machine$double.eps
        }
      }
      lower <- c(lower, min(startx_all[,p+2]),  1)-.Machine$double.eps
      upper <- c(upper, max(startx_all[,p+2]), 70)+.Machine$double.eps
      est <- NULL
      for (j in seq(num.coef))
        est[[j]] <-  suppressWarnings(optimx::optimx(startx_all[j,], obj_fun, y = est.y, x = est.x, z = est.lag.y, weight = weight, tau = tau, lower = lower, upper = upper, scheme=c("L-BFGS-B")))
      estim <- NULL
      for (j in seq(num.coef))
        estim <- rbind(estim,est[[j]]$value)
      
      est <- est[[which.min(estim[,dim(estim)[2]])]]
      coef <- est[-((length(est)-7):(length(est)))]
      rownames(coef) <- ""
      ar_lags <- NULL
      if(!is.null(dim(est.lag.y))){
        for (i in 1:dim(est.lag.y)[2]) 
          ar_lags <- c(ar_lags, paste0("AR-",i))
      }
      colnames(coef) <- c("(Intercept)",ar_lags,"beta","k")
    } else {
      max_iter <- options$prof_max_iter
      if (is.null(max_iter)){
        max_iter <- 50
      }
      if(loss=="mse"){
        which_loss <- 1
        tau <- 0
      } else if(loss=="als"){
        which_loss <- 2
      } else if(loss=="rq"){
        which_loss <- 3
      }
      if(!is.null(est.lag.y)){
        coef <- midasar_pr(as.vector(est.y),as.matrix(est.lag.y),as.matrix(est.x),as.double(1),as.double(tau),as.double(which_loss),as.double(max_iter))
      } else {
        coef <- midas_pr(as.vector(est.y),as.matrix(est.x),as.double(1),as.double(tau),as.double(which_loss),as.double(max_iter))
      }
      ar_lags <- NULL
      if(!is.null(est.lag.y)){
        for (i in 1:dim(est.lag.y)[2]) 
          ar_lags <- c(ar_lags, paste0("AR-",i))
      }
      coef <- data.frame(matrix(coef,nrow=1))
      rownames(coef) <- ""
      colnames(coef) <- c("(Intercept)",ar_lags,"beta","k")
    }#end of profiling
    
  } else if (polynomial=="step_fun"){
    step_idx <- options$step_idx
    hf_lags <- dim(est.x)[2]
    if (is.null(step_idx)) {
      if(hf_lags/3<=6){
        step_idx <- seq(3,hf_lags,by=3)
        message("step_fun is used without specifying step function indices. program guessed the covariate is monthly, and therefore used quarterly/monthly step function. set 'step_idx' in options if this is not a good choice")
      }
      if(hf_lags/22>=1){
        step_idx <- seq(1,5,22) 
        message("step_fun is used without specifying step function indices. program guessed the covariate is daily, and therefore used momthly/daily step function. set 'step_idx' in options if this is not a good choice")
      }
    }
    if(max(step_idx)>hf_lags){
      step_idx <- step_idx[-length(step_idx)]
      message("last 'step_idx' entry was removed as the index is larger than the total number of 'x' lags. Please reset 'step_idx' if this is not desired option")
    }
    if(max(step_idx)>hf_lags)
      stop("step index was wrongly set. Please set based on 'x.lag' and re-run.")
    # create covariates:
    x_cov <- x_names <- NULL
    lag_idx <- seq(1,hf_lags,by=1)
    n <- dim(est.x)[1]
    for (s_j in seq(length(step_idx))){
      x_cov <- cbind(x_cov, rowMeans(matrix(est.x[,which(lag_idx<=step_idx[s_j])],nrow=n)))
      x_names <- c(x_names, paste0("Step poly-", s_j))
    }
    if (loss=="mse"){
      if(is.null(est.lag.y)){
        est <- lm(est.y~x_cov)
      } else {
        est <- lm(est.y~est.lag.y+x_cov)
      }
      coef <- as.numeric(est$coefficients)
    } 
    if (loss=="als"){
      if(is.null(est.lag.y)){
        est <- fastals(est.y,x_cov,as.double(1),as.double(tau),as.double(1e3),as.double(1e-7))
      } else {
        est <- fastals(est.y,cbind(est.lag.y,x_cov),as.double(1),as.double(tau),as.double(1e3),as.double(1e-7)) 
      }
      coef <- as.numeric(est)
    }
    if (loss=="rq"){
      if(is.null(est.lag.y)){
        est <- quantreg::rq(est.y~x_cov,tau = tau)
      } else {
        est <- quantreg::rq(est.y~est.lag.y+x_cov, tau = tau)
      }
      coef <- as.numeric(est$coefficients)
    }
    
    ar_lags <- NULL
    if(!is.null(est.lag.y)){
      for (i in 1:dim(est.lag.y)[2]) 
        ar_lags <- c(ar_lags, paste0("AR-",i))
    }
    coef <- data.frame(matrix(coef,nrow=1))
    rownames(coef) <- ""
    colnames(coef) <- c("(Intercept)",ar_lags,x_names)
  } else if (polynomial=="umidas_w") {
    total_param <- dim(est.lag.y)[2]+dim(est.x)[2]
    if (total_param>length(est.y))
      stop("number of parameters exceed the sample size. change estimation start or end date")
    if (loss=="mse"){
      if(is.null(est.lag.y)){
        est <- lm(est.y~est.x)
      } else {
        est <- lm(est.y~est.lag.y+est.x)
      }
      coef <- as.numeric(est$coefficients)
    } 
    if (loss=="als"){
      if(is.null(est.lag.y)){
        est <- fastals(est.y,est.x,as.double(1),as.double(tau),as.double(1e3),as.double(1e-7))
      } else {
        est <-  fastals(est.y,cbind(est.lag.y,est.x),as.double(1),as.double(tau),as.double(1e3),as.double(1e-7))
      }
      coef <- as.numeric(est)
    }
    if (loss=="rq"){
      if(is.null(est.lag.y)){
        est <- quantreg::rq(est.y~est.lag.y+est.x,tau = tau)
      } else {
        est <- quantreg::rq(est.y~est.lag.y+est.x,tau = tau)
      }
      coef <- as.numeric(est$coefficients)
    }
    ar_lags <- NULL
    if(!is.null(est.lag.y)){
      for (i in 1:dim(est.lag.y)[2]) 
        ar_lags <- c(ar_lags, paste0("AR-",i))
    }
    x_names <- NULL
    for (i in 1:dim(est.x)[2]) 
      x_names <- c(x_names, paste0("UMIDAS-",i))
    
    coef <- data.frame(matrix(coef,nrow=1))
    rownames(coef) <- ""
    colnames(coef) <- c("(Intercept)",ar_lags,x_names)
    
  } else if (polynomial=="legendre_w") {
    legendre_degree <- options$legendre_degree
    if (is.null(legendre_degree)) {
      legendre_degree <- 3
      message("Legendre polynomials are used without specifying the degree. default was set (degree=3).  set 'legendre_degree' in options if this is not a good choice")
    }
    w <- lb(legendre_degree, a = 0, b = 1, jmax = dim(est.x)[2])
    xw <- est.x%*%w
    
    if (loss=="mse"){
      if(is.null(est.lag.y)){
        est <- lm(est.y~xw)
      } else {
        est <- lm(est.y~est.lag.y+xw)
      }
      coef <- as.numeric(est$coefficients)
    } 
    if (loss=="als"){
      if(is.null(est.lag.y)){
        est <- fastals(est.y,xw,as.double(1),as.double(tau),as.double(1e3),as.double(1e-7))
      } else {
        est <- fastals(est.y,cbind(est.lag.y,xw),as.double(1),as.double(tau),as.double(1e3),as.double(1e-7))
      }
      coef <- as.numeric(est)
    }
    if (loss=="rq"){
      if(is.null(est.lag.y)){
        est <- quantreg::rq(est.y~xw,tau = tau)
      } else {
        est <- quantreg::rq(est.y~est.lag.y+xw,tau = tau)
      }
      coef <- as.numeric(est$coefficients)
    }
    ar_lags <- NULL
    if(!is.null(est.lag.y)){
      for (i in 1:dim(est.lag.y)[2]) 
        ar_lags <- c(ar_lags, paste0("AR-",i))
    }
    x_names <- NULL
    for (i in 1:(legendre_degree+1))
      x_names <- c(x_names, paste0("Legendre poly-",i-1))
    
    coef <- data.frame(matrix(coef,nrow=1))
    rownames(coef) <- ""
    colnames(coef) <- c("(Intercept)",ar_lags,x_names)
    
  }
  return(coef)
}

#' MIDAS regression prediction function 
#' 
#' @description 
#'  Predicts from single variate MIDAS model estimated based on specification supplied to \code{midas_estimate} function.
#' @details 
#'  For specficiation details, see \code{midas_dl} or \code{midas_ardl} function descriptions for more details.
#' @param params parameter vector from \code{midas_estimate}.
#' @param x out-of-sample predictor variable data.
#' @param ylag out-of-sample lagged depedent variable data.
#' @param polynomial polynomial specification.
#' @param ... optional parameters to feed into other functions. 
#' step_idx - index for step function polynomial specification (warning: if left unspecified, the program computes index the same way as in the estimation function), 
#' legendre_degree - the degree of Legendre polynomials (warning: if left unspecified, the program sets it to 3, the same way as in the estimation function).
#' @return returns prediction value.
#' @export midas_forecast
#' @keywords internal
midas_forecast <- function(params,x,ylag,polynomial,...){
  options <- list(...)
  n <- dim(x)[1]
  d <- dim(x)[2]
  iota <- rep(1,times=n)
  c <- arpred <- midaspred <- 0 
  c <- iota*as.numeric(params[1])
  params_c <- params[-1]
  
  if(!is.null(ylag)){
    ylag <- matrix(ylag, nrow = n)
    ard <- dim(ylag)[2]
    param_ar <- as.numeric(params_c[1:ard])
    param_midas <- as.numeric(params_c[-c(1:ard)])
    arpred <- ylag%*%param_ar
  } else {
    param_midas <- as.numeric(params_c)
  }
  if(polynomial%in%c("beta_w","rbeta_w","expalmon_w")) {
    # nonlinear MIDAS schemes
    if (polynomial=="beta_w")
      weight <- beta_w
    if (polynomial=="expalmon_w")
      weight <-  expalmon_w
    if (polynomial=="rbeta_w")
      weight <- rbeta_w
    midaspred <- param_midas[1]*x%*%weight(param_midas[-1], d)
  } else {
    # linear in parameters MIDAS schemes
    if (polynomial=="step_fun"){
      step_idx <- options$step_idx
      hf_lags <- dim(x)[2]
      if (is.null(step_idx)) {
        if(hf_lags/3<=6){
          step_idx <- seq(3,hf_lags,by=3)
        }
        if(hf_lags/22>=1){
          step_idx <- c(1,5,22)
        }
      }
      x_m <- NULL
      lag_idx <- seq(1,hf_lags,by=1)
      n <- dim(x)[1]
      for (s_j in seq(length(step_idx))){
        x_m <- cbind(x_m, rowMeans(matrix(x[,which(lag_idx<=step_idx[s_j])],nrow=n)))
      }
      midaspred <- x_m%*%param_midas
    }
    if (polynomial=="legendre_w"){
      legendre_degree <- options$legendre_degree
      if (is.null(legendre_degree)) {
        legendre_degree <- 3
      }
      w <- lb(legendre_degree, a = 0, b = 1, jmax = dim(x)[2])
      xw <- x%*%w
      midaspred <- xw%*%param_midas
    }
    if (polynomial=="umidas_w"){
      midaspred <-x%*%param_midas
    }
  }
  if (is.null(arpred))
    arpred <- 0
  
  pred <- c + arpred + midaspred
  return(pred)
}

#' MIDAS regression function for initial values
#' 
#' @description 
#'  Computes initial values for different MIDAS regression model specifications.
#' @details 
#'  For a given loss function and MIDAS weight function specification, computes \code{num.evals} number of random initial values and evaluates the objective function at those parameters. 
#'  Function retains \code{num.coef} of best in terms of fit initial starting values, which are then feed into optimization algorithms.
#' @param params parameter vector from \code{midas_estimate}.
#' @param X predictor variable lags in MIDAS data format.
#' @param z autoregressive lags of response variable (default NULL, DL-MIDAS model is computed).
#' @param loss polynomial specification.
#' @param weight MIDAS weight function (depends on \code{polynomial}).
#' @param polynomial MIDAS polynomial specification.
#' @param num.evals number of objective function evaluations using random starting parameter values.
#' @param num.coef number of best coefficients to use as starting values in nonlinear optimization.
#' @param seed value used in set.seed for randomly drawing initial starting values.
#' @param tau quantile level for als and rq regressions. 
#' @return returns \code{num.coef} number of initial parameter values.
#' @export get_start_midas
#' @keywords internal
get_start_midas <- function(y, X, z = NULL, loss = c("mse", "als", "rq"), weight, polynomial, num.evals=1000, num.coef=10, seed=NULL, tau=tau) {
  if (is.null(seed))
    seed <- 100
  if (polynomial=="rbeta_w")
    nw <- 1
  if (polynomial=="beta_w" ||  polynomial=="expalmon_w")
    nw <- 2
  
  loss <- match.arg(loss)
  if (loss=="als"){
    if(is.null(tau)){
      tau <- 0.95
      message("MIDAS als regression is used without specifying the quantile level. 'tau' was set set to default value (0.95). set 'tau' in options if this is not a good choice")
    }
  }
  if (loss=="rq"){
    if(is.null(tau)){
      tau <- 0.95
      message("MIDAS quantile regression is used without specifying the quantile level. 'tau' was set set to default value (0.95). set 'tau' in options if this is not a good choice")
    }
  }
  set.seed(seed)
  d <- ncol(X)
  model <- na.omit(cbind(y, X, z))
  y <- model[, 1]
  XX <- model[, -1]
  n <- nrow(model)
  if (is.null(z)) {
    fit_xb <- function(param_eval,XX,nw) {
      n <- dim(XX)[1]
      iota <- rep(1,times=n)
      xb <- param_eval[1]*iota + param_eval[2]*XX%*%weight(param_eval[-c(1:2)], d)
      xb
    }
    p.eval <- matrix(NA,ncol=1+1+nw,nrow=num.evals)
    
    if(loss=="mse"){
      # get slope via OLS:
      fit_ls <- as.numeric(lm(y~rowMeans(X))$coef)
    } 
    if(loss=="als"){
      fit_ls <- as.numeric(fastals(y,matrix(rowMeans(X),nrow=dim(X)[1],ncol=1),as.double(1),as.double(tau), as.double(1e3), as.double(1e-6)) )  
    }
    if(loss=="rq"){
      # get slope via rq:
      fit_ls <- as.numeric(quantreg::rq(y~rowMeans(X), tau = tau)$coef)
    }
    # fill in with intercept:
    p.eval[,1] <- rep(fit_ls[1],times=num.evals)+fit_ls[1]*runif(num.evals,min=-0.1,max=0.1)
    # fill in slope:
    p.eval[,2] <- rep(fit_ls[2],times=num.evals)+fit_ls[2]*runif(num.evals,min=-0.1,max=0.1)
    if (polynomial=="rbeta_w")
      p.eval[,3] <- rep(5,times=num.evals) + runif(num.evals,min=-4,max=14)
    if (polynomial=="beta_w")
      p.eval[,c(3:4)] <- matrix(rep(5,times=num.evals*2) + runif(num.evals*2,min=-4,max=14),nrow=num.evals,ncol=2)
    if (polynomial=="expalmon_w"){
      p.eval[,3] <- rep(0,times=num.evals) + runif(num.evals,min=-0.1,max=0.3)
      p.eval[,4] <- rep(0,times=num.evals) + runif(num.evals,min=-0.3,max=-0.01)
    }
    
  } else {
    z <- as.matrix(z,nrow=n)
    p <- ncol(z)
    fit_xb <- function(param_eval,XX,nw) {
      n <- dim(XX)[1]
      iota <- rep(1,times=n)
      idx_m <- (length(param_eval)-nw):length(param_eval)
      param_midas <- param_eval[idx_m]
      param_ar <- param_eval[-idx_m]
      lar <- length(param_ar[-1])
      X_ar <- matrix(XX[,1:lar],nrow=n)
      X_m <- matrix(XX[,-c(1:lar)],nrow=n)
      xb <- param_ar[1]*iota + X_ar%*%param_ar[-1] + param_midas[1]*X_m%*%weight(param_midas[-1], d)
      xb
    }
    p.eval <- matrix(NA,ncol=1+p+1+nw,nrow=num.evals)
    if(loss=="mse"){
      # get slope via OLS:
      fit_ls <- as.numeric(lm(y~z+rowMeans(X))$coef)
    } 
    if(loss=="als"){
      fit_ls <- as.numeric(fastals(y,cbind(z,matrix(rowMeans(X),nrow=dim(X)[1],ncol=1)),as.double(1),as.double(tau), as.double(1e3), as.double(1e-6)) )  
    }
    if(loss=="rq"){
      # get slope via rq:
      fit_ls <- as.numeric(quantreg::rq(y~z+rowMeans(X))$coef,tau = tau)
    }
    # fill in with intercept:
    p.eval[,1] <- rep(fit_ls[1],times=num.evals)+fit_ls[1]*runif(num.evals,min=-0.1,max=0.1)
    # fill in lags:
    for (n_lag in seq(p))
      p.eval[,n_lag+1] <- rep(fit_ls[n_lag+1],times=num.evals)+fit_ls[n_lag+1]*runif(num.evals,min=-0.1,max=0.1)
    # fill in slope:
    p.eval[,p+2] <- rep(fit_ls[p+2],times=num.evals)+fit_ls[p+2]*runif(num.evals,min=-0.1,max=0.1)
    # fill in midas weights parameters
    if (polynomial=="rbeta_w")
      p.eval[,p+3] <- rep(5,times=num.evals) + runif(num.evals,min=-4,max=14)
    if (polynomial=="beta_w")
      p.eval[,c((p+3):(p+4))] <- matrix(rep(5,times=num.evals*2) + runif(num.evals*2,min=-4,max=14),nrow=num.evals,ncol=2)
    if (polynomial=="expalmon_w"){
      p.eval[,p+3] <- rep(0,times=num.evals) + runif(num.evals,min=-0.1,max=0.3) 
      p.eval[,p+4] <- rep(0,times=num.evals) + runif(num.evals,min=-0.3,max=-0.01)
    }
    
  }
  if (loss=="mse"){
    fn0 <- function(param_eval,XX,nw) {
      1/n*sum((y - fit_xb(param_eval,XX,nw))^2)
    }
  }
  if (loss=="als"){
    fn0 <- function(param_eval,XX,nw) {
      r  <- y - fit_xb(param_eval,XX,nw)
      1/n*sum(r^2*abs( tau - as.numeric((r < 0)) ))
    }
  }
  if (loss=="rq"){
    fn0 <- function(param_eval,XX,nw) {
      r  <- y - fit_xb(param_eval,XX,nw)
      1/n*sum(r*( tau - as.numeric((r < 0)) ))
    }
  }
  f.eval <- apply(p.eval,1,fn0, XX, nw)
  all <- cbind(f.eval,p.eval)
  all <- all[order(all[,1]),]
  coefs <- all[1:num.coef,-1]
  return(coefs)
}

#' MIDAS weights plot function
#' 
#' @description 
#'  Based on specification in \code{obj}, plots a basic R figure of estimated MIDAS weights.
#' @details 
#'  MIDAS regression pecifcation is picked up from obj, see \code{midas_dl} or \code{midas_ardl} function descriptions for more details. 
#' @param obj midas_ardl or midas_dl object with parameter estimates and model specfication inputs.
#' @return returns R figure of estimated MIDAS weights.
#' @export plot_weights
plot_weights <- function(obj){
  polynomial <- obj$polynomial
  params <- obj$params
  dlag <- obj$info$x.lag
  if(polynomial%in%c("beta_w","rbeta_w","expalmon_w")){
    if (polynomial%in%"rbeta_w"){
      d <- 2
      weight <- rbeta_w
    } else {
      d <- 3
      if (polynomial%in%"beta_w")
        weight <- beta_w
      if (polynomial%in%"expalmon_w")
        weight <- expalmon_w
    }
    params_midas <- as.numeric(params[(length(params)-d+1):length(params)])
    w <- params_midas[1]*weight(params_midas[-1],dlag)
  } else {
    if(polynomial%in%"umidas_w"){
      w <- as.numeric(params[(length(params)-dlag+1):length(params)])
    }
    
    if (polynomial%in%"step_fun"){
      step_idx <- obj$info$step_idx
      if (is.null(step_idx))
        stop("information about estimated step function scheme weights is missing")
      
      w <- numeric(max(step_idx))
      d <- length(step_idx)
      params_step <- as.numeric(params[(length(params)-d+1):length(params)])
      lag_idx <- seq(1,dlag,by=1)
      for (j in seq(length(step_idx),1,by=-1))
        w[which(lag_idx<=step_idx[j])] <- rep(params_step[j],times=step_idx[j])
    }
    
    if (polynomial%in%"legendre_w"){
      legendre_degree <- obj$info$legendre_degree
      if (is.null(legendre_degree))
        stop("information about estimated Legender polynomial scheme weights is missing")
      
      wlb <- lb(legendre_degree, jmax = dlag)
      d <- legendre_degree+1
      params_lb <- as.numeric(params[(length(params)-d+1):length(params)])
      w <- as.numeric(wlb%*%params_lb)
    }
  }
  grid <- 1:dlag
  title <- paste0("Estimated weights for ",polynomial," scheme")
  plot(grid, w, main = title , type = 'l', xlab = 'Lag', ylab = ' Weight')
}

