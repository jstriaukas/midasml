#' MIDAS ML regression prediction function
#' 
#' @description 
#'  Predicts from a high-dimensional MIDAS model.  
#' @details 
<<<<<<< HEAD
#'  \ifelse{html}{\out{Based on desired computation of the tuning parameter &lambda; (<code>ic/cv</code>), the MIDAS ML regression is fit and predictions are computed based on the chosen scheme (<code>rolling/expanding</code>). The regression function that is fit is <br><br> <center> RSS(&alpha;, &beta;) + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where X (<code>x_in</code> and <code>x_out</code>) contains autoregressive lags and MIDAS covariates of <code>sort_midasml</code> class and RSS(&alpha;, &beta;) is the mse loss function. The group index vector should be used from data sorting function <code>qtarget.sort_midasml</code>. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)||&beta;||<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{Based on desired computation of the tuning parameter \eqn{\lambda} (\code{ic/cv}), the MIDAS ML regression is fit and predictions are computed based on the chosen scheme (\code{rolling/expanding}). The regression function that is fit is \cr \cr  \deqn{||y - \alpha \iota - X \beta||^2_T + 2\lambda  \Omega_\gamma(\beta)}, \cr where \eqn{X} (\code{x_in} and \code{x_out}) contains autoregressive lags and MIDAS covariates of \code{sort_midasml} class. The group index vector should be used from data sorting function \code{qtarget.sort_midasml}. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)||\beta||_2,1,} a convex combination of LASSO and group LASSO penalty functions.}
=======
<<<<<<< HEAD
#'  \ifelse{html}{\out{Based on desired computation of the tuning parameter &lambda; (<code>ic/cv</code>), the MIDAS ML regression is fit and predictions are computed based on the chosen scheme (<code>rolling/expanding</code>). The regression function that is fit is <br><br> <center> RSS(&alpha;, &beta;) + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where X (<code>x_in</code> and <code>x_out</code>) contains autoregressive lags and MIDAS covariates of <code>sort_midasml</code> class and RSS(&alpha;, &beta;) is the mse loss function. The group index vector should be used from data sorting function <code>qtarget.sort_midasml</code>. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)||&beta;||<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{Based on desired computation of the tuning parameter \eqn{\lambda} (\code{ic/cv}), the MIDAS ML regression is fit and predictions are computed based on the chosen scheme (\code{rolling/expanding}). The regression function that is fit is \cr \cr  \deqn{||y - \alpha \iota - X \beta||^2_T + 2\lambda  \Omega_\gamma(\beta)}, \cr where \eqn{X} (\code{x_in} and \code{x_out}) contains autoregressive lags and MIDAS covariates of \code{sort_midasml} class. The group index vector should be used from data sorting function \code{qtarget.sort_midasml}. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)||\beta||_2,1,} a convex combination of LASSO and group LASSO penalty functions.}
=======
#'  \ifelse{html}{\out{Based on desired computation of the tuning parameter &lambda; (<code>ic/cv</code>), the MIDAS ML regression is fit and predictions are computed based on the chosen scheme (<code>rolling/expanding</code>). The regression function that is fit is <br><br> <center> RSS(&alpha;, &beta;) + 2&lambda;  &Omega;<sub>&gamma;</sub>(&beta;), </center> <br> where X contains autoregressive lags and MIDAS covariates of <code>sort_midasml</code> class and RSS(&alpha;, &beta;) is the mse loss function. The group index vector should be used from data sorting function <code>qtarget.sort_midasml</code>. The penalty function &Omega;<sub>&gamma;</sub>(.) is applied on &beta; coefficients and is <br> <br> <center> &Omega;<sub>&gamma;</sub>(&beta;) = &gamma; |&beta;|<sub>1</sub> + (1-&gamma;)||&beta;||<sub>2,1</sub>, </center> <br> a convex combination of LASSO and group LASSO penalty functions.}}{Based on desired computation of the tuning parameter \eqn{\lambda} (\code{ic/cv}), the MIDAS ML regression is fit and predictions are computed based on the chosen scheme (\code{rolling/expanding}). The regression function that is fit is \cr \cr  \deqn{||y - \alpha \iota - X \beta||^2_T + 2\lambda  \Omega_\gamma(\beta)}, \cr where \eqn{X} contains autoregressive lags and MIDAS covariates of \code{sort_midasml} class. The group index vector should be used from data sorting function \code{qtarget.sort_midasml}. The penalty function \eqn{\Omega_\gamma(.)} is applied on \eqn{\beta} coefficients and is \deqn{\Omega_\gamma(\beta) = \gamma |\beta|_1 + (1-\gamma)||\beta||_2,1,} a convex combination of LASSO and group LASSO penalty functions.}
>>>>>>> 275e165e23cdf61c5be22beec7badf7ddf049b09
>>>>>>> 194273b8ffc16c6fb69144fc04b4d1799adb8048
#' @usage 
#' midasml_forecast(y_in, y_out, x_in, x_out, group_index,
#'   gamma_w, y_out_dates, scheme, verbose = FALSE, ...)
#' @param y_in response variable (in-sample). 
#' @param y_out response variable (out-of-sample).
#' @param x_in predictor variables (in-sample).
#' @param x_out predictor variables (out-of-sample).
#' @param group_index  group membership of each covariate.
#' @param gamma_w sg-LASSO mixing parameter. \ifelse{html}{\out{<code>gamma_w = 1</code>}}{\code{gamma_w = 1}} is LASSO and \ifelse{html}{\out{<code>gamma_w = 0</code>}}{\code{gamma_w = 0}} group LASSO.
#' @param y_out_dates out-of-sample dates.
#' @param scheme prediction scheme. Choices are: \ifelse{html}{\out{<code>expand</code>}}{\code{expand}} - expanding window scheme, \ifelse{html}{\out{<code>rolling</code>}}{\code{rolling}} - rolling window scheme.
#' @param verbose flag to print information. 
#' @param ... optional parameters to feed into \ifelse{html}{\out{<code>reg_sgl</code>}}{\code{reg_sgl}}. 
#' @return out-of-sample predictions.
#' @export midasml_forecast
midasml_forecast <- function(y_in, y_out, x_in, x_out, group_index, gamma_w, y_out_dates, scheme, verbose=FALSE, ...){
  y_hat <- numeric(length(y_out))
  if(scheme=="expand"){# expanding window scheme:
    for (i in 1:length(y_out)) {
      if(verbose)
        message(paste0("\n", "estimating for the quarter: ", y_out_dates[i]))
      if(i==1){
        fit <- reg_sgl(X = x_in, y = y_in, index = group_index, gamma_w = gamma_w, verbose = verbose, ...)
      } else {
        fit <- reg_sgl(X = rbind(x_in,x_out[i-1,]), y = c(y_in, y_out[i]), 
                              index = group_index, gamma_w = gamma_w, verbose = verbose, ...)
      }
      y_hat[i] <- predict.reg_sgl(fit, newX = as.matrix(x_out[i,]))$pred
    }
  }
  if(scheme=="rolling"){
    for (i in 1:length(y_out)) {
      if(verbose)
        message(paste0("estimating for the quarter: ", y_out_dates[i]))
      if(i==1){
        fit <- reg_sgl(X = x_in, y = y_in, index = group_index, gamma_w = gamma_w, verbose = verbose, ...)
      } else {
        fit <- reg_sgl(X = rbind(x_in[-c(1:(i-1)),],x_out[i-1,]), y = c(y_in[-c(1:(i-1)),], y_out[i]), 
                       index = group_index, gamma_w = gamma_w, verbose = verbose, ...)
      }
      y_hat[i] <- predict.reg_sgl(fit, newX = as.matrix(x_out[i,]))$pred
    }
  }
  return(y_hat)
}
#' High-dimensional mixed frequency data sort function
#' 
#' @description 
#'  Sorts high-dimensional mixed frequency data for quarterly target variable.  
#' @details 
#'  \ifelse{html}{\out{Sorts high-dimensional mixed frequency data for quarterly target variable read to be inputed into <code>midasml_forecast</code> function to produce nowcasts & forecasts.}}{Sorts high-dimensional mixed frequency data for quarterly target variable read to be inputed into \code{midasml_forecast} function to produce nowcasts & forecasts.}
#' @usage
#' qtarget.sort_midasml(y.data, x.macro.data = NULL, x.real.time = NULL,
#'   x.quarterly_group = NULL, x.lag = NULL, legendre_degree, horizon, 
#'   macro_delay = 1, est.start, est.end, standardize = TRUE, group_ar_lags = FALSE, 
#'   disp.flag = TRUE)
#' @param y.data response variable data. 
#' @param x.macro.data macro data which is not real-time, i.e. is used with publication delay defined in \ifelse{html}{\out{<code>macro_delay</code>}}{\code{macro_delay}}.
#' @param x.real.time real-time data.
#' @param x.quarterly_group quarterly data currently taken as real-time data.
#' @param x.lag single value or vector of size of the total number of variables defining the number of lags for each high-frequency variable in \ifelse{html}{\out{<code>x.macro.data, x.real.time</code>}}{\code{x.macro.data, x.real.time}}.
#' @param legendre_degree single value or vector of size of the total number of variables defining the polynomial degree for each each high-frequency variable in \ifelse{html}{\out{<code>x.macro.data, x.real.time</code>}}{\code{x.macro.data, x.real.time}}.
#' @param horizon forecast horizon relative to \ifelse{html}{\out{<code>y.data</code>}}{\code{y.data}} date column in high-frequency time units.
#' @param macro_delay number of months that macro series in \ifelse{html}{\out{<code>x.macro.data</code>}}{\code{x.macro.data}} are delayed.
#' @param est.start estimation start date, taken as the first ... .
#' @param est.end estimation end date, taken as the last ... . Remainig data after this date is dropped to out-of-sample evaluation data. 
#' @param standardize TRUE/FALSE to standardize high-frequneyc covariates in high-frequency units.
#' @param group_ar_lags TRUE/FALSE to group AR lags.
#' @param disp.flag display flag to indicate whether or not to display obtained MIDAS data structure in console.
#' @return MIDAS covariates and group memberships based on desired specification.
#' @export qtarget.sort_midasml
qtarget.sort_midasml <- function(y.data, x.macro.data = NULL, x.real.time = NULL, x.quarterly_group = NULL, x.lag = NULL, legendre_degree, horizon, macro_delay = 1, est.start, est.end, standardize = TRUE, group_ar_lags = FALSE, disp.flag = TRUE){
  dim_macro <- dim_real.time <- dim_quarterly <- 0
  if(!is.null(x.macro.data))
    dim_macro <- dim(x.macro.data)[2]-1
  if(!is.null(x.real.time))
    dim_real.time <- dim(x.real.time)[2]-1
  
  dim_x <- dim_macro + dim_real.time + 1
  if(is.null(x.lag))
    stop("x.lag variable, which defines the lag structure of each covariate, must be specified.")
  if(length(x.lag)!=1 && length(x.lag)!=(dim_macro + dim_real.time))
    stop(paste0("x.lag variable length must be either of size 1 (the same lag structure for x.macro.data and/or x.real.time) or must be of the length size equal to the total number of high-frequency covariates: ",dim_macro + dim_real.time,"."))
  if(length(x.lag)==1)
    x.lag <- rep(x.lag,times=(dim_macro + dim_real.time))
  
  if(length(legendre_degree)!=1 && length(legendre_degree)!=dim_x){
    message(paste0("Legendre polynomial degree must be specified the same for all covariates (one number) or a seperate value for each (the length size equals to the total number of covariates). the length of legendre_degree: ", length(legendre_degree), ", number of covaraites that are inputed: ", dim_x,". Legendre degree is set to (for all covariates): ",  legendre_degree[1], " - the first entry in legendre_degree"))
    legendre_degree <- legendre_degree[1]
  }
  if (length(legendre_degree)==1)
    legendre_degree <- rep(legendre_degree,times=dim_x)
  
  if(macro_delay!=1)
    message("typically macro series is published with 1 month lag. please check if your inputed macro data has different publication lag. the program does not stop, you need to re-run by re-setting macro_delay input.")
  
  # storage
  x_str_out <- x_str <- NULL
  x_unstr_out <- x_unstr <- NULL
  x_average_out <- x_average <- NULL
  group_index <- group_index_un <- group_index_av <- 0
  ref_in <- ref_out <- NULL
  # computing the macro data if inputed
  if(!is.null(x.macro.data)){
    x.macro_date <- x.macro.data$DATE
    x.macro_data <- x.macro.data[,-1]
    data.refdate <- y.data$DATE
    lubridate::month(data.refdate) <- lubridate::month(data.refdate)-macro_delay
    lubridate::month(est.start) <- lubridate::month(est.start)-macro_delay
    lubridate::month(est.end) <- lubridate::month(est.end)-macro_delay
    for (j_macro in seq(dim_macro)){
      if(standardize){
        j_data <- scale(x.macro_data[,j_macro], center = TRUE, scale = TRUE)
      } else {
        j_data <- scale(x.macro_data[,j_macro], center = TRUE, scale = TRUE)
      }
      # get MIDAS structure:
      tmp <- mixed_freq_data_single(data.refdate = data.refdate, data.x = j_data, data.xdate = x.macro_date,
                                    x.lag[j_macro], horizon, est.start, est.end, disp.flag = disp.flag)
      if(j_macro==1){
        ref_in <- tmp$est.refdate
        ref_out <- tmp$out.refdate
        lubridate::month(ref_in) <- lubridate::month(ref_in)+macro_delay
        lubridate::month(ref_out) <- lubridate::month(ref_out)+macro_delay
      }
      # get Legendre weights:
      tmp_w <- lb(legendre_degree[j_macro],a=0,b=1,jmax=x.lag[j_macro])
      # aggregate in-sample:
      x_str <- cbind(x_str, tmp$est.x%*%tmp_w)
      x_str_out <- cbind(x_str_out, tmp$out.x%*%tmp_w)
      # store unrestricted case:
      x_unstr <- cbind(x_unstr, tmp$est.x)
      x_unstr_out <- cbind(x_unstr_out, tmp$out.x)
      # store averages:
      x_average <- cbind(x_average, rowMeans(tmp$est.x))
      x_average_out <- cbind(x_average_out, rowMeans(tmp$out.x))
      # get group indices
      group_index <- c(group_index, rep(max(group_index)+1,times=legendre_degree[j_macro]+1))
      group_index_un <- c(group_index_un, rep(max(group_index_un)+1,times=x.lag[j_macro]))
      group_index_av <- c(group_index_av, max(group_index_av)+1)
    }
  }
  # computing the real-time data if inputed
  if(!is.null(x.real.time)){
    x.real.time_date <- x.real.time$DATE
    x.real.time_data <- x.real.time[,-1]
    data.refdate <- y.data$DATE
    for (j_real.time in seq(dim_real.time)){
      if(standardize){
        j_data <- scale(x.real.time_data[,j_real.time], center = TRUE, scale = TRUE)
      } else {
        j_data <- scale(x.real.time_data[,j_real.time], center = TRUE, scale = TRUE)
      }
      # get MIDAS structure:
      tmp <- mixed_freq_data_single(data.refdate = data.refdate, data.x = j_data, data.xdate = x.real.time_date,
                                    x.lag[j_real.time+dim_macro], horizon, est.start, est.end, disp.flag = disp.flag)
      
      if(j_real.time==1){
        ref_in <- tmp$est.refdate
        ref_out <- tmp$out.refdate
      }
      # get Legendre weights:
      tmp_w <- lb(legendre_degree[j_real.time+dim_macro],a=0,b=1,jmax=x.lag[j_real.time+dim_macro])
      # store aggregated case:
      x_str <- cbind(x_str, tmp$est.x%*%tmp_w)
      x_str_out <- cbind(x_str_out, tmp$out.x%*%tmp_w)
      # store unrestricted case:
      x_unstr <- cbind(x_unstr, tmp$est.x)
      x_unstr_out <- cbind(x_unstr_out, tmp$out.x)
      # store averages:
      x_average <- cbind(x_average, rowMeans(tmp$est.x))
      x_average_out <- cbind(x_average_out, rowMeans(tmp$out.x))
      # store group indices:
      group_index <- c(group_index, rep(max(group_index)+1,times=legendre_degree[j_real.time+dim_macro]+1))
      group_index_un <- c(group_index_un, rep(max(group_index_un)+1,times=x.lag[j_real.time+dim_macro]))
      group_index_av <- c(group_index_av, max(group_index_av)+1)
    }
  }  
  
  if(is.null(ref_in) || is.null(ref_out))
    stop("ref dates were not computed. likely that both macro and real time datasets where not inputed. at least one dataset must be inputed.")
  
  # computing the quarterly data if inputed
  if(!is.null(x.quarterly_group)){
    # sort quarterly group of covariates and lags
    x.quarterly_group_in <- x.quarterly_group[as.Date(x.quarterly_group$DATE)%in%ref_in,]
    x.quarterly_group_out <- x.quarterly_group[as.Date(x.quarterly_group$DATE)%in%ref_out,]
    
    # store aggregated case:
    x_str <- cbind(x_str, as.matrix(x.quarterly_group_in[,-1])%*%lb(legendre_degree[dim_x],a=0,b=1,jmax=dim(x.quarterly_group_in[,-1])[2]))
    x_str_out <- cbind(x_str_out, as.matrix(x.quarterly_group_out[,-1])%*%lb(legendre_degree[dim_x],a=0,b=1,jmax=dim(x.quarterly_group_out[,-1])[2]))
    # store unrestricted case: 
    x_unstr <- cbind(x_unstr, x.quarterly_group_in[,-1])
    x_unstr_out <- cbind(x_unstr_out, x.quarterly_group_out[,-1])
    # store averages:
    x_average <- cbind(x_average, rowMeans(x.quarterly_group_in[,-1]))
    x_average_out <- cbind(x_average_out, rowMeans(x.quarterly_group_out[,-1]))
    # store group indices: 
    group_index <- c(group_index, rep(max(group_index)+1,times=legendre_degree[dim_x]+1))
    group_index_un <- c(group_index_un, rep(max(group_index_un)+1,times=dim(x.quarterly_group_in[,-1])[2]))
    group_index_av <- c(group_index_av, max(group_index_av)+1)
  }
  # drop initializing zero:
  group_index <- group_index[-1]
  group_index_un <- group_index_un[-1]
  group_index_av <- group_index_av[-1]
  
  # sort quarterly group of covariates and lags
  y.data_in <- y.data[as.Date(y.data$DATE)%in%ref_in,]
  y.data_out <- y.data[as.Date(y.data$DATE)%in%ref_out,]
  
  # augment lags if inputed
  if(dim(y.data)[2]>2){
    y.lags_in <- y.data_in[,-c(1,2)]
    y.lags_out <- y.data_out[,-c(1,2)]
    if(group_ar_lags)
      group_ar <- rep(1,times=dim(y.lags_in)[2])
    if(!group_ar_lags)
      group_ar <- 1:dim(y.lags_in)[2]
    
    x_str <- cbind(y.lags_in, x_str)
    x_str_out <- cbind(y.lags_out, x_str_out)
    x_unstr <- cbind(y.lags_in, x_unstr)
    x_unstr_out <- cbind(y.lags_out, x_unstr_out)
    x_average <- cbind(y.lags_in, x_average)
    x_average_out <- cbind(y.lags_out, x_average_out)
    
    group_index <- c(group_ar,group_index)
    group_index_un <- c(group_ar,group_index_un)
    group_index_av <- c(group_ar,group_index_av)
  }
  y_in <- y.data_in$Y
  y_in_dates <- y.data_in$DATE
  y_out <- y.data_out$Y
  y_out_dates <- y.data_out$DATE
  
  output <- list(y_in = y_in, y_in_dates = y_in_dates, y_out = y_out, y_out_dates = y_out_dates, 
                 x_str = x_str, x_str_out = x_str_out, x_unstr = x_unstr, x_unstr_out = x_unstr_out, x_average = x_average, x_average_out = x_average_out, 
                 group_index = group_index, group_index_un = group_index_un, group_index_av = group_index_av)
  return(output)
}

  
  
  