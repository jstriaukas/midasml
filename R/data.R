#' US real GDP data with several high-frequency predictors
#'
#' @name us_rgdp
#' @description 
#'  US real GDP, Chicago National Activity Index, Nonfarm payrolls and ADS Index
#' @docType data
#' @usage data(us_rgdp)
#' @format A \code{\link{list}} object.a
#' @keywords datasets
#' @source 
#'   \href{https://fred.stlouisfed.org/series/GDPC1}{rgdp} \cr
#'   \href{https://www.chicagofed.org/research/data/cfnai/current-data}{cfnai} \cr
#'   \href{https://fred.stlouisfed.org/series/PAYEMS}{payems} \cr
#'   \href{https://www.philadelphiafed.org/research-and-data/real-time-center/business-conditions-index}{ads} 
#' @examples
#' \donttest{
#' data(us_rgdp)
#' us_rgdp$rgdp # - GDP data
#' us_rgdp$cfnai # - CFNAI predictor data
#' us_rgdp$payems # - Nonfarm payrolls predictor data 
#' us_rgdp$ads # - ADS predictor data
#' }
NULL

#' SNP500 returns
#'
#' @name market_ret
#' @description 
#'  SNP500 returns
#' @docType data
#' @usage data(market_ret)
#' @format A \code{\link{data.frame}} object.a
#' @keywords datasets
#' @source 
#'   \code{market_ret} - \href{https://fred.stlouisfed.org/series/SP500}{FRED} 
#' @examples
#' \donttest{
#' data(market_ret)
#' market_ret$snp500ret
#' }
NULL

#' Real GDP vintages
#'
#' @name rgdp_vintages
#' @description 
#'  Real GDP vintages
#' @docType data
#' @usage data(rgdp_vintages)
#' @format A \code{\link{list}} objects
#' @keywords datasets
#' @source 
#'   \href{https://alfred.stlouisfed.org/}{ALFRED} \cr
#' @examples
#' \donttest{
#' data(rgdp_vintages)
#' rgdp_vintages$date # dates
#' rgdp_vintages$time_series # series, q-q annual rate
#' rgdp_vintages$realtime_period # real time dates 
#' }
NULL


#'  ALFRED monthly and quarterly series vintages
#'
#' @name alfred_vintages
#' @description 
#'  ALFRED monthly and quarterly series vintages
#' @docType data
#' @usage data(alfred_vintages)
#' @format A \code{\link{list}} objects
#' @keywords datasets
#' @source 
#'   \href{https://alfred.stlouisfed.org/}{ALFRED} \cr
#' @examples
#' \donttest{
#' data(alfred_vintages)
#' i <- 1
#' alfred_vintages[[i]] # ith variable
#' }
NULL


#' Real GDP release dates
#'
#' @name rgdp_dates
#' @description 
#'  Real GDP release dates
#' @docType data
#' @usage data(rgdp_dates)
#' @format A \code{\link{list}} objects
#' @keywords datasets
#' @source 
#'   \href{https://alfred.stlouisfed.org/}{ALFRED} \cr
#' @examples
#' \donttest{
#' data(rgdp_dates)
#' rgdp_dates$Quarter_q # reference quarters in quarters
#' rgdp_dates$Quarter_m # reference quarters in months
#' rgdp_dates$Quarter_d # reference quarters in days
#' rgdp_dates$`First release` # first release date for the reference
#' rgdp_dates$`Second release` # second release date for the reference
#' rgdp_dates$`Third release` # third release date for the reference
#' }
NULL




