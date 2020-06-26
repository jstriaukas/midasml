#' US real GDP data with several high-frequency predictors
#'
#' @name us_rgdp
#' @description 
#'  US real GDP, Chicago National Activity Index, Nonfarm payrolls and ADS Index
#' @docType data
#' @usage data(us_rgdp)
#' @format A \code{\link{data.frame}} object.a
#' @keywords datasets
#' @source 
#'   \code{rgdp} - \href{https://fred.stlouisfed.org/}{FRED} \cr
#'   \code{cfnai} - \href{https://www.chicagofed.org/research/data/cfnai/current-data}{Chicago Fed website} \cr
#'   \code{payems} - \href{https://fred.stlouisfed.org/}{FRED} \cr
#'   \code{ads} - \href{https://www.philadelphiafed.org/research-and-data/real-time-center/business-conditions-index}{Philadelphia Fed website} 
#' @examples
#' data(us_rgdp)
#' us_rgdp$rgdp # - GDP data
#' us_rgdp$cfnai # - CFNAI predictor data
#' us_rgdp$payems # - Nonfarm payrolls predictor data 
#' us_rgdp$ads # - ADS predictor data
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
#'   \code{rgdp} - \href{https://fred.stlouisfed.org/}{FRED} 
#' @examples
#' data(market_ret)
#' market_ret$snp500ret
NULL

#' GDP nowcasting using midasML approach example data
#'
#' @name macro_midasml
#' @description 
#'  US real GDP, FRED-MD monthly dataset, SPF survey and textual analysis data.
#' @docType data
#' @usage data(market_ret)
#' @format A \code{\link{list}} object
#' @keywords datasets
#' @source 
#'   \code{macro_midasml$rgdp} - \href{https://fred.stlouisfed.org/}{FRED} \cr
#'   \code{macro_midasml$md.data} - \href{https://research.stlouisfed.org/econ/mccracken/fred-databases/}{St. Louis Fed M. W. McCracken website} \cr
#'   \code{macro_midasml$text.data} - \href{http://www.structureofnews.com/}{The Structure of Economic News website} \cr
#'   \code{macro_midasml$survey.data} - \href{https://philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters}{Philadelphia Fed SPF data} 
#' @examples
#' data(macro_midasml)
#' macro_midasml$rgdp.data # GDP data
#' macro_midasml$md.data # FRED-MD data
#' macro_midasml$text.data # textual analysis data
#' macro_midasml$survey.data # SPF data
NULL


