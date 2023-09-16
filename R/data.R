#' Corn and wheat prices from 1986 to 2014
#'
#' The to prices (in U.S. Dollars) per bushel and the log returns of corn and wheat from 1986-01-03 to 2014-10-10.
#' Each observation corresponds to the price on that day, but not all days are present in this dataset.
#'
#' @format
#' A data frame with 7,253 rows and 5 columns:
#' \describe{
#'   \item{date}{The date of the observation.}
#'   \item{corn.price, wheat.price}{ The price (in U.S. Dollars) per bushel of corn and wheat, respectively.}
#'   \item{corn.log.return, wheat.log.return }{The log returns for corn and wheat, respectively.}
#' }
#' @source \link{https://www.macrotrends.net/charts/commodities}
"cornWheat"

#' Hospital admissions by chicken pox in Brazil
#'
#' Monthly hospital admissions by chicken pox in Brazil from January 2010 to December 2019.
#'
#' @format
#' A data frame with 120 rows and 6 columns:
#' \describe{
#'   \item{date}{The date of the observations.}
#'   \item{< 5 year, 5 to 9 years, 10 to 14 years, 15 to 49 years, 50 years or more}{The number of admissions for each age group.}
#' }
#' @source \link{https://datasus.saude.gov.br/informacoes-de-saude-tabnet/}
"chickenPox"
