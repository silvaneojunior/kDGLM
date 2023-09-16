#' summary.dlm_distr
#'
#' summary method for class dlm_distr.
#'
#' @param object A fitted_dlm object.
#' @param ... Not used.
#'
#' @export
#' @keywords internal
#'
#' @family {auxiliary functions for a creating outcomes}
#' @family {Reports for dlm_distr objects.}
summary.dlm_distr <- function(object, ...) {
  report_distr(object)
}

#' summary.fitted_dlm
#'
#' summary method for class fitted_dlm.
#'
#' @param object A fitted_dlm object.
#' @param ... Arguments passed to report_dlm
#'
#' @export
#' @keywords internal
#'
#' @family {auxiliary visualization functions for the fitted.dlm class}
summary.fitted_dlm <- function(object, ...) {
  report_dlm(object, ...)
}

#' summary.dlm_block
#'
#' summary method for class dlm_block
#'
#' @param object A dlm_block object.
#' @param ... Arguments passed to report_dlm
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for structural blocks}
summary.dlm_block <- function(object, ...) {
  report_block(object, ...)
}

#' summary.searched_dlm
#'
#' summary method for class searched_dlm
#'
#' @param object A searched_dlm object.
#' @param ... Arguments passed to report_searched_dlm
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted.dlm class}
summary.searched_dlm <- function(object, ...) {
  report_searched_dlm(object, ...)
}

#' plot.fitted_dlm
#'
#' plot method for class fitted_dlm
#'
#' @param x A fitted_dlm object.
#' @param ... Arguments passed to show_fit.
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted_dlm class}
plot.fitted_dlm <- function(x, ...) {
  show_fit(x, ...)$plot
}

#' print.fitted_dlm
#'
#' @param x A fitted_dlm object.
#' @param ... Arguments passed to summary.fitted_dlm
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted_dlm class}
print.fitted_dlm <- function(x, ...) {
  summary.fitted_dlm(x, ...)
}

#' print.dlm_distr
#'
#' @param x A dlm_distr object.
#' @param ... Arguments passed to summary.dlm_distr
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted.dlm class}
print.dlm_distr <- function(x, ...) {
  summary.dlm_distr(x, ...)
}

#' print.dlm_block
#'
#' @param x A dlm_block object.
#' @param ... Arguments passed to summary.dlm_block
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for structural blocks}
print.dlm_block <- function(x, ...) {
  summary.dlm_block(x, ...)
}

#' print.searched_dlm
#'
#' @param x A searched_dlm object.
#' @param ... Arguments passed to summary.searched_dlm
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted.dlm class}
print.searched_dlm <- function(x, ...) {
  summary.searched_dlm(x, ...)
}

#' coef.fitted_dlm
#'
#' coef method for class fitted_dlm
#'
#' @param object A fitted_dlm object.
#' @param ... Arguments passed to eval_past.
#'
#' @importFrom stats coef
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for fitted_dlm objects}
coef.fitted_dlm <- function(object, ...) {
  eval_past(object, ...)
}

#' coefficients.fitted_dlm
#'
#' coefficients method for class fitted_dlm
#'
#' @param object A fitted_dlm object.
#' @param ... Arguments passed to eval_past.
#'
#' @importFrom stats coefficients
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for fitted_dlm objects}
coefficients.fitted_dlm <- function(object, ...) {
  eval_past(object, ...)
}

#' forecast.fitted_dlm
#'
#' forecast method for class fitted_dlm.
#'
#' @param object A fitted_dlm object.
#' @param ... Arguments passed to forecast_dlm
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for fitted_dlm objects}
forecast.fitted_dlm <- function(object, ...) {
  forecast_dlm(object, ...)
}

#' +.fitted_dlm
#'
#' Define add operator for class dlm_block
#'
#' @param e1 A dlm_block.
#' @param e2 A dlm_block.
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for structural blocks}
`+.dlm_block` <- function(e1, e2) {
  block_superpos(e1, e2)
}

#' *.fitted_dlm
#'
#' Define product operator for class dlm_block
#'
#' @param e1 A dlm_block (if e2 is an Integer) or an Integer (if e2 is a dlm_block).
#' @param e2 An Integer (if e1 is an dlm_block) or a dlm_block (if e1 is an Integer).
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for structural blocks}
`*.dlm_block` <- function(e1, e2) {
  if (is.numeric(e2)) {
    return(block_mult(e1, e2))
  } else {
    return(block_mult(e2, e1))
  }
}

#' Forecasting from an object
#'
#' The functions allow producing forecasts based on the provided object. This method is the same as the generics package (version 0.1.3.9000).
#'
#' @param object A model for which forecasts are required.
#' @param ... Other arguments passed to methods
#'
#' @section Methods:
#' \Sexpr[stage=render,results=rd]{generics:::methods_rd("forecast")}
#'
#' @author Hadley Wickham, \email{hadley@@rstudio.com}; Max Kuhn, \email{max@@rstudio.com}; Davis Vaughan, \email{davis@@rstudio.com}; RStudio.
#' @seealso [generics package](https://github.com/r-lib/generics)
#'
#' @details
#' MIT License
#' Copyright (c) 2020 RStudio
#'
#' Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#'
#' The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#'
#' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#'
#'
#' @export
#' @keywords internal
forecast <- function(object, ...) {
  UseMethod("forecast")
}
