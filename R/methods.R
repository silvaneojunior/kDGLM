#' summary.dlm_distr
#'
#' summary method for class dlm_distr.
#'
#' @param object A fitted_dlm object.
#' @param ... Not used.
#'
#' @export
#' @keywords internal
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
summary.fitted_dlm <- function(object, ...) {
  report_dlm(object, ...)
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
print.dlm_distr <- function(x, ...) {
  summary.dlm_distr(x, ...)
}

#' print.searched_dlm
#'
#' @param x A searched_dlm object.
#' @param ... Arguments passed to summary.searched_dlm
#'
#' @export
#' @keywords internal
print.searched_dlm <- function(x, ...) {
  summary.searched_dlm(x, ...)
}

#' effects.fitted_dlm
#'
#' effects method for class fitted_dlm
#'
#' @param object A fitted_dlm object.
#' @param ... Arguments passed to plot_lat_var.
#'
#' @importFrom stats effects
#'
#' @export
#' @keywords internal
effects.fitted_dlm <- function(object, ...) {
  plot_lat_var(object, ...)
}

#' fitted.values.fitted_dlm
#'
#' fitted.values method for class fitted_dlm
#'
#' @param object A fitted_dlm object.
#' @param ... Arguments passed to eval_past
#'
#' @export
#' @keywords internal
fitted.values.fitted_dlm <- function(object, ...) {
  eval_past(object, ...)
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
`+.dlm_block` <- function(e1, e2) {
  block_merge(e1, e2)
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
`*.dlm_block` <- function(e1, e2) {
  if (is.numeric(e2)) {
    return(block_mult(e1, e2))
  } else {
    return(block_mult(e2, e1))
  }
}
