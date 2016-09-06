
#' Wrapper for BTYD::dc.check.model.params with additional check for parameter
#' names if these are present
#'
#' @param printnames Names to print parameter errors
#' @param params model parameters
#' @param func Function calling dc.check.model.params.safe
#' @return stops program if there is something wrong with the parameters
dc.check.model.params.safe <- function(printnames, params, func) {
  # first do basic checks
  dc.check.model.params(printnames, params, func)
  # then check for names, if these are present
  if (!is.null(names(params))) {
    idx <- names(params) != ""
    if (any(printnames[idx] != names(params)[idx])) {
      stop("Error in ", func, ": Parameter names do not match - ", paste0(printnames, collapse = ","), " != ", paste0(names(params), 
        collapse = ","), call. = FALSE)
    }
  }
}
