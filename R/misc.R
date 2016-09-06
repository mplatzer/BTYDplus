
#' Wrapper for BTYD::dc.check.model.params with additional check for parameter
#' names if these are present
#'
dc.check.model.params.safe <- function(printnames, params, func) {
  # first do basic checks
  dc.check.model.params(printnames, params, func)
  # then check for names, if these are present
  if (!is.null(names(params)) && any(printnames != names(params)))
    stop("Error in ", func, ": Parameter names do not match - ", 
         paste0(printnames, collapse = ","), " != ",
         paste0(names(params), collapse = ","), 
         call. = FALSE)
}
