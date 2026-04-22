# utils.R -- Utility functions for `tdcm` package.

#' Emit `tdcm`-related message
#'
#' The `tdcm_emit` function is used internally by the `tdcm` package to ensure
#' that messages, warnings, and errors are emitted in a way that can be easily
#' suppressed by users of the package.
#'
#' @param text The message text to emit, an object that can be coerced to a
#' `"character"` via [[base::as.character()]].
#'
#' @param label The label that precedes the message. The default `label` is
#' `"INFO:"`.
#'
#' @param func The function to use to emit the message. The default `func` is
#' [[base::message()]].
#'
#' @keywords internal
#' @noRd
tdcm_emit <- function(text, label = "INFO:", func = base::message, ...) {
  text <- paste("[tdcm]", label, base::as.character(text))
  func(text, ...)
} # tdcm_emit

#' @describeIn tdcm_emit Emit `tdcm`-related warning message
#' @keywords internal
#' @noRd
tdcm_warn <- function(text, ...) {
  tdcm_emit(text, label = "WARN:", func = base::warning, ...)
} # tdcm_warn

#' @describeIn tdcm_emit Emit `tdcm`-related stop message
#' @keywords internal
#' @noRd
tdcm_stop <- function(text, ...) {
  tdcm_emit(text, label = "STOP:", func = base::stop, ...)
} # tdcm_stop

#' Test for a String Value
#'
#' Internal function to test if an object is a string.
#'
#' @param x An R object to be tested.
#'
#' @param stop.when.false logical. If `TRUE`, call `tdcm_stop` when `x` is
#'     is determined to not be a string. If `FALSE` (the default), this
#'     function simply returns `FALSE`.
#'
#' @return `TRUE` if `x` is a string, else `FALSE`.
#'
#' @details The `is.string` function defines a *string* as a character vector
#'     that is length `1`. A character vector of other any other length is
#'     not considered a string by `is.string`.
#'
#' @keywords internal
#' @noRd
tdcm_is_string <- function(x, stop.when.false = FALSE) {
  x.type <- typeof(x)
  x.length <- length(x)
  if (is.character(x) & x.length == 1) {
    return(TRUE)
  } else {
    if (stop.when.false) {
      tdcm_stop(
        paste(
          "x must be a string (i.e., a character vector of length 1), but",
          "typeof(x) is", x.type, "and",
          "length(x) is", paste0(x.length, ").")
        ) # paste
      ) # tdcm_stop
    } # if
    return(FALSE)
  } # if
} # tdcm_is_string

