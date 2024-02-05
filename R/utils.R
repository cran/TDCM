# utils.R -- Utility functions for `tdcm` package.

#' Emit `tdcm`-related message
#'
#' @description
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
tdcm_emit <- function(text, label = "INFO:", func = base::message, ...) {
  text <- paste("[tdcm]", label, base::as.character(text))
  func(text, ...)
} # tdcm_emit

#' @describeIn tdcm_emit Emit `tdcm`-related warning message
#' @keywords internal
tdcm_warn <- function(text, ...) {
  tdcm_emit(text, label = "WARN:", func = base::warning, ...)
} # tdcm_warn

#' @describeIn tdcm_emit Emit `tdcm`-related stop message
#' @keywords internal
tdcm_stop <- function(text, ...) {
  tdcm_emit(text, label = "STOP:", func = base::stop, ...)
} # tdcm_stop
