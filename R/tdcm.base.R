#' Utility function for base estimation in TDCMs.
#'
#' TODO update title and this description.
#'
#' @param data item response data
#' @param q.matrix.induced The induced Q-matrix.
#' @param rule A string or a vector of itemwise condensation rules that specific DCM to estimate.
#'
#' @return An object of class `gdina` returned by the internal call to [CDM::gdina()].
#'
#' @keywords internal
#' @noRd
tdcm.base <- function(data, q.matrix.induced, rule) {
  tdcm.1 <- suppressWarnings(CDM::gdina(
    data,
    q.matrix.induced,
    rule = rule,
    linkfct = "logit",  # use logit link function for all items
    method = "ML",      # directly maximize the log-likelihood function
    mono.constr = TRUE, # meet monotonicity constraints in estimation
    progress = FALSE,   # do NOT print the iteration progress
    maxit = 1           # TODO why only 1 iteration?
  )) # tdcm.1
  return(tdcm.1)
} # tdcm.base

