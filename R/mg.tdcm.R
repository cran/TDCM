#' Estimating the multigroup transition diagnostic classification model (TDCM)
#'
#' This function estimates the multigroup TDCM (Madison & Bradshaw, 2018).
#'
#' @param data a required \eqn{N \times T \times I} matrix. For each time point, binary item responses are in the columns.
#'
#' @param q.matrix a required \eqn{I \times A} matrix indicating which items measure which attributes.
#'
#' @param num.time.points the number of time points (i.e., measurement/testing occasions), integer \eqn{\ge 2}.
#'
#' @param rule the specific DCM to be employed. Currently accepts “GDINA”, “ACDM”, “DINA”, “GDINA1”, “GDINA2”, and so on. Default is “GDINA”, which is implemented with a logit link to estimate the LCDM. The “ACDM” rule will estimate the LCDM with only main effects. The “DINA” rule will estimate the DINA model. “GDINA1” will estimate the LCDM with only main effects, equivalent to “ACDM”. “GDINA2” will estimate the LCDM with up to two-way interaction effects. If rule is entered as a single string, that DCM will be assumed for each item. If entered as a vector, a DCM can be specified for each item.
#'
#' @param groups A required vector of group identifiers for multiple group estimation.
#'
#' @param group.invariance logical indicator for whether item parameter invariance should be assumed equal for all groups. Default = T. If specified as false, item parameters are not assumed equal for groups.
#'
#' @param item.invariance logical indicator for whether item parameter invariance should be constrained to be equal at each time point. Default = T. If specified as false, item parameters are not assumed equal over time.
#'
#' @param progress An optional logical indicating whether the function should print the progress of estimation.
#'
#' @return An object of class \code{gdina} with entries as indicated in the \pkg{CDM} package. For the TDCM-specific results (e.g., growth, transitions), results are summarized using the \code{\link{mg.tdcm.summary}} function.
#'
#' @note
#' Currently, this function only accepts a single Q-matrix.
#'
#' @export
#'
#' @references
#' Madison, M. J., & Bradshaw, L. (2018). Evaluating intervention effects in a diagnostic classification model framework. \emph{Journal of Educational Measurement, 55}(1), 32-51.
#'
#'
#' @examples
#' \donttest{
#' ## Example 4: G = 2, T = 2, A = 4
#' data(data.tdcm04, package = "TDCM")
#' data <- data.tdcm04$data
#' q.matrix <- data.tdcm04$q.matrix
#' groups <- data.tdcm04$groups
#'
#' # Estimate full multigroup TDCM with invariance assumed.
#' mg.model <- TDCM::mg.tdcm(data, q.matrix, num.time.points = 2, groups = groups)
#'
#' # summarize results
#' results <- TDCM::mg.tdcm.summary(mg.model, num.time.points = 2)
#'
#' # plot results
#' TDCM::tdcm.plot(results)
#' }
mg.tdcm <- function(
    data,
    q.matrix,
    num.time.points,
    rule = "GDINA",
    groups,
    group.invariance = TRUE,
    item.invariance = TRUE,
    progress = FALSE
) {

  if (progress) {
    tdcm_emit("Preparing data for mg.tdcm()...")
  } # if

  # Initial Data Sorting
  n.items <- ncol(data) # Total Items
  items <- n.items / num.time.points # Items per time point
  N <- nrow(data) # Number of Examinees
  n.att <- ncol(q.matrix) # Number of Attributes
  group.invariance <- group.invariance
  item.invariance <- item.invariance
  num.groups <- length(unique(groups))
  colnames(data) <- paste("Item", 1:n.items)

  qnew <- matrix(0, ncol = n.att * num.time.points, nrow = n.items)
  for (z in 1:num.time.points) {
    for (i in 1:items) {
      for (j in 1:n.att) {
        qnew[i + ((z - 1) * items), j + ((z - 1) * n.att)] <- q.matrix[i, j]
      }
    }
  }

  if (progress) {
    tdcm_emit("Estimating the TDCM in mg.tdcm()... This may take a few minutes...")
  } # if

  complexity <- n.att * num.time.points * num.groups
  if (complexity < 20) {
    esttime <- round(stats::runif(1, 30, 50), 0)
  } else {
    esttime <- round(stats::runif(1, 10, 25), 0)
  }

  # Case 1: all invariance
  if (group.invariance == TRUE & item.invariance == TRUE) {
    # base model, 1 iteration for design matrix
    tdcm.1 <- CDM::gdina(data, qnew,
                         linkfct = "logit", method = "ML", mono.constr = TRUE,
                         group = groups, progress = FALSE, maxit = 1, rule = rule)

    # build design matrix
    c0 <- tdcm.1$coef
    c.0 <- nrow(c0)
    designmatrix <- diag(nrow = c.0 / num.time.points, ncol = c.0 / num.time.points)
    delta.designmatrix <- matrix(rep(t(designmatrix), num.time.points), ncol = ncol(designmatrix), byrow = TRUE)

    # estimate mg tdcm
    tdcm <- CDM::gdina(data, qnew,
                       group = groups, linkfct = "logit", method = "ML",
                       delta.designmatrix = delta.designmatrix, rule = rule, progress = FALSE)
  }

  # Case 2: group invariance, no time invariance
  else if (group.invariance == TRUE & item.invariance == FALSE) {
    # estimate mg tdcm
    tdcm <- CDM::gdina(data, qnew,
                       group = groups, linkfct = "logit", method = "ML", progress = FALSE,
                       rule = rule)
  }

  # Case 3: time invariance, no group invariance
  else if (group.invariance == FALSE & item.invariance == TRUE) {
    # base model, 1 iteration for design matrix
    tdcm.1 <- CDM::gdina(data, qnew,
                         linkfct = "logit", method = "ML", mono.constr = TRUE,
                         group = groups, progress = FALSE, maxit = 1, rule = rule, invariance = FALSE)

    # build design matrix
    c0 <- tdcm.1$coef
    c.0 <- nrow(c0)
    designmatrix <- diag(nrow = c.0 / num.time.points, ncol = c.0 / num.time.points)
    delta.designmatrix <- matrix(rep(t(designmatrix), num.time.points), ncol = ncol(designmatrix), byrow = TRUE)

    # estimate mg tdcm
    tdcm <- CDM::gdina(data, qnew,
                       group = groups, linkfct = "logit", method = "ML", progress = FALSE,
                       delta.designmatrix = delta.designmatrix, rule = rule, invariance = FALSE)
  }

  # Case 4: no group or time invariance
  else {
    # estimate mg tdcm
    tdcm <- CDM::gdina(data, qnew,
                       group = groups, linkfct = "logit", method = "ML", progress = FALSE,
                       rule = rule, invariance = FALSE)
  }

  tdcm$group.invariance <- group.invariance
  tdcm$item.invariance <- item.invariance

  # set progress value in result object
  tdcm$progress <- progress

  if (progress) {
    tdcm_emit(
      sprintf(
        "%s %s",
        "Done estimating the TDCM in mg.tdcm().",
        "Use mg.tdcm.summary() to display results."
      ) # sprintf
    ) # tdcm_emit
  } # if

  return(tdcm)
}
