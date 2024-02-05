#' Estimating the Transition Diagnostic Classification Model (TDCM)
#'
#' `tdcm()` is used to estimate the transition diagnostic classification model (TDCM; Madison &
#' Bradshaw, 2018a), which is a longitudinal extension of the log-linear cognitive diagnosis model
#' (LCDM; Henson, Templin, & Willse, 2009). It allows for the specification of many specific DCMs
#' via the `rule` option. For the multigroup TDCM, see [TDCM::mg.tdcm()].
#'
#' @param data A required \eqn{N \times T \times I} data matrix containing binary item responses.
#' For each time point, the binary item responses are in the columns.
#'
#' @param q.matrix A required \eqn{I \times A} matrix indicating which items measure which
#' attributes. If there are multiple Q-matrices, then they must have the same number of attributes
#' and must be stacked on top of each other for estimation (to specify multiple Q-matrices, see
#' `num.q.matrix`, `num.items`, and `anchor`).
#'
#' @param num.time.points A required integer \eqn{\ge 2} specifying the number of time points (i.e.,
#' measurement / testing occasions).
#'
#' @param invariance logical. If `TRUE` (the default), the item parameter invariance will be
#' constrained to be equal at each time point. If `FALSE`, item parameters are not assumed to be
#' equal over time.
#'
#' @param rule A string or a vector of the specific DCM to be employed. Currently accepts the
#' same values as `rule` in [CDM::gdina()]. The default is `"GDINA"`, which is implemented with a
#' logit link to estimate the LCDM. If `rule` is supplied as a single string, then that DCM will
#' be assumed for each item. If entered as a vector, a DCM can be specified for each item.
#'
#' @param num.q.matrix An optional integer specifying the number of Q-matrices. For many
#' applications, the same assessment is administered at each time point and this number is 1 (the
#' default). If there are different Q-matrices for each time point, then this argument must be
#' specified and should be equal to the number of time points. For example, if there are three time
#' points, and the Q-matrix for each time point is different, then `num.q.matrix = 3`. If there are
#' three time points, and the Q-matrix is different only for time point 3, then `num.q.matrix` is
#' still specified as `3`.
#'
#' @param num.items An optional integer specifying the number of Q-matrices (the default is `1`).
#' when there are multiple Q-matrices, the number of items in each Q-matrix is specified as a
#' vector of length `T`. For example, if there are three time points, and the Q-matrices for each
#' time point have 8, 10, and 12 items, respectively, then `num.items = c(8, 10, 12)`. Default is an
#' empty vector to indicate there is only one Q-matrix.

#' @param anchor When there are different tests at each time point, this optional anchor argument is
#'  a vector of pairs of item numbers indicating which items are the same across time points and
#'  should be held invariant. For example, if there are three Q-matrices with 10 items each, and
#'  Items 1, 11, and 21 are the same, and Items 14 and 24 are the same, then
#'  `anchor = c(1,11,1,21,14,24)`. Default is an empty vector to indicate there is only one
#'  Q-matrix.
#'
#' @param progress logical. If `FALSE` An optional logical indicating whether the function should
#' print the progress of estimation.
#'
#' @details Estimation of the TDCM via the \pkg{CDM} package (George, et al., 2016), which is based on an
#' EM algorithm as described in de la Torre (2011). The estimation approach is further detailed in
#' Madison et al. (2023).

#' @return An object of class \code{gdina} with entries as described in [CDM::gdina()]. To see a
#' TDCM-specific summary of the object (e.g., growth, transitions), use [TDCM::tdcm.summary()].
#'
#' @references
#' de la Torre, J. (2011). The generalized DINA model framework. *Psychometrika*, 76, 179-199.
#'
#' George, A. C., Robitzsch, A., Kiefer, T., Gross, J., & Ünlü , A. (2016). The R package CDM for
#' cognitive diagnosis models. *Journal of Statistical Software*, 74(2), 1-24.
#'
#' Henson, R., Templin, J., & Willse, J. (2009). Defining a family of cognitive diagnosis models
#' using log linear models with latent variables. *Psychometrika*, 74, 191-21.
#'
#' Madison, M. J., & Bradshaw, L. (2018a). Assessing growth in a diagnostic classification model
#' framework. *Psychometrika*, 82(4), 963-990.
#'
#' Madison, M. J., & Bradshaw, L. (2018b). Evaluating intervention effects in a diagnostic
#' classification model framework. *Journal of Educational Measurement*, 55(1), 32-51.
#'
#' Madison, M.J., Chung, S., Kim, J., & Bradshaw, L. (2024). Approaches to estimating longitudinal
#' diagnostic classification models. *Behaviormetrika*, 51, 7-19.
#' https://doi.org/10.1007/s41237-023-00202-5
#'
#' Rupp, A. A., Templin, J., & Henson, R. (2010).
#' *Diagnostic Measurement: Theory, Methods, and Applications*. New York: Guilford.
#'
#' @examples
#' \donttest{
#' ## Example 1: T = 2, A = 4
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#'
#' # Estimate full TDCM with invariance assumed.
#' model1 <- TDCM::tdcm(data, q.matrix, num.time.points = 2)
#'
#' # Summarize results with tdcm.summary().
#' results <- TDCM::tdcm.summary(model1, num.time.points = 2)
#' results$item.parameters
#' results$growth
#' results$transition.probabilities
#' }
#' @export
tdcm <- function(
    data,
    q.matrix,
    num.time.points,
    invariance = TRUE,
    rule = "GDINA",
    num.q.matrix = 1,
    num.items = c(),
    anchor = c(),
    progress = FALSE
) {

  if (num.q.matrix == 1) {

    if (progress) {
      tdcm_emit("Preparing data for tdcm()...")
    } # if

    # Initial Data Sorting
    n.items <- ncol(data) # Total Items
    items <- n.items / num.time.points # Items per time point
    N <- nrow(data) # Number of Examinees
    n.att <- ncol(q.matrix) # Number of Attributes

    # give names to items
    colnames(data) <- paste("Item", 1:n.items)
    rownames(q.matrix) <- paste("Item", 1:items)

    # build stacked Q-matrix
    qnew <- matrix(0, ncol = n.att * num.time.points, nrow = n.items)
    for (z in 1:num.time.points) {
      for (i in 1:items) {
        for (j in 1:n.att) {
          qnew[i + ((z - 1) * items), j + ((z - 1) * n.att)] <- q.matrix[i, j]
        } # for
      } # for
    } # for

    if (progress) {
      tdcm_emit("Estimating the TDCM in tdcm()...")
    } # if

    if (invariance == FALSE) {
      # If NOT invariant ~ no designmatrix
      tdcm <- suppressWarnings(CDM::gdina(
        data,
        qnew,
        linkfct = "logit",
        method = "ML",
        rule = rule,
        progress = FALSE
      )) # tdcm
      tdcm$invariance <- FALSE
    } else {
      # if invariance = T, then constrain item params in design matrix
      tdcm.1 <- tdcm.base(data, qnew, rule)
      c0 <- tdcm.1$coef
      c.0 <- nrow(c0)
      designmatrix <- diag(nrow = c.0 / num.time.points, ncol = c.0 / num.time.points)
      delta.designmatrix <- matrix(rep(t(designmatrix), num.time.points), ncol = ncol(designmatrix), byrow = TRUE)
      tdcm <- suppressWarnings(CDM::gdina(
        data,
        qnew,
        linkfct = "logit",
        method = "ML",
        progress = FALSE,
        delta.designmatrix = delta.designmatrix,
        rule = rule
      )) # tdcm
    } # if
  } else { # multiple Q-matrices
    tdcm <- tdcm.mq(
      data = data,
      q.matrix = q.matrix,
      num.time.points = num.time.points,
      invariance = FALSE,
      rule = rule,
      num.q.matrix = num.q.matrix,
      num.items = num.items,
      anchor = anchor
    ) # tdcm
  } # if

  # set progress value in result object
  tdcm$progress <- progress

  if (progress) {
    tdcm_emit(
      sprintf(
        "%s %s",
        "Done estimating the TDCM in tdcm().",
        "Use tdcm.summary() to display results."
      ) # sprintf
    ) # tdcm_emit
  } # if

  return(tdcm)

} # tdcm
