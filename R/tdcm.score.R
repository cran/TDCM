#' @title DCM scoring function.
#'
#' @description Function to score responses with fixed item parameters from a previously calibrated LCDM.
#'
#' @param calibration.model the previously calibrated model; an object of class \code{gdina}.
#'
#' @param newdata a required \eqn{N \times I} matrix. Binary item responses are in the columns.
#'
#' @param q.matrix a required \eqn{I \times A} matrix indicating which items measure which attributes.
#'
#' @param attr.prob.fixed optional argument for attribute profile proportions. Default is uniform distribution of profiles.
#'
#' @param progress An optional logical indicating whether the function should print the progress of estimation.
#'
#' @details Obtain classifications for new responses to items that were previously calibrated. The calibrate-and-score approach is further detailed in Madison et al. (2023).
#'
#' @return An object of class \code{gdina} with entries as indicated in the CDM package.
#'
#' @export
#'
#' @examples
#' ## Example 1: T = 2, A = 4
#' data(data.tdcm01, package = "TDCM")
#' dat1 <- data.tdcm01$data
#' qmat1 <- data.tdcm01$q.matrix
#' pre <- dat1[, 1:20]
#' post <- dat1[, 21:40]
#'
#' # calibrate LCDM with post-test data
#' m1 <- CDM::gdina(data = pre, q.matrix = qmat1, linkfct = "logit", method = "ML")
#'
#' # score pre-test responses
#' m2 <- TDCM::tdcm.score(m1, newdata = post, q.matrix = qmat1)
#' summary(m2)
#' m2$pattern
#'
#' @references
#'
#' George, A. C., Robitzsch, A., Kiefer, T., Gross, J., & Ünlü , A. (2016). The R package CDM for cognitive diagnosis models. \emph{Journal of Statistical Software, 74}(2), 1-24.
#'
#' Henson, R., Templin, J., & Willse, J. (2009). Defining a family of cognitive diagnosis models using log linear models with latent variables. \emph{Psychometrika, 74}, 191-21.
#'
#' Madison, M.J., Chung, S., Kim, J., & Bradshaw, L. (2023). Approaches to estimating longitudinal diagnostic classification models. \emph{Behaviormetrika}.
#'
tdcm.score <- function(calibration.model, newdata, q.matrix,
                       attr.prob.fixed = NULL, progress = TRUE) { # open function

  dist <- NULL

  # compute number of attributes
  natt <- ncol(q.matrix)

  if (is.null(attr.prob.fixed)) {
    dist <- rep(1 / (2^natt), 2^natt)
  }
  else {dist <- attr.prob.fixed}

  if (progress == TRUE) {
    mod <- CDM::gdina(data = newdata, q.matrix = q.matrix,
                 delta.fixed = calibration.model$delta, linkfct = "logit",
                 attr.prob.fixed = dist, progress = FALSE)
    print("Scoring responses with item parameters from previously calibrated model.",
          quote = FALSE)
    print("Scoring complete. Check results with the CDM summary function.",
          quote = FALSE)

  }
  else {
    mod <- CDM::gdina(data = newdata, q.matrix = q.matrix,
                 delta.fixed = calibration.model$delta, linkfct = "logit",
                 attr.prob.fixed = dist, progress = FALSE)
  }

  return(mod)
}
