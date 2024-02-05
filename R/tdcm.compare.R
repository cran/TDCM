#' @title Comparing the fit of two TDCMs
#'
#' @description
#' Provides a comparison of two TDCMs. Can be used to compare different measurement models or assess measurement invariance over time or over groups in the multigroup TDCM case. Only accepts two models.
#'
#'
#' @param model1 a \code{gdina} object returned from the \code{\link{tdcm}} or \code{\link{mg.tdcm}} function.
#' @param model2 a second \code{gdina} object returned from the \code{\link{tdcm}} or \code{\link{mg.tdcm}} function
#'
#' @note
#' - Currently, this function currently accepts two models for comparison.
#' - Both models must be fit to the same item responses and Q-matrix.
#' - The function will provide results for two non-nested models. Please ensure that models are nested before interpreting the likelihood ratio test for nested models.
#' - The likelihood ratio test is not valid for some model comparisons (e.g., LCDM vs DINA) because of model constraints.
#'
#' @return
#' This function returns a data frame with model fit statistics (AIC/BIC) and results from a likelihood ratio or deviance test.
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Example 1: T = 2, A = 4
#' data(data.tdcm01, package = "TDCM")
#' dat1 <- data.tdcm01$data
#' qmat1 <- data.tdcm01$q.matrix
#'
#' # estimate TDCM with invariance assumed and full LCDM
#' m1 <- TDCM::tdcm(dat1, qmat1, num.time.points = 2, invariance = TRUE, rule = "GDINA")
#'
#' # estimate TDCM with invariance not assumed
#' m2 <- TDCM::tdcm(dat1, qmat1, num.time.points = 2, invariance = FALSE, rule = "GDINA")
#'
#' # compare models to assess measurement invariance.
#' TDCM::tdcm.compare(m1, m2)
#' }
tdcm.compare <- function(model1, model2) {
  # Model names
  m1name <- substitute(model1)
  m2name <- substitute(model2)

  # log likelihoods
  ll1 <- round(model1$loglike, 2)
  ll2 <- round(model2$loglike, 2)

  # deviance
  dev1 <- round(model1$deviance, 2)
  dev2 <- round(model2$deviance, 2)

  # deviance
  np1 <- model1$Npars
  np2 <- model2$Npars

  # AIC/BIC
  aic1 <- round(model1$AIC, 2)
  aic2 <- round(model2$AIC, 2)
  bic1 <- round(model1$BIC, 2)
  bic2 <- round(model2$BIC, 2)

  # chi2
  chi1 <- round(-2 * (ll1 - ll2), 4)
  if (chi1 < 0) {
    chi1 <- -chi1
  }
  chi2 <- NA

  # df
  df1 <- np2 - np1
  if (df1 < 0) {
    df1 <- -df1
  }
  df2 <- NA

  # pvalue
  p1 <- round(1 - stats::pchisq(chi1, df1), 4)
  p2 <- NA
  results <- matrix(c(m1name, m2name, ll1, ll2, dev1, dev2, np1, np2, aic1, aic2, bic1, bic2, chi1, chi2, df1, df2, p1, p2), nrow = 2)
  colnames(results) <- c("Model", "loglike", "Deviance", "Npars", "AIC", "BIC", "Chisq", "df", "p")
  rownames(results) <- c(1, 2)
  results <- as.data.frame(results)
  return(results)
}
