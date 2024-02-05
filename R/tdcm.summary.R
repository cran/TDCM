#' @title TDCM results compiler and summarizer.
#'
#' @description Function to summarize results from TDCM analyses.
#'
#' @details
#' Provides a summary of TDCM results including item parameters, attribute posterior
#' probabilities, transition posterior probabilities, classifications, growth,
#' transition probabilities, attribute correlations, several transition reliability metrics,
#' and model fit. Includes longitudinal DCM reliability metrics developed by
#' Schellman and Madison (2021).
#'
#'
#' @param model a \code{gdina} object returned from the \code{\link{tdcm}} function.
#'
#' @param num.time.points the number of time points (i.e., measurement/testing occasions), integer \eqn{\ge 2}.
#'
#' @param transition.option option for reporting results. \code{= 1} compares the first time point to the last. \code{= 2} compares the first time point to every other time point. \code{= 3} compares successive time points. Default = 1.
#'
#' @param classthreshold probability threshold for establishing proficiency from examinee posterior probabilities. Default is .50, which maximizes overall classification accuracy. It can be set to a lower value to minimize false negatives (i.e., misclassifying proficient examinees as non-proficient) or set to a higher value to minimize false positives (i.e., misclassifying non-proficient examinees as proficient).
#'
#' @param attribute.names optional vector of attribute names to include in results output.
#'
#' @return A list with the following items:
#' \itemize{
#'    \item \code{$item.parameters}: LCDM item parameter estimates from the specified DCM.
#'
#'    \item \code{$growth}: proficiency proportions for each time point and each attribute
#'
#'    \item \code{$transition.probabilities}: conditional attribute proficiency transition probability matrices
#'
#'    \item \code{$posterior.probabilities}: examinee marginal attribute posterior probabilities of proficiency
#'
#'    \item \code{$transition.posteriors}: examinee marginal attribute transition posterior probabilities
#'
#'    \item \code{$most.likely.transitions}: examinee most likely transitions for each attribute and transition
#'
#'    \item \code{$classifications}: examinee classifications determined by the specified threshold applied to the posterior probabilities
#'
#'    \item \code{$reliability}: estimated transition reliability metrics for each attribute for the specified transitions. “pt bis” = longitudinal point biserial metric; “info gain” = longitudinal information gain metric; “polychor” = longitudinal tetrachoric metric; “ave max tr” = average maximum transition posterior metric; “P(t>k)” = proportion of examinee marginal attribute transition posteriors greater than k; “wt pt bis” = weighted longitudinal point biserial; “wt info gain” = weighted longitudinal information gain.
#'
#'    \item \code{$att.corr}: estimated attribute correlation matrix
#'
#'    \item \code{$model.fit}: Several model fit indices and tests are output including item root mean square error of approximation (RMSEA; von Davier, 2005), mean RMSEA, bivariate item fit statistics (Chen et al., 2013), and absolute fit statistics such as mean absolute deviation for observed and expected item correlations (MADcor; DiBello, Roussons, & Stout, 2007), and standardized root mean square root of squared residuals (SRMSR; Maydeu-Olivares, 2013)
#'
#' }
#'
#'
#' @references Chen, J., de la Torre, J. ,& Zhang, Z. (2013). Relative and absolute fit evaluation in cognitive diagnosis modeling. \emph{Journal of Educational Measurement, 50}, 123-140.
#'
#' DiBello, L. V., Roussos, L. A., & Stout, W. F. (2007). \emph{Review of cognitively diagnostic assessment and a summary of psychometric models}. In C. R. Rao and S. Sinharay (Eds.), Handbook of Statistics, Vol. 26 (pp.979–1030). Amsterdam: Elsevier.
#'
#' Johnson, M. S., & Sinharay, S. (2020). The reliability of the posterior probability of skill attainment in diagnostic classification models. \emph{Journal of Educational Measurement, 47}(1), 5 – 31.
#'
#' Madison, M. J. (2019). Reliably assessing growth with longitudinal diagnostic classification models. \emph{Educational Measurement: Issues and Practice, 38}(2), 68-78.
#'
#' Maydeu-Olivares, A. (2013). Goodness-of-fit assessment of item response theory models (with discussion). \emph{Measurement: Interdisciplinary Research and Perspectives, 11}, 71-137.
#'
#' Schellman, M., & Madison, M. J. (2021, July). \emph{Estimating the reliability of skill transition in longitudinal DCMs}. Paper presented at the 2021 International Meeting of the Psychometric Society.
#'
#' Templin, J., & Bradshaw, L. (2013). Measuring the reliability of diagnostic classification model examinee estimates. \emph{Journal of Classification, 30}, 251-275.
#'
#' von Davier M. (2008). A general diagnostic model applied to language testing data. \emph{The British journal of mathematical and statistical psychology, 61}(2), 287–307.

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
#' # summarize results with tdcm.summary function
#' results1 <- TDCM::tdcm.summary(m1, num.time.points = 2)
#' results1$item.parameters
#' results1$growth
#' results1$transition.probabilities
#' results1$reliability
#' head(results1$most.likely.transitions)
#' results1$model.fit$Item.RMSEA
#' }
tdcm.summary <- function(model, num.time.points, transition.option = 1, classthreshold = .50,
                         attribute.names = c()) {

  num.items <- length(model$itemfit.rmsea) # total items
  items <- num.items / num.time.points # Items per time point
  num.atts <- ncol(model$attribute.patt.splitted) / num.time.points # Number of attribute measured

  # posterior probabilities
  postprobs <- matrix(NA, nrow = model$N, ncol = num.atts * num.time.points)
  postprobs <- round(model$pattern[, 6:(6 + num.atts * num.time.points - 1)], 3)
  colnames(postprobs) <- t(outer(c(paste("T", 1:num.time.points, sep = "")), c(paste("A", 1:num.atts, sep = "")), FUN = "paste0"))

  # estimated classifications
  estclass <- data.frame((postprobs > classthreshold) * 1)
  A <- num.atts # number of attributes

  # extract posteriors and profile proportions
  posteriors <- model$posterior
  est_baserates <- data.frame(model$attribute.patt$class.prob)

  if (model$invariance == FALSE) {
    invariance <- FALSE
  } else {
    invariance <- TRUE
  }

  ### Growth Block ###
  if (transition.option == 1) {
    o1 <- summary.option1(model = model, num.atts = num.atts, num.time.points = num.time.points, attribute.names = attribute.names)
    trans <- o1$trans
    growth <- o1$growth
    rel <- o1$rel
  } else if (transition.option == 2) {
    o2 <- summary.option2(model = model, num.atts = num.atts, num.time.points = num.time.points, attribute.names = attribute.names)
    trans <- o2$trans
    growth <- o2$growth
    rel <- o2$rel
  } else {
    o3 <- summary.option3(model = model, num.atts = num.atts, num.time.points = num.time.points, attribute.names = attribute.names)
    trans <- o3$trans
    growth <- o3$growth
    rel <- o3$rel
  }

  ### Parameters Block ###
  param <- summary.param(model = model, num.atts = num.atts, num.items = num.items, num.time.points = num.time.points)

  # attribute correlations
  cor <- model$polychor
  row.names(cor[[1]]) <- colnames(estclass)
  colnames(cor[[1]]) <- colnames(estclass)
  cor <- round(cor[[1]], 3)

  # model fit
  mf <- CDM::modelfit.cor.din(model)
  mf2 <- list(model$itemfit.rmsea)
  names(mf2[[1]]) <- paste("Item", 1:length(mf2[[1]]))
  mf3 <- list(model$mean.rmsea)
  mf$itempairs$item1 <- stringr::str_replace(mf$itempairs$item1, "V", "Item ")
  mf$itempairs$item2 <- stringr::str_replace(mf$itempairs$item2, "V", "Item ")
  mf <- append(mf, mf2)
  mf <- append(mf, mf3)
  mf <- append(mf, c(model$loglike, model$deviance, model$AIC, model$BIC, model$CAIC, model$Npars))

  names(mf) <- c(
    "Global.Fit.Stats", "Item.Pairs", "Global.Fit.Tests", "Global.Fit.Stats2", "Item.RMSEA", "Mean.Item.RMSEA",
    "loglike", "deviance", "AIC", "BIC", "CAIC", "Npars")

  if (model$progress == TRUE) {
    print("Routine finished. Check results.", quote = FALSE)
  }

  newList <- list(
    "item.parameters" = param, "growth" = growth, "transition.probabilities" = trans,
    "posterior.probabilities" = postprobs, "classifications" = estclass, "reliability" = rel$metrics,
    "transition.posteriors" = rel$transposts, "most.likely.transitions" = rel$mostlikelytransitions,
    "option" = transition.option, "att.corr" = cor, "model.fit" = mf, "numgroups" = 1)
  return(newList)
}
