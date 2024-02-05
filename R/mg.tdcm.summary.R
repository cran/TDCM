#' Multigroup TDCM results compiler and summarizer.
#'
#' A function to compile results from calibration of the multigroup TDCM (Madison & Bradshaw, 2018).
#'
#' @details
#' Provides a summary of multigroup TDCM results including item parameters, attribute posterior probabilities, transition posterior probabilities, classifications, group-wise growth, group-wise transition probabilities, attribute correlations, several transition reliability metrics, and model fit. Includes longitudinal versions of reliability metrics developed by Templin and Bradshaw (2013) and Johnson and Sinharay (2020).
#'
#' @param model a \code{gdina} object returned from the \code{\link{mg.tdcm}} function.
#'
#' @param num.time.points the number of time points (i.e., measurement/testing occasions), integer \eqn{\ge 2}.
#'
#' @param transition.option option for reporting results. \code{= 1} compares the first time point to the last. \code{= 2} compares the first time point to every other time point. \code{= 3} compares successive time points. Default \code{= 1}.
#'
#' @param classthreshold probability threshold for establishing proficiency from examinee posterior probabilities. Default is .50, which maximizes overall classification accuracy. It can be set to a lower value to minimize false negatives (i.e., misclassifying proficient examinees as non-proficient) or set to a higher value to minimize false positives (i.e., misclassifying non-proficient examinees as proficient).
#'
#' @param attribute.names optional vector of attribute names to include in plots.
#'
#' @param group.names optional vector of group names to include in plots.
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
#' }
#'
#' @references Chen, J., de la Torre, J. ,& Zhang, Z. (2013). Relative and absolute fit evaluation in cognitive diagnosis modeling. \emph{Journal of Educational Measurement, 50}, 123-140.
#'
#' DiBello, L. V., Roussos, L. A., & Stout, W. F. (2007). \emph{Review of cognitively diagnostic assessment and a summary of psychometric models}. In C. R. Rao and S. Sinharay (Eds.), Handbook of Statistics, Vol. 26 (pp.979–1030). Amsterdam: Elsevier.
#'
#' Johnson, M. S., & Sinharay, S. (2020). The reliability of the posterior probability of skill attainment in diagnostic classification models. \emph{Journal of Educational Measurement, 47}(1), 5 – 31.
#'
#' Madison, M. J. (2019). Reliably assessing growth with longitudinal diagnostic classification models. \emph{Educational Measurement: Issues and Practice, 38}(2), 68-78.
#'
#' Madison, M. J., & Bradshaw, L. (2018). Evaluating intervention effects in a diagnostic classification model framework. \emph{Journal of Educational Measurement, 55}(1), 32-51.
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
#' ## Example 4: G = 2, T = 2, A = 4
#' data(data.tdcm04, package = "TDCM")
#' dat4 <- data.tdcm04$data
#' qmat4 <- data.tdcm04$q.matrix
#' group4 <- data.tdcm04$groups
#'
#' # estimate mgTDCM with invariance assumed and full LCDM
#' mg1 <- TDCM::mg.tdcm(dat4, qmat4,
#'   num.time.points = 2, rule = "GDINA",
#'   group = group4, group.invariance = TRUE, item.invariance = TRUE)
#'
#' # summarize results
#' results1 <- TDCM::mg.tdcm.summary(mg1, num.time.points = 2)
#'
#' # plot results
#' TDCM::tdcm.plot(results1)
#'
#' # estimate mgTDCM without group invariance
#' mg2 <- TDCM::mg.tdcm(dat4, qmat4,
#'   num.time.points = 2, rule = "GDINA",
#'   group = group4, group.invariance = FALSE, item.invariance = TRUE)
#'
#'
#' # compare models to assess group invariance
#' TDCM::tdcm.compare(mg1, mg2)
#' }
mg.tdcm.summary <- function(model, num.time.points, transition.option = 1, classthreshold = .50,
                            attribute.names = c(), group.names = c()) {
  if (model$progress == TRUE) {
    print("Summarizing results...", quote = FALSE)
  }

  # Initial model specifications
  N <- sum(model$N)
  num.groups <- model$G
  num.atts <- ncol(model$attribute.patt.splitted) / num.time.points # Number of Attribute Tested
  if (model$group.invariance == TRUE) {
    num.items <- length(model$itemfit.rmsea) # total items
  } else {
    num.items <- length(model$itemfit.rmsea) / num.groups # total items
  }
  items <- num.items / num.time.points # Items per time point

  # compute posterior probabilities for attributes
  postprobs <- matrix(NA, nrow = N, ncol = num.atts * num.time.points)
  for (i in 1:(num.atts * num.time.points)) {
    for (j in 1:N) {
      postprobs[j, i] <- sum(model$posterior[j, which(model$attribute.patt.splitted[, i] == 1)])
    }
  }
  colnames(postprobs) <- t(outer(c(paste("T", 1:num.time.points, sep = "")), c(paste("A", 1:num.atts, sep = "")), FUN = "paste0"))
  postprobs <- round(postprobs, 3)

  # estimated classifications
  estclass <- data.frame((postprobs > classthreshold) * 1)
  A <- num.atts # number of attributes

  # extract posteriors and profile proportions
  posteriors <- model$posterior
  est_baserates <- data.frame(model$attr.prob)


  ### Growth Block ###
  if (transition.option == 1) {
    mgo1 <- mg.summary1(
      model = model, num.atts = num.atts, num.time.points = num.time.points, num.groups = num.groups,
      attribute.names = attribute.names, group.names = group.names)
    all.trans <- mgo1$trans
    all.growth <- mgo1$growth
    rel <- mgo1$rel
  } # end if transition.option == 1, first to last

  else if (transition.option == 2) {
    mgo2 <- mg.summary2(
      model = model, num.atts = num.atts, num.time.points = num.time.points, num.groups = num.groups,
      attribute.names = attribute.names, group.names = group.names)
    all.trans <- mgo2$trans
    all.growth <- mgo2$growth
    rel <- mgo2$rel
  } # end if transition.option == 2, first to each

  else {
    mgo3 <- mg.summary3(
      model = model, num.atts = num.atts, num.time.points = num.time.points, num.groups = num.groups,
      attribute.names = attribute.names, group.names = group.names)
    all.trans <- mgo3$trans
    all.growth <- mgo3$growth
    rel <- mgo3$rel
  } # end elseif transition.option = 3, successive

  ### Parameter Block ###

  # Case 1: all invariance
  if (model$group.invariance == TRUE & model$item.invariance == TRUE) {
    p1 <- mg.summary.param1(model = model, num.time.points = num.time.points, num.atts = num.atts, num.items = num.items)
    param <- p1
  }

  # Case 2: group invariance, no time invariance
  else if (model$group.invariance == TRUE & model$item.invariance == FALSE) {
    p2 <- mg.summary.param2(model = model, num.time.points = num.time.points, num.atts = num.atts, num.items = num.items)
    param <- p2
  }


  # Case 3: time invariance, no group invariance
  else if (model$group.invariance == FALSE & model$item.invariance == TRUE) {
    p3 <- mg.summary.param3(model = model, num.time.points = num.time.points, num.atts = num.atts, num.items = num.items)
    param <- p3
  } # end case 3


  # Case 4: no group or time invariance
  else {
    p4 <- mg.summary.param3(model = model, num.time.points = num.time.points, num.atts = num.atts, num.items = num.items)
    param <- p4
  }

  # attribute correlations
  cor <- model$polychor
  for (i in 1:num.groups) {
    row.names(cor[[i]]) <- colnames(estclass)
    colnames(cor[[i]]) <- colnames(estclass)
    cor[[i]] <- round(cor[[i]], 3)
  }

  if (length(group.names) == num.groups) {
    names(cor) <- group.names
  } else {
    names(cor) <- paste("Group ", 1:num.groups, sep = "")
  }

  # add model fit
  if (model$group.invariance == TRUE) {
    mf1 <- model$itemfit.rmsea
    mf2 <- model$mean.rmsea
  } else {
    mf1 <- model$itemfit.rmsea
    mf2 <- model$mean.rmsea
  }
  if (model$progress == TRUE) {
    print("Routine finished. Check results.", quote = FALSE)
  }

  newList <- list(
    "item.parameters" = param, "growth" = all.growth, "transition.probabilities" = all.trans, "posterior.probabilities" = postprobs,
    "classifications" = estclass, "option" = transition.option, "reliability" = rel$metrics,
    "transition.posteriors" = rel$transposts, "att.corr" = cor, "most.likely.transitions" = rel$mostlikelytransitions, "item.rmsea" = mf1, "mean.item.rmsea" = mf2, "numgroups" = model$G)

  return(newList)
}
