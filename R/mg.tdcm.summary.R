#' Multigroup TDCM results compiler and summarizer
#'
#' @description Function to summarize results obtained with the \code{\link{mg.tdcm}} function. It includes information regarding the item parameters, attribute posterior probabilities, transition posterior probabilities, attribute mastery classifications, growth, growth effects,transition probabilities, attribute correlations, model fit statistics, and several transition reliability metrics developed by Templin and Bradshaw (2013) and Johnson and Sinharay (2020).
#'
#' @param model A \code{tdcm} object returned from the \code{\link{mg.tdcm}} function.
#'
#' @param transition.option An optional argument to specify how transition probabilities should be reported for each attribute in a Q-matrix across time points.
#' - \code{transition.option = 1} (default): Summarizes the transition probabilities by comparing the first and last time point.
#' - \code{transition.option = 2}: Summarizes the transition by comparing the first time point to every subsequent time point.
#' - \code{transition.option = 3}: summarizes the transition probabilities by comparing each consecutive time point sequentially.
#'
#' @param classthreshold A numeric value between 0 and 1 specifying the probability threshold for determining examinees' proficiency based on the posterior probabilities.
#' - The default value is \code{.50}, which optimizes overall classification accuracy.
#' - Lower values reduce the probability of false negatives, such that fewer mastery examinees are misclassified as non-proficient.
#' - Higher values reduce the probability of false positives, such that fewer non-master examinees are misclassified as proficient.
#'
#' @param attribute.names An optional character \code{vector} specifying the attribute names to be included in the plots. By default, \code{attribute.names=NULL}, which uses the generic attribute labels from the Q-matrix.
#'
#' @param group.names An optional character \code{vector} specifying the group names to be included in the plots. By default, \code{group.names=NULL}, which uses a generic group label based on the number of groups in the model (e.g., Group 1, Group 2, etc).
#'
#' @return A list with the following items:
#' \itemize{
#'    \item \code{$item.parameters}: Item parameter estimates (logit) from the specified DCM.
#'
#'    \item \code{$growth}: Proficiency proportions for each time point and each attribute.
#'
#'    \item \code{$growth.effects}: It includes three growth effect size metrics for each
#'    attribute and specified transitions:
#'    1. **Growth**: Difference in proficiency proportions between the later and earlier time point.
#'    2. **Odds Ratio**: Ratio between the proficiency odds at the later time point and the proficiency odds at the earlier time point.
#'    3. **Cohen's h** (Cohen, 1988): Arcsine-transformed difference in proficiency proportions.
#'
#'    Note that the \code{growth.effect} output directly depend on the option specified in \code{transition.option}.
#'
#'    *Example*:
#'
#'    Suppose a test measures two attributes at three time points. Because there are more than two time points, the growth effect output is calculated based on the option specified in \code{transition.option}.
#'    - If \code{transition.option=1}, the growth effect for Attribute 1 and 2 is computed between time point 1 (earlier) and time point 3 (latter).
#'    - If \code{transition.option=2}, the growth effect for Attribute 1 and 2 is computed between:
#'      - Time point 1 (earlier) and time point 2 (latter).
#'      - Time point 1 (earlier) and time point 3 (latter).
#'    - If \code{transition.option=3}, the growth effect for Attribute 1 and 2 is obtained between:
#'      - Time point 1 (earlier) and time point 2 (latter).
#'      - Time point 2 (earlier), and time point 3 (latter).
#'
#'    \item \code{$transition.probabilities}: Conditional attribute proficiency transition probability matrices.
#'
#'    \item \code{$posterior.probabilities}: Examinee marginal attribute posterior probabilities of proficiency.
#'
#'    \item \code{$transition.posteriors}: Examinee marginal attribute transition posterior probabilities.
#'
#'    \item \code{$most.likely.transitions}: Examinee most likely transitions for each attribute and transition.
#'
#'    \item \code{$classifications}: Examinee classifications determined by the specified threshold applied
#'    to the posterior probabilities.
#'
#'    \item \code{$reliability}: Estimated transition reliability metrics for each attribute for the specified transitions option specified. It includes seven metrics:
#'    - **pt bis**: Longitudinal point biserial metric, which reflects the ratio between the estimated attribute proficiency base rates with the attribute proficiency posterior probabilities.
#'    - **info gain**: Longitudinal information gain metric. It quantifies how much additional information is gained regarding an attribute's transition over time.
#'    - **polychor**: Longitudinal tetrachoric metric. It quantifies how consistently an examinee transitions between mastery states between two time points.
#'    - **ave max tr**: Average maximum transition posterior metric. It quantities how likely an examinee is classified into a specific transition state over time.
#'    - **P(t > k)**: Proportion of examinees whose marginal attribute transition posteriors exceed a threshold *k*. The thresholds used are *k* = 0.6, 0.7, 0.8, and 0.9, representing the proportion of examinees with attribute transition posterior probabilities greater than these values. For example, if P(t>.6) = 0.90, 90% of examinees have a posterior probability greater than 0.6.
#'    - **wt pt bis**: Weighted longitudinal point biserial. A variation of the longitudinal point biserial metric that computes the correlation between true attribute transition classification and observed marginal transition probabilities. It assigns greater weight to more prevalent attribute transitions based on each attributes' transition base rate, ensuring that transitions occurring more frequently in the data contribute more significantly to the computed reliability value.
#'    - **wt info gain**: Weighted longitudinal information gain. A variation of the longitudinal information gain that quantifies the additional information provided by the attribute transition posterior probabilities in predicting examinees' true transition status. It assigns greater weight to more prevalent attribute transitions, ensuring that transitions occurring more frequently in the data contribute more significantly to the computed reliability value.
#'
#'    \item \code{$att.corr}: Estimated attribute correlation matrix.
#'
#'    \item \code{$model.fit}: Several model fit indices and tests are output including:
#'    - Item root mean square error of approximation (RMSEA; von Davier, 2005).
#'    - Mean RMSEA.
#'    - Bivariate item fit statistics (Chen et al., 2013).
#'    - Absolute fit statistics such as mean absolute deviation for observed.
#'    - Expected item correlations (MADcor; DiBello, Roussos, & Stout, 2007).
#'    - Standardized root mean square root of squared residuals (SRMSR; Maydeu-Olivares, 2013).
#'
#' }
#'
#' @references Chen, J., de la Torre, J. ,& Zhang, Z. (2013). Relative and absolute fit evaluation in cognitive
#' diagnosis modeling. \emph{Journal of Educational Measurement, 50}, 123-140.
#'
#' DiBello, L. V., Roussos, L. A., & Stout, W. F. (2007). \emph{Review of cognitively diagnostic assessment and
#' a summary of psychometric models}. In C. R. Rao and S. Sinharay (Eds.), Handbook of Statistics, Vol. 26
#' (pp.979–1030). Amsterdam: Elsevier.
#'
#' Johnson, M. S., & Sinharay, S. (2020). The reliability of the posterior probability of skill attainment
#' in diagnostic classification models. \emph{Journal of Educational Measurement, 47}(1), 5 – 31.
#'
#' Madison, M. J. (2019). Reliably assessing growth with longitudinal diagnostic classification models.
#' \emph{Educational Measurement: Issues and Practice, 38}(2), 68-78.
#'
#' Madison, M. J., & Bradshaw, L. (2018). Evaluating intervention effects in a diagnostic classification model
#' framework. \emph{Journal of Educational Measurement, 55}(1), 32-51.
#'
#' Maydeu-Olivares, A. (2013). Goodness-of-fit assessment of item response theory models (with discussion).
#' \emph{Measurement: Interdisciplinary Research and Perspectives, 11}, 71-137.
#'
#' Schellman, M., & Madison, M. J. (2024). Estimating the reliability of skill transition in longitudinal DCMs.
#'  \emph{Journal of Educational and Behavioral Statistics}.
#'
#' Templin, J., & Bradshaw, L. (2013). Measuring the reliability of diagnostic classification model examinee
#' estimates. \emph{Journal of Classification, 30}, 251-275.
#'
#' von Davier M. (2008). A general diagnostic model applied to language testing data. \emph{The British journal
#' of mathematical and statistical psychology, 61}(2), 287–307.
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' ### ADD EXAMPLE WITH DIFFERENT TRANSITION OPTION --> there is no dataset for this
#'
#' ############################################################################
#' # Example 1: Multigroup TDCM assuming time and group invariance
#' ############################################################################
#'
#' # Load data: G = 2, T = 2, A = 4, I = 20
#' data(data.tdcm04, package = "TDCM")
#' data <- data.tdcm04$data
#' q.matrix <- data.tdcm04$q.matrix
#' groups <- data.tdcm04$groups
#'
#' # Estimate model
#' mg.model1 <- TDCM::mg.tdcm(data, q.matrix, num.time.points = 2,
#'                            groups = groups, time.invariance = TRUE,
#'                            group.invariance = TRUE)
#'
#'#----------------------------------------------------------------------------
#'# With different thresholds
#'#----------------------------------------------------------------------------
#' ## a) If classthreshold = 0.5 (default)
#'
#' # Summarize results
#' results1 <- TDCM::mg.tdcm.summary(mg.model1, transition.option = 1)
#' head(results1$posterior.probabilities)
#'#       T1A1  T1A2  T1A3  T1A4  T2A1  T2A2  T2A3  T2A4
#'# [1,] 0.000 0.000 0.000 0.002 0.343 0.009 0.961 0.908
#'# [2,] 0.000 0.013 0.000 0.013 0.009 0.006 0.001 0.006
#'# [3,] 0.029 0.013 0.000 0.000 0.006 0.004 0.848 0.004
#'# [4,] 0.020 0.015 0.006 0.651 0.999 0.998 0.973 1.000
#'# [5,] 0.660 0.978 0.002 0.006 0.137 0.779 0.025 1.000
#'# [6,] 0.999 0.998 0.998 1.000 0.001 1.000 0.990 0.894
#'
#' head(results1$classifications)
#'#      T1A1 T1A2 T1A3 T1A4 T2A1 T2A2 T2A3 T2A4
#'#    1    0    0    0    0    0    0    1    1
#'#    2    0    0    0    0    0    0    0    0
#'#    3    0    0    0    0    0    0    1    0
#'#    4    0    0    0    1    1    1    1    1
#'#    5    1    1    0    0    0    1    0    1
#'#    6    1    1    1    1    0    1    1    1
#'
#' ## b) If classthreshold = 0.7
#'
#' # Summarize results
#' results2 <- TDCM::mg.tdcm.summary(mg.model1, transition.option = 1, classthreshold = 0.7)
#' head(results2$posterior.probabilities)
#'#       T1A1  T1A2  T1A3  T1A4  T2A1  T2A2  T2A3  T2A4
#'# [1,] 0.000 0.000 0.000 0.002 0.343 0.009 0.961 0.908
#'# [2,] 0.000 0.013 0.000 0.013 0.009 0.006 0.001 0.006
#'# [3,] 0.029 0.013 0.000 0.000 0.006 0.004 0.848 0.004
#'# [4,] 0.020 0.015 0.006 0.651 0.999 0.998 0.973 1.000
#'# [5,] 0.660 0.978 0.002 0.006 0.137 0.779 0.025 1.000
#'# [6,] 0.999 0.998 0.998 1.000 0.001 1.000 0.990 0.894
#'
#' head(results2$classifications)
#'#       T1A1 T1A2 T1A3 T1A4 T2A1 T2A2 T2A3 T2A4
#'#    1    0    0    0    0    0    0    1    1
#'#    2    0    0    0    0    0    0    0    0
#'#    3    0    0    0    0    0    0    1    0
#'#    4    0    0    0    0    1    1    1    1
#'#    5    0    1    0    0    0    1    0    1
#'#    6    1    1    1    1    0    1    1    1
#'
#' ## c) If classthreshold = 0.3
#'
#' # Summarize results
#' results3 <- TDCM::mg.tdcm.summary(mg.model1, transition.option = 1, classthreshold = 0.3)
#' head(results3$posterior.probabilities)
#'#       T1A1  T1A2  T1A3  T1A4  T2A1  T2A2  T2A3  T2A4
#'# [1,] 0.000 0.000 0.000 0.002 0.343 0.009 0.961 0.908
#'# [2,] 0.000 0.013 0.000 0.013 0.009 0.006 0.001 0.006
#'# [3,] 0.029 0.013 0.000 0.000 0.006 0.004 0.848 0.004
#'# [4,] 0.020 0.015 0.006 0.651 0.999 0.998 0.973 1.000
#'# [5,] 0.660 0.978 0.002 0.006 0.137 0.779 0.025 1.000
#'# [6,] 0.999 0.998 0.998 1.000 0.001 1.000 0.990 0.894
#'
#' head(results3$classifications)
#'#      T1A1 T1A2 T1A3 T1A4 T2A1 T2A2 T2A3 T2A4
#'#   1    0    0    0    0    1    0    1    1
#'#   2    0    0    0    0    0    0    0    0
#'#   3    0    0    0    0    0    0    1    0
#'#   4    0    0    0    1    1    1    1    1
#'#   5    1    1    0    0    0    1    0    1
#'#   6    1    1    1    1    0    1    1    1
#'
#' }


mg.tdcm.summary <- function(model, transition.option = 1, classthreshold = .50,
                            attribute.names = c(), group.names = c()) {
  if (model$progress == TRUE) {
    print("Summarizing results...", quote = FALSE)
  }

  # Initial model specifications
  num.time.points = model$numtimepoints
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
  colnames(postprobs) <- t(outer(c(paste("T", 1:num.time.points, sep = "")), c(paste("A", 1:num.atts, sep = "")),
                                 FUN = "paste0"))
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
    growth.effects <- mgo1$growth.effects
  } # end if transition.option == 1, first to last

  else if (transition.option == 2) {
    mgo2 <- mg.summary2(
      model = model, num.atts = num.atts, num.time.points = num.time.points, num.groups = num.groups,
      attribute.names = attribute.names, group.names = group.names)
    all.trans <- mgo2$trans
    all.growth <- mgo2$growth
    rel <- mgo2$rel
    growth.effects <- mgo2$growth.effects

  } # end if transition.option == 2, first to each

  else {
    mgo3 <- mg.summary3(
      model = model, num.atts = num.atts, num.time.points = num.time.points, num.groups = num.groups,
      attribute.names = attribute.names, group.names = group.names)
    all.trans <- mgo3$trans
    all.growth <- mgo3$growth
    rel <- mgo3$rel
    growth.effects <- mgo3$growth.effects

  } # end elseif transition.option = 3, successive

  ### Parameter Block ###

  # Case 1: all invariance
  if (model$group.invariance == TRUE & model$time.invariance == TRUE) {
    p1 <- mg.summary.param1(model = model, num.time.points = num.time.points, num.atts = num.atts, num.items =
                              num.items)
    param <- p1
  }

  # Case 2: group invariance, no time invariance
  else if (model$group.invariance == TRUE & model$time.invariance == FALSE) {
    p2 <- mg.summary.param2(model = model, num.time.points = num.time.points, num.atts = num.atts, num.items =
                              num.items)
    param <- p2
  }


  # Case 3: time invariance, no group invariance
  else if (model$group.invariance == FALSE & model$time.invariance == TRUE) {
    p3 <- mg.summary.param3(model = model, num.time.points = num.time.points, num.atts = num.atts, num.items =
                              num.items)
    param <- p3
  } # end case 3


  # Case 4: no group or time invariance
  else {
    p4 <- mg.summary.param4(model = model, num.time.points = num.time.points, num.atts = num.atts, num.items =
                              num.items)
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
  # if (model$group.invariance == TRUE) {
  #   mf1 <- model$itemfit.rmsea
  #   mf2 <- model$mean.rmsea
  # } else {
  #   mf1 <- model$itemfit.rmsea
  #   mf2 <- model$mean.rmsea
  # }


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
    "item.parameters" = param, "growth" = all.growth,
    "growth.effects" = growth.effects,
    "transition.probabilities" = all.trans,
    "posterior.probabilities" = postprobs,
    "classifications" = estclass, "option" = transition.option, "reliability" = rel$metrics,
    "transition.posteriors" = rel$transposts, "att.corr" = cor, "most.likely.transitions" =
      rel$mostlikelytransitions, "model.fit" = mf, "numgroups" = model$G)

  return(newList)
}
