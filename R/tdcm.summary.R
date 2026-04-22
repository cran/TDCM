#' @title TDCM results compiler and summarizer.
#'
#' @description Function to summarize results obtained with the \code{\link{tdcm}} function. It includes information regarding the item parameters, attribute posterior probabilities, transition posterior probabilities, attribute mastery classifications, growth, growth effects,transition probabilities, attribute correlations, model fit statistics, and several transition reliability metrics developed by Schellman and Madison (2024).
#'
#' @param model A \code{tdcm} object returned from the \code{\link{tdcm}} function.
#'
#' @param transition.option An optional argument to specify how growth and transition probabilities should be reported for each attribute across time points.
#' - \code{transition.option = 1} (default): Summarizes the transition probabilities by comparing the first and last time point.
#' - \code{transition.option = 2}: Summarizes the transition by comparing the first time point to every subsequent time point.
#' - \code{transition.option = 3}: Summarizes the transition probabilities by comparing each consecutive time point sequentially.
#'
#' @param classthreshold A numeric value between 0 and 1 specifying the probability threshold for determining examinees' proficiency based on the posterior probabilities.
#' - The default value is \code{.50}, which optimizes overall classification accuracy.
#' - Lower values reduce the probability of false negatives, such that fewer proficient examinees are misclassified as non-proficient.
#' - Higher values reduce the probability of false positives, such that fewer non-proficient examinees are misclassified as proficient.
#'
#' @param attribute.names An optional character \code{vector} specifying the attribute names to be included in the results outputs. By default, \code{attribute.names=NULL}, which uses the generic attribute labels from the Q-matrix.
#'
#' @return A list with the following items:
#' \itemize{
#'    \item \code{$item.parameters}: Item parameter estimates from the specified DCM.
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
#'    - If \code{transition.option=1}, the growth effect for Attribute 1 and 2 is computed between Time Point 1 (first) and Time Point 3 (last).
#'    - If \code{transition.option=2}, the growth effect for Attribute 1 and 2 is computed between:
#'      - Time Point 1 (first) and Time Point 2 (latter).
#'      - Time Point 1 (first) and Time Point 3 (latter).
#'    - If \code{transition.option=3}, the growth effect for Attribute 1 and 2 is obtained between:
#'      - Time Point 1 (earlier) and Time Point 2 (next).
#'      - Time Point 2 (earlier), and Time Point 3 (next).
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
#'    \item \code{$reliability}: Estimated transition reliability metrics for each attribute for the specified transitions option specified (Madison, 2019; Schellman & Madison, 2024). It includes seven metrics:
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
#'    - Mean item RMSEA.
#'    - Bivariate item fit statistics (Chen et al., 2013).
#'    - Absolute fit statistics such as mean absolute deviation for observed.
#'    - Expected item correlations (MADcor; DiBello, Roussos, & Stout, 2007).
#'    - Standardized root mean square root of squared residuals (SRMSR; Maydeu-Olivares, 2013).
#'
#' }
#'
#' @references Chen, J., de la Torre, J., & Zhang, Z. (2013). Relative and absolute fit evaluation in cognitive
#' diagnosis modeling. \emph{Journal of Educational Measurement, 50}, 123-140.
#'
#' Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral Sciences} (2nd ed.). Hillsdale, NJ:
#' Lawrence Erlbaum Associates, Publishers.
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
#'############################################################################
#'# Example 1: TDCM with full measurement invariance and equal Q-matrix
#'############################################################################
#'
#' # Load data: T = 2, A = 4
#' data(data.tdcm01, package = "TDCM")
#' dat1 <- data.tdcm01$data
#' qmat1 <- data.tdcm01$q.matrix
#'
#' # Estimate model
#'
#' model1 <- TDCM::tdcm(dat1, qmat1, num.time.points = 2, invariance = TRUE, rule = "LCDM")
#'
#' # summarize results with tdcm.summary function
#' results1 <- TDCM::tdcm.summary(model1, transition.option = 1)
#' results1$item.parameters
#' results1$growth
#' results1$growth.effects
#' results1$transition.probabilities
#' results1$reliability
#' head(results1$most.likely.transitions)
#' results1$model.fit$Item.RMSEA
#'
#'############################################################################
#'# Example 2: TDCM with full measurement invariance and different Q-matrices
#'############################################################################
#'
#'# Load dataset: T=3, A=2
#'data(data.tdcm03, package = "TDCM")
#'data <- data.tdcm03$data
#'q1 <- data.tdcm03$q.matrix.1
#'q2 <- data.tdcm03$q.matrix.2
#'q3 <- data.tdcm03$q.matrix.3
#'q <- data.tdcm03$q.matrix.stacked
#'
#'# Estimate model
#'model2 <- TDCM::tdcm(data, q, num.time.points = 3,
#'                     rule = "LCDM",
#'                     num.q.matrix = 3,
#'                     num.items = c(10,10,10))
#'
#'#----------------------------------------------------------------------------
#'# With different transition options
#'#----------------------------------------------------------------------------
#'
#'## a) If transition.option = 1
#'
#'results2_option1 <- TDCM::tdcm.summary(model2)
#'results2_option1$transition.probabilities
#'#, , Attribute 1: Time 1 to Time 3
#'#
#'#       T3 [0] T3 [1]
#'#T1 [0]  0.202  0.798
#'#T1 [1]  0.146  0.854
#'#
#'#, , Attribute 2: Time 1 to Time 3
#'#
#'#       T3 [0] T3 [1]
#'#T1 [0]  0.325  0.675
#'#T1 [1]  0.257  0.743
#'
#'results2_option1$reliability
#'#            pt bis info gain polychor ave max tr P(t>.6) P(t>.7) P(t>.8)
#'#Attribute 1  0.550     0.387    0.737      0.830   0.888   0.780   0.643
#'#Attribute 2  0.665     0.474    0.808      0.851   0.899   0.801   0.694
#'
#'# b) If transition.option = 2
#'
#'# Summary with transition.option = 2
#'results2_option2 <- TDCM::tdcm.summary(model2, transition.option = 2)
#'results2_option2$transition.probabilities
#'#, , Attribute 1: Time 1 to Time 2
#'#
#'#.    [0]   [1]
#'#[0] 0.510 0.490
#'#[1] 0.424 0.576
#'#
#'#, , Attribute 2: Time 1 to Time 2
#'#
#'#     [0]   [1]
#'#[0] 0.456 0.544
#'#[1] 0.334 0.666
#'#
#'#, , Attribute 1: Time 1 to Time 3
#'#
#'#     [0]   [1]
#'#[0] 0.202 0.798
#'#[1] 0.146 0.854
#'#
#'#, , Attribute 2: Time 1 to Time 3
#'#
#'#     [0]   [1]
#'#[0] 0.325 0.675
#'#[1] 0.257 0.743
#'
#'results2_option2$reliability
#'#, , T1 to T2
#'#
#'#.           pt bis info gain polychor ave max tr P(t>.6) P(t>.7) P(t>.8)
#'#Attribute 1  0.586     0.444    0.770      0.796   0.828   0.710   0.581
#'#Attribute 2  0.692     0.503    0.838      0.853   0.885   0.799   0.713
#'#
#'#, , T1 to T3
#'#
#'#.           pt bis info gain polychor ave max tr P(t>.6) P(t>.7) P(t>.8)
#'#Attribute 1  0.550     0.387    0.737      0.830   0.888   0.780   0.643
#'#Attribute 2  0.665     0.474    0.808      0.851   0.899   0.801   0.694
#'
#'## c) If transition.option = 3
#'
#'results2_option3 <- TDCM::tdcm.summary(model2, transition.option = 3)
#'results2_option3$transition.probabilities
#'#, , Attribute 1: Time 1 to Time 2
#'#
#'#    [0]   [1]
#'#[0] 0.510 0.490
#'#[1] 0.424 0.576
#'#
#'#, , Attribute 2: Time 1 to Time 2
#'#
#'#     [0]   [1]
#'#[0] 0.456 0.544
#'#[1] 0.334 0.666
#'#
#'#, , Attribute 1: Time 2 to Time 3
#'#
#'#     [0]   [1]
#'#[0] 0.183 0.817
#'#[1] 0.188 0.812
#'#
#'#, , Attribute 2: Time 2 to Time 3
#'#
#'#     [0]   [1]
#'#[0] 0.361 0.639
#'#[1] 0.262 0.738
#'
#'results2_option3$reliability
#'#, , T1 to T2
#'#
#'#.           pt bis info gain polychor ave max tr P(t>.6) P(t>.7) P(t>.8)
#'#Attribute 1  0.586     0.444    0.770      0.796   0.828   0.710   0.581
#'#Attribute 2  0.692     0.503    0.838      0.853   0.885   0.799   0.713
#'#
#'#, , T2 to T3
#'#
#'#.            pt bis info gain polychor ave max tr P(t>.6) P(t>.7) P(t>.8)
#'#Attribute 1  0.537     0.396    0.724      0.801   0.841   0.724   0.578
#'#Attribute 2  0.691     0.502    0.861      0.853   0.880   0.799   0.714
#'
#'#----------------------------------------------------------------------------
#'# With different thresholds
#'#----------------------------------------------------------------------------
#'
#' ## a) If classthreshold = 0.5 (default)
#'
#'results2_1 <- TDCM::tdcm.summary(model2)
#'head(results2_1$posterior.probabilities)
#'#    T1A1  T1A2  T2A1  T2A2  T3A1  T3A2
#'# 1 0.068 0.882 0.961 0.967 1.000 1.000
#'# 2 0.001 0.010 0.845 0.749 0.070 0.402
#'# 3 0.005 0.683 0.816 0.395 0.987 0.133
#'# 4 0.007 0.988 0.996 0.998 0.993 0.997
#'# 5 0.001 0.000 0.205 0.019 0.999 0.814
#'# 6 0.011 0.001 0.630 0.004 0.900 0.077
#'
#'head(results2_1$classifications)
#'#   T1A1 T1A2 T2A1 T2A2 T3A1 T3A2
#'# 1    0    1    1    1    1    1
#'# 2    0    0    1    1    0    0
#'# 3    0    1    1    0    1    0
#'# 4    0    1    1    1    1    1
#'# 5    0    0    0    0    1    1
#'# 6    0    0    1    0    1    0
#'
#' ## b) If classthreshold = 0.7
#'
#'results2_2 <- TDCM::tdcm.summary(model2, classthreshold = 0.7)
#'head(results2_2$posterior.probabilities)
#'#    T1A1  T1A2  T2A1  T2A2  T3A1  T3A2
#'# 1 0.068 0.882 0.961 0.967 1.000 1.000
#'# 2 0.001 0.010 0.845 0.749 0.070 0.402
#'# 3 0.005 0.683 0.816 0.395 0.987 0.133
#'# 4 0.007 0.988 0.996 0.998 0.993 0.997
#'# 5 0.001 0.000 0.205 0.019 0.999 0.814
#'# 6 0.011 0.001 0.630 0.004 0.900 0.077
#'
#'head(results2_2$classifications)
#'#   T1A1 T1A2 T2A1 T2A2 T3A1 T3A2
#'# 1    0    1    1    1    1    1
#'# 2    0    0    1    1    0    0
#'# 3    0    0    1    0    1    0
#'# 4    0    1    1    1    1    1
#'# 5    0    0    0    0    1    1
#'# 6    0    0    0    0    1    0
#'
#' ## c) If classthreshold = 0.3
#'
#'results2_3 <- TDCM::tdcm.summary(model2, classthreshold = 0.3)
#'head(results2_3$posterior.probabilities)
#'#    T1A1  T1A2  T2A1  T2A2  T3A1  T3A2
#'# 1 0.068 0.882 0.961 0.967 1.000 1.000
#'# 2 0.001 0.010 0.845 0.749 0.070 0.402
#'# 3 0.005 0.683 0.816 0.395 0.987 0.133
#'# 4 0.007 0.988 0.996 0.998 0.993 0.997
#'# 5 0.001 0.000 0.205 0.019 0.999 0.814
#'# 6 0.011 0.001 0.630 0.004 0.900 0.077
#'
#'head(results2_3$classifications)
#'#   T1A1 T1A2 T2A1 T2A2 T3A1 T3A2
#'# 1    0    1    1    1    1    1
#'# 2    0    0    1    1    0    1
#'# 3    0    1    1    1    1    0
#'# 4    0    1    1    1    1    1
#'# 5    0    0    0    0    1    1
#'# 6    0    0    1    0    1    0
#' }
tdcm.summary <- function(model, transition.option = 1, classthreshold = .50,
                         attribute.names = c()) {

  num.time.points = model$numtimepoints
  num.items <- length(model$itemfit.rmsea) # total items
  items <- num.items / num.time.points # Items per time point
  num.atts <- ncol(model$attribute.patt.splitted) / num.time.points # Number of attribute measured

  # posterior probabilities
  postprobs <- matrix(NA, nrow = model$N, ncol = num.atts * num.time.points)
  postprobs <- round(model$pattern[, 6:(6 + num.atts * num.time.points - 1)], 3)
  colnames(postprobs) <- t(outer(c(paste("T", 1:num.time.points, sep = "")), c(paste("A", 1:num.atts, sep = "")),
                                 FUN = "paste0"))

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
    o1 <- summary.option1(model = model, num.atts = num.atts, num.time.points = num.time.points,
                          attribute.names = attribute.names)
    trans <- o1$trans
    growth <- o1$growth
    rel <- o1$rel
    growth.effects <- o1$growth.effects

  } else if (transition.option == 2) {
    o2 <- summary.option2(model = model, num.atts = num.atts, num.time.points = num.time.points,
                          attribute.names = attribute.names)
    trans <- o2$trans
    growth <- o2$growth
    rel <- o2$rel
    growth.effects <- o2$growth.effects


  } else {
    o3 <- summary.option3(model = model, num.atts = num.atts, num.time.points = num.time.points,
                          attribute.names = attribute.names)
    trans <- o3$trans
    growth <- o3$growth
    rel <- o3$rel
    growth.effects <- o3$growth.effects

  }


  ### Parameters Block ###
  param <- summary.param(model = model, num.atts = num.atts, num.items = num.items, num.time.points =
                           num.time.points)

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
    "growth.effects" = growth.effects,
    "posterior.probabilities" = postprobs, "classifications" = estclass, "reliability" = rel$metrics,
    "transition.posteriors" = rel$transposts, "most.likely.transitions" = rel$mostlikelytransitions,
    "option" = transition.option, "att.corr" = cor, "model.fit" = mf, "numgroups" = 1)

  return(newList)

}
