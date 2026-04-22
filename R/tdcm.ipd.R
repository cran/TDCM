#' Assessing item parameter drift (IPD) in the Transition Diagnostic Classification Model (TDCM)
#'
#' The `tdcm.ipd()` function assesses item parameter drift (IPD) in the TDCM (e.g., Madison & Bradshaw, 2018) by applying the Wald test
#' for differential item functioning (de la Torre, 2011; Hou, de la Torre & Nandakumar, 2014).
#' The _p_-values are also calculated by a Holm adjustment for multiple comparisons. In the case of two time
#' points, an effect size of item parameter drift (labeled as UA in the ipd.stats value)
#' is defined as the weighted absolute difference of item response functions.
#'
#' @param model A \code{tdcm} object returned from the \code{\link{tdcm}} function.
#'
#' @return A list with the following items:
#' \itemize{
#'    \item \code{$ipd.stats}: Data frame containing results of item-wise Wald tests.
#'
#'    \item \code{$coef}: Data frame containing item parameter estimates for each time point.
#'
#'    \item \code{$estimates}: List of \eqn{\lambda} vectors containing all item parameter estimates.
#'
#'    \item \code{$item.probs.time}: List with predicted item response probabilities for each time point.
#'
#'  }
#'
#' @references
#' de la Torre, J. (2011). The Generalized DINA model framework. _Psychometrika, 76_, 179–199.
#' <doi:10.1007/s11336-011-9207-7>.
#'
#' Hou, L., de la Torre, J., & Nandakumar, R. (2014). Differential item functioning assessment
#' in cognitive diagnostic modeling: Application of the Wald test to investigate DIF
#' in the DINA model. _Journal of Educational Measurement, 51_, 98-125. <doi:10.1111/jedm.12036>.
#'
#' Madison, M. J., & Bradshaw, L. (2018a). Assessing growth in a diagnostic classification model
#' framework. _Psychometrika_, **83**(4), 963-990.
#' <doi:10.1007/s11336-018-9638-5>.
#'
#'
#' @examples
#' \donttest{
#'
#' ############################################################################
#'# Example 1: TDCM with full measurement invariance
#'############################################################################
#'
#'# Load dataset: T=2, A=4
#'data(data.tdcm01, package = "TDCM")
#'data <- data.tdcm01$data
#'q.matrix <- data.tdcm01$q.matrix
#'# Estimate model
#'model1 <- TDCM::tdcm(data, q.matrix, num.time.points = 2, invariance = TRUE,
#'                     rule = "LCDM", num.q.matrix = 1)
#'
#'# Run IPD analysis
#'ipd = tdcm.ipd(model1)
#'ipd$ipd.stats
#'ipd$coef
#'ipd$parameters
#'
#' }
#' @export
tdcm.ipd <- function(
    model
) {

  #model specifications
  n <- model$N     #sample size
  num.time.points <- model$numtimepoints     #number of time points
  num.items <- length(model$itemfit.rmsea) / model$numtimepoints #number of items
  q <- model$qt.matrix   #q.matrix

  #create data
  responses <- model$data    #item responses
  colnames(responses) <- NULL
  names(responses) <- rep(c(paste("Item", 1:num.items, sep = "")), num.time.points)
  gr <- c(rep(1:num.time.points, each = n))    #groups for mg gdina estimation

  #format item response data in mg form
  datg <- responses[, 1:num.items]
  for(i in 2:num.time.points){

    datg <- rbind(datg, responses[, (num.items * (i - 1) + 1):(i * num.items)])

  }

  #estimate mg model
  mg1 <- CDM::gdina(datg, q.matrix = q, group = gr, rule = model$rule,
                    linkfct = model$linkfct, method = "ML",
                    reduced.skillspace = F, progress = F)
  ipd0 <- CDM::gdina.dif(mg1)

  ipd1 <- ipd0$difstats
  ipd2 <- ipd0$coef
  names(ipd2) <- NULL
  #column names for estimates and SEs
  cn <- c()
  for(i in 1:num.time.points){

    cn <- append(cn, c(paste("Est_Time", i, sep = ""), paste("SE_Time", i, sep = "")))

  }
  names(ipd2) <- c("Link", "Item", "Item Number",
                   "Param Type", "Rule", "Est", "SE", "Attr", cn)

  ipd2$Link <- rep(model$linkfct, nrow(ipd2))
  ipd2$Rule <- rep("LCDM", nrow(ipd2))

  ipd3 <- ipd0$delta_all
  #ipd4 <- ipd0$varmat_all
  ipd5 <- ipd0$prob.exp.group
  names(ipd5) <- c(paste("Time", 1:num.time.points, sep = ""))

  ipd <- list("ipd.stats" = ipd1,
              "coef" = ipd2, "estimates" = ipd3,
              "item.probs.time" = ipd5
              )

  return(ipd)

} # tdcm.ipd
