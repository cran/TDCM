#' TDCM Condensation Rules
#'
#' A *condensation rule* is a formula that states how different attributes
#' combine to form an observed or latent response (Rupp, Templin, & Henson,
#' 2010). The \pkg{TDCM} package includes support for `"LCDM"`, `"DINA"`,
#' `"DINO"`, `"CRUM"`, `"RRUM"`, `"LCDM1"` for the LCDM with only main effects,
#' `"LCDM2"` for the LCDM with two-way interactions, `"LCDM3"`, and so on.
#' Evaluate `TDCM::tdcm.rules$TDCM` for a complete list of condensation rules
#' supported by the \pkg{TDCM} package.
#'
#' @references
#' Rupp, A. A., Templin, J., & Henson, R. (2010).
#' _Diagnostic Measurement: Theory, Methods, and Applications_. New York:
#' Guilford. ISBN: 9781606235430.
#'
#' @examples
#' TDCM::tdcm.rules$TDCM
#'
#' @order 1
#' @export
tdcm.rules <- {
  rules <- rbind.data.frame(
    cbind(  "LCDM",   "GDINA"),
    cbind(  "CRUM",    "ACDM"),
    cbind(  "DINA",    "DINA"),
    cbind(  "DINO",    "DINO"),
    cbind(  "RRUM",    "RRUM"),
    cbind( "LCDM1",  "GDINA1"),
    cbind( "LCDM2",  "GDINA2"),
    cbind( "LCDM3",  "GDINA3"),
    cbind( "LCDM4",  "GDINA4"),
    cbind( "LCDM5",  "GDINA5"),
    cbind( "LCDM6",  "GDINA6"),
    cbind( "LCDM7",  "GDINA7"),
    cbind( "LCDM8",  "GDINA8"),
    cbind( "LCDM9",  "GDINA9"),
    cbind("LCDM10", "GDINA10")
  ) # rules
  colnames(rules) <- c("TDCM", "CDM")
  rules # do NOT use return(rules), because tdcm.rules is NOT a function
} # tdcm.rules

#' Map a TDCM condensation rule to its corresponding CDM condensation rule.
#'
#' The `tdcm.rule.map` function is used internally by functions like
#' [TDCM::tdcm()], [TDCM::mg.tdcm()], `tdcm.rule.map()`, etc.
#'
#' @param tdcm.rule A string naming a TDCM condensation rule. The default value
#'     is `"LDCM"`.
#'
#' @return A string naming the [CDM::CDM-package] condensation rule that
#'     corresponds to supplied TDCM condensation rule.
#'
#' @keywords internal
#' @noRd
tdcm.rule.map <- function(tdcm.rule = "LDCM") {
  cdm.rule <- NULL
  if (tdcm_is_string(tdcm.rule, stop.when.false = TRUE)) {
    cdm.rule <- TDCM::tdcm.rules[tdcm.rules$TDCM == tdcm.rule, ]$CDM
  } # if
  is.na(cdm.rule) <- length(cdm.rule) == 0
  return(cdm.rule)
} # tdcm.rule.map

#' Utility function for for base estimation in TDCMs.
#'
#' The `tdcm.rule.as.cdm.rule` function is used internally by functions like
#' [TDCM::tdcm()], [TDCM::mg.tdcm()], etc.
#'
#' @param rule A string or a vector of itemwise condensation rules that describe
#'     a specific DCM to estimate when using functions like [TDCM::tdcm],
#'     [CDM::gdina()], etc. If `rule` is a vector, then the condensation rules
#'     can be specified itemwise. The default is `"LDCM"` for all items.
#'
#' @return An object of class `gdina` returned by the internal call to
#'     [CDM::gdina()].
#'
#' @keywords internal
#' @noRd
tdcm.rule.as.cdm.rule <- function(tdcm.rule) {
  cdm.rule <- sapply(
    tdcm.rule,
    tdcm.rule.map,
    simplify = TRUE,
    USE.NAMES = FALSE
  ) # cdm.rule
  return(cdm.rule)
} # tdcm.rule.as.cdm.rule
