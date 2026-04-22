#' The Transition Diagnostic Classification Model Framework
#'
#' The _R_ package \pkg{TDCM} enables longitudinal diagnostic classification model (DCM) analyses
#' within the transition diagnostic classification model (TDCM) framework.
#'
#' @details Diagnostic classification models (DCMs) are psychometric models designed to classify
#' examinees according to their proficiency on a set of categorical latent traits, known as
#' attributes. Longitudinal DCMs have been developed as psychometric options for modeling changes
#' in attribute proficiency over time.
#'
#' The \pkg{TDCM} implements estimation of the transition DCM (TDCM; Madison & Bradshaw, 2018a),
#' a longitudinal extension of the log-linear cognitive diagnosis model (LCDM; Henson, Templin, &
#' Willse, 2009). As the LCDM subsumes many other DCMs, many other DCMs can be estimated
#' longitudinally via the TDCM. The package includes functions to estimate the single-group and
#' multigroup TDCM, summarize results of interest including item parameters, growth proportions,
#' transition probabilities, transitional reliability, attribute correlations, model fit, and
#' growth plots.
#'
#' For more details and examples, see \code{vignette("TDCM", package = "TDCM")}.
#'
#' @author
#' The \pkg{TDCM} package was created by:
#' * Matthew J. Madison, University of Georgia, mjmadison@uga.edu
#' * Sergio Haab, University of Iowa, sergio-haab@uiowa.edu
#' * Constanza Mardones-Segovia, University of Georgia, constanza.mardonessegovia@uga.edu
#' * Minjeong Jeon, University of California - Los Angeles, mjjeon@ucla.edu
#' * Michael E. Cotterell, University of Georgia, mepcott@uga.edu
#'
#' More information about the \pkg{TDCM} package authors, copyright holders, funders, etc. is
#' available in the package's `DESCRIPTION` file. You can see this information by executing
#' `utils::packageDescription("TDCM")`.
#'
#' @references
#' de la Torre, J. (2011). The Generalized DINA model framework. _Psychometrika, 76_, 179–199.
#' <doi:10.1007/s11336-011-9207-7>.
#'
#' George, A. C., Robitzsch, A., Kiefer, T., Gross, J., & Ünlü , A. (2016). The _R_ package **CDM**
#' for cognitive diagnosis models. _Journal of Statistical Software,74_(2), 1-24.
#' <doi:10.18637/jss.v074.i02>
#'
#' Henson, R., Templin, J., & Willse, J. (2009). Defining a family of cognitive diagnosis models
#' using log-linear models with latent variables. _Psychometrika, 74_, 191-21.
#' <doi:10.1007/s11336-008-9089-5>.
#'
#' Johnson, M. S., & Sinharay, S. (2020). The reliability of the posterior probability of skill
#' attainment in diagnostic classification models. _Journal of Educational Measurement, 47_(1),
#' 5–31. <doi:10.3102/1076998619864550>.
#'
#' Kaya, Y., & Leite, W. (2017). Assessing change in latent skills across time With longitudinal
#' cognitive diagnosis modeling: An evaluation of model performance.
#' _Educational and Psychological Measurement, 77_(3), 369–388.
#' <doi:10.1177/0013164416659314>.
#'
#' Li, F., Cohen, A., Bottge, B., & Templin, J. (2015). A latent transition analysis model for
#' assessing change in cognitive skills. _Educational and Psychological Measurement, 76_(2),
#' 181–204.
#' <doi:10.1177/0013164415588946>.
#'
#' Madison, M. J. (2019). Reliably assessing growth with longitudinal diagnostic classification
#' models. _Educational Measurement: Issues and Practice, 38_(2), 68-78.
#' <doi:10.1111/emip.12243>.
#'
#' Madison, M. J., & Bradshaw, L. (2018a). Assessing growth in a diagnostic classification model
#' framework. _Psychometrika, 83_(4), 963-990.
#' <doi:10.1007/s11336-018-9638-5>.
#'
#' Madison, M. J., & Bradshaw, L. (2018b). Evaluating intervention effects in a diagnostic
#' classification model framework. _Journal of Educational Measurement, 55_(1), 32-51.
#' <doi:10.1111/jedm.12162>.
#'
#' Madison, M.J., Chung, S., Kim, J., & Bradshaw, L.P. (2024) Approaches to estimating longitudinal
#' diagnostic classification models. _Behaviormetrika, 51_(7), 7-19.
#' <doi:10.1007/s41237-023-00202-5>.
#'
#' Madison, M. J., Jeon, M., Cotterell, M., Haab, S., & Zor, S. (2025). TDCM: An R package for estimating
#' longitudinal diagnostic classification models. _Multivariate Behavioral Research_, 1–10.
#' <doi:10.1080/00273171.2025.2453454>.
#'
#' Ravand, H., Effatpanah, F., Kunina-Habenicht, O., & Madison, M. J. (2025). A didactic illustration
#' of writing skill growth through a longitudinal diagnostic classification model. _Frontiers in Psychology_.
#' <doi:10.3389/fpsyg.2024.1521808>.
#'
#' Rupp, A. A., Templin, J., & Henson, R. (2010).
#' _Diagnostic Measurement: Theory, Methods, and Applications_. New York: Guilford.
#' ISBN: 9781606235430.
#'
#' Schellman, M., & Madison, M. J. (2024). Estimating the reliability of skill transition in longitudinal DCMs.
#'  _Journal of Educational and Behavioral Statistics_.<doi:10.3102/10769986241256032>.
#'
#' Templin, J., & Bradshaw, L. (2013). Measuring the reliability of diagnostic classification model
#' examinee estimates. _Journal of Classification, 30_, 251-275.
#' <doi:10.1007/s00357-013-9129-4>.
#'
#' Wang. S., Yang. Y., Culpepper, S. A., & Douglas, J. (2018). Tracking skill acquisition with
#' cognitive diagnosis models: A higher-order, hidden Markov model With covariates.
#' _Journal of Educational and Behavioral Statistics, 43_(1), 57-87.
#' <doi:10.3102/1076998617719727>.
#'
#' @keywords internal
#' @aliases TDCM-package NULL
#' @docType package
#' @rdname TDCM-package
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
