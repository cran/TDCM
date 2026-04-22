#' @title One-parameter log-linear cognitive diagnosis model.
#'
#' @description Function to estimate the 1-PLCDM (Madison et al., 2023; Maas et al., 2023).
#'
#' @param data a required \eqn{N \times I} matrix. Binary item responses are in the columns.
#'
#' @param q.matrix a required \eqn{I \times A} matrix indicating which items measure which attributes.
#'
#' @param progress An optional logical indicating whether the function should print the progress of estimation.
#'
#' @note
#' Currently, this model cannot be embedded within the TDCM via the \code{rule} argument.
#'
#' @details Estimates the single-attribute and multi-attribute 1-PLCDM described in Madison et al. (2024).
#' Example shows that attribute subscores are sufficient statistics for classifications.
#'
#' @return An object of class \code{gdina} with entries as indicated in the CDM package.
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Example 1: A = 4
#' data(data.tdcm05)
#' dat5 <- data.tdcm05$data
#' qmat5 <- data.tdcm05$q.matrix
#'
#' # calibrate LCDM
#' m1 <- CDM::gdina(dat5, qmat5, linkfct = "logit", method = "ML")
#'
#' # calibrate 1-PLCDM
#' m2 <- TDCM::oneplcdm(dat5, qmat5)
#' summary(m2)

#' #demonstrate 1-PLCDM sum score sufficiency for each attribute
#' subscores <- cbind(rowSums(dat5[, 1:5]), rowSums(dat5[, 6:10]),
#' rowSums(dat5[, 11:15]), rowSums(dat5[, 16:20]))
#' colnames(subscores) <- c("Att1", "Att2", "Att3", "Att4")
#' proficiency <- cbind(m2$pattern[, 6] > .50, m2$pattern[, 7] > .50,
#' m2$pattern[, 8] > .50, m2$pattern[, 9] > .5) * 1
#' table(subscores[, 1], proficiency[, 1])
#' table(subscores[, 2], proficiency[, 2])
#' table(subscores[, 3], proficiency[, 3])
#' table(subscores[, 4], proficiency[, 4])
#'
#' #plot sum score sufficiency for each attribute
#' posterior1pl <- m2$pattern[, 6:9]
#' posteriorlcdm <- m1$pattern[, 6:9]

#' oldpar <- par(mfrow = c(2, 2))
#' for (i in 1:4) {
#'  plot(subscores[, i], posteriorlcdm[, i], pch = 19,las = 1, cex.lab = 1.5,
#'  xlab = "Sum Scores", ylab = "P(proficiency)",
#'  cex.main = 1.5, col = "grey", xaxt = "n", yaxt = "n", cex = 1.2,
#'  main = paste("Attribute ", i, sep = ""))
#'  graphics::axis(side = 1, at = c(0, 1, 2, 3, 4, 5), )
#'  graphics::axis(side = 2, at = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0), las = 1)
#'  graphics::points(subscores[, i], posterior1pl[, i], col = "black", pch = 18, cex = 1.5)
#'  graphics::abline(a = .50, b = 0, col = "red")
#'  graphics::legend("bottomright", c("1-PLCDM", "LCDM"), col = c("black", "grey"),
#'  pch = c(18 ,19), box.lwd = 0, box.col = "white", bty = 'n')
#' }
#' par(oldpar)
#' }
#'
#' @references
#'
#' George, A. C., Robitzsch, A., Kiefer, T., Gross, J., & Ünlü , A. (2016). The R package CDM for
#' cognitive diagnosis models. \emph{Journal of Statistical Software, 74}(2), 1-24.
#'
#' Henson, R., Templin, J., & Willse, J. (2009). Defining a family of cognitive diagnosis models
#' using log linear models with latent variables. \emph{Psychometrika, 74}, 191-21.
#'
#' Madison, M.J., Wind, S., Maas, L., Yamaguchi, K. & Haab, S. (2024). A one-parameter diagnostic
#' classification model with familiar measurement properties. \emph{Journal of Educational Measurement}.
#'
#' Maas, L., Madison, M. J., & Brinkhuis, M. J. (2024). Properties and performance of the one-parameter
#' log-linear cognitive diagnosis model. \emph{Frontiers}.


oneplcdm <- function(data, q.matrix, progress = TRUE) { # open function

  # check q.matrix simple
  s <- rowSums(q.matrix)
  if (sum(s) == nrow(q.matrix)) { # open if q.matrix is simple

    # print line
    if (progress == TRUE) {
      print("Estimating 1-PLCDM...", quote = FALSE)
    }

    # estimate full lcdm
    m1 <- CDM::gdina(data, q.matrix, linkfct = "logit", method = "ML", progress = FALSE, maxit = 1)

    # how many items and attributes
    I <- nrow(q.matrix)
    num.atts <- ncol(q.matrix)

    if (num.atts == 1) { # open single attribute case
      c0 <- m1$coef
      dd <- diag(nrow(c0))
      dd[seq(4, 2 * I, by = 2), 2] <- 1
      dd[, seq(4, 2 * I, by = 2)] <- 0
      dd <- dd[, -seq(4, 2 * I, by = 2)]
      m2 <- CDM::gdina(data, q.matrix, linkfct = "logit", method = "ML", progress = FALSE, delta.designmatrix = dd)

      if (progress == TRUE) {
        print("Estimation is complete. Use the CDM summary function to display results.", quote = FALSE)
      }
    } # end single attribute

    if (num.atts > 1) { # open multiattribute case
      c0 <- m1$coef
      delta.designmatrix <- matrix(0, nrow = nrow(c0), ncol = nrow(c0))
      diag(delta.designmatrix) <- 1
      all <- c()
      for (i in 1:num.atts) {
        x <- (which(q.matrix[, i] == 1)) * 2
        y <- x[1]
        x <- x[-1]
        delta.designmatrix[x, y] <- 1
        delta.designmatrix[, x] <- 0
        all <- append(all, x)
      }
      all <- sort(all)
      delta.designmatrix <- delta.designmatrix[, -all]

      m2 <- CDM::gdina(data, q.matrix,
                  linkfct = "logit", method = "ML",
                  delta.designmatrix = delta.designmatrix, HOGDINA = 0,
                  progress = FALSE)

      if (progress == TRUE) {
        print("Estimation is complete. Use the CDM summary function to display results.", quote = FALSE)
      }
    } # end single attribute

    return(m2)
  } # close if not simple Q-matrix

  else {
    stop("Q-matrix has complex items. The 1-PLCDM can only be employed for single attribute assessment or
         for a multi-attribute assessment with a simple structure Q-matrix.")
  }

} # close function
