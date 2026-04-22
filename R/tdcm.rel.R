#' Utility function to compute transition reliability.
#'
#' @details
#' Includes longitudinal DCM reliability metrics developed by Schellman and Madison (2021).
#'
#' @param model gdina object from tdcm estimation
#' @param num.atts number of attributes
#' @param num.time.points number of time points
#' @param transition.option transition.option specified in summary function
#' @param attribute.names optional attribute names
#'
#' @return Several reliability metrics.
#'
#' @references
#' Schellman, M., & Madison, M. J. (2024). Estimating the reliability of skill transition in longitudinal DCMs.
#' \emph{Journal of Educational and Behavioral Statistics}.
#'
#' @keywords internal
#' @noRd
tdcm.rel <- function(model, num.atts, num.time.points, transition.option, attribute.names = c()) {
  # sample size
  if (model$G == 1) {
    N <- model$N
  } else {
    N <- sum(model$N)
  }

  ################################################################
  # transition.option = 1, first to last
  if (transition.option == 1) {
    A <- num.atts # number of attributes

    # extract posteriors and profile proportions
    posteriors <- model$posterior
    if (model$G == 1) { # single group
      est_baserates <- data.frame(model$attr.prob)
      transposts <- array(NA, dim = c(N, 4, A))
      transclass <- array(NA, dim = c(N, 4, A))
    } else { # multiple group
      est_baserates <- data.frame(rowSums(model$N * model$attr.prob) / sum(model$N))
      transposts <- array(NA, dim = c(sum(N), 4, A))
      transclass <- array(NA, dim = c(sum(N), 4, A))
    }
    # enumerate all profiles, reverse order from mplus
    classes <- data.frame(model$attribute.patt.splitted)

    # store results for each metric and each attribute
    relresults <- array(NA, dim = c(A, 12, 1))
    baserates <- matrix(NA, nrow = 1, ncol = 4)

    for (k in 1:A) {
      # rows for the transitions
      tr00 <- which(classes[, k] == 0 & classes[, k + (num.time.points - 1) * (A)] == 0)
      tr01 <- which(classes[, k] == 0 & classes[, k + (num.time.points - 1) * (A)] == 1)
      tr10 <- which(classes[, k] == 1 & classes[, k + (num.time.points - 1) * (A)] == 0)
      tr11 <- which(classes[, k] == 1 & classes[, k + (num.time.points - 1) * (A)] == 1)

      # base rates for the transitions
      br00 <- sum(est_baserates[tr00, 1])
      br01 <- sum(est_baserates[tr01, 1])
      br10 <- sum(est_baserates[tr10, 1])
      br11 <- sum(est_baserates[tr11, 1])

      baserates[1, 1] <- br00
      baserates[1, 2] <- br01
      baserates[1, 3] <- br10
      baserates[1, 4] <- br11

      # posteriors for the transitions
      if (num.atts == 1 & num.time.points == 2) {
        post00 <- posteriors[, tr00]
        post01 <- posteriors[, tr01]
        post10 <- posteriors[, tr10]
        post11 <- posteriors[, tr11]
      } else {
        post00 <- rowSums(posteriors[, tr00])
        post01 <- rowSums(posteriors[, tr01])
        post10 <- rowSums(posteriors[, tr10])
        post11 <- rowSums(posteriors[, tr11])
      }


      #########################################################################
      # average most likely transition (Madison, 2019)
      alltransposts <- data.frame(cbind(post00, post01, post10, post11))
      transposts[, , k] <- round(as.matrix(alltransposts), 3)
      transclass[, , k] <- apply(alltransposts, 1, which.max)
      transclass[, , k][transclass[, , k] == 1] <- "    00"
      transclass[, , k][transclass[, , k] == 2] <- "    01"
      transclass[, , k][transclass[, , k] == 3] <- "    10"
      transclass[, , k][transclass[, , k] == 4] <- "    11"

      pmax <- matrix(NA, nrow = nrow(alltransposts))
      for (j in 1:nrow(alltransposts)) {
        pmax[j] <- max(alltransposts[j, ])
      }
      relresults[k, 5, 1] <- mean(pmax)
      relresults[k, 6, 1] <- mean(pmax > .6)
      relresults[k, 7, 1] <- mean(pmax > .7)
      relresults[k, 8, 1] <- mean(pmax > .8)
      relresults[k, 9, 1] <- mean(pmax > .9)

      #########################################################################
      # Longitudinal point biserial metric

      pb00 <- mean(post00^2 - br00^2) / (br00 * (1 - br00))
      pb01 <- mean(post01^2 - br01^2) / (br01 * (1 - br01))
      pb10 <- mean(post10^2 - br10^2) / (br10 * (1 - br10))
      pb11 <- mean(post11^2 - br11^2) / (br11 * (1 - br11))

      if (is.na(pb00) | pb00 == Inf) {
        pb <- (pb01 + pb10 + pb11) / 3
        relresults[k, 1, 1] <- pb
        pb_w <- pb01 * br01 + pb10 * br10 + pb11 * br11
        relresults[k, 10, 1] <- pb_w
      } else if (is.na(pb01) | pb01 == Inf) {
        pb <- (pb00 + pb10 + pb11) / 3
        relresults[k, 1, 1] <- pb
        pb_w <- pb00 * br00 + pb10 * br10 + pb11 * br11
        relresults[k, 10, 1] <- pb_w
      } else if (is.na(pb10) | pb10 == Inf) {
        pb <- (pb00 + pb01 + pb11) / 3
        relresults[k, 1, 1] <- pb
        pb_w <- pb00 * br00 + pb01 * br01 + pb11 * br11
        relresults[k, 10, 1] <- pb_w
      } else if (is.na(pb11) | pb11 == Inf) {
        pb <- (pb00 + pb01 + pb10) / 3
        relresults[k, 1, 1] <- pb
        pb_w <- pb00 * br00 + pb01 * br01 + pb10 * br10
        relresults[k, 10, 1] <- pb_w
      } else {
        pb <- (pb00 + pb01 + pb10 + pb11) / 4
        relresults[k, 1, 1] <- pb
        pb_w <- pb00 * br00 + pb01 * br01 + pb10 * br10 + pb11 * br11
        relresults[k, 10, 1] <- pb_w
      }

      if (is.na(pb00)) {
        pb00 <- 0
      }
      if (is.na(pb01)) {
        pb01 <- 0
      }
      if (is.na(pb10)) {
        pb10 <- 0
      }
      if (is.na(pb11)) {
        pb11 <- 0
      }

      #########################################################################
      # Longitudinal parallel forms metric

      est_baserates2 <- est_baserates
      est_baserates2[est_baserates2 == 0] <- 0.00001

      # updated base rates for the transitions
      br00_2 <- sum(est_baserates2[tr00, 1])
      br01_2 <- sum(est_baserates2[tr01, 1])
      br10_2 <- sum(est_baserates2[tr10, 1])
      br11_2 <- sum(est_baserates2[tr11, 1])

      num00 <- NULL
      num01 <- NULL
      num10 <- NULL
      num11 <- NULL

      for (i in 1:nrow(classes)) {
        num00_comp <- (1 / est_baserates2[i, 1]) * (mean(post00 * posteriors[, i]))^2
        num00 <- c(num00, num00_comp)

        num01_comp <- (1 / est_baserates2[i, 1]) * (mean(post01 * posteriors[, i]))^2
        num01 <- c(num01, num01_comp)

        num10_comp <- (1 / est_baserates2[i, 1]) * (mean(post10 * posteriors[, i]))^2
        num10 <- c(num10, num10_comp)

        num11_comp <- (1 / est_baserates2[i, 1]) * (mean(post11 * posteriors[, i]))^2
        num11 <- c(num11, num11_comp)
      }

      num00 <- sum(num00) - br00_2^2
      denom00 <- mean((post00)^2 - br00_2^2)
      pf00 <- num00 / denom00

      num01 <- sum(num01) - br01_2^2
      denom01 <- mean((post01)^2 - br01_2^2)
      pf01 <- num01 / denom01

      num10 <- sum(num10) - br10_2^2
      denom10 <- mean((post10)^2 - br10_2^2)
      pf10 <- num10 / denom10

      num11 <- sum(num11) - br11_2^2
      denom11 <- mean((post11)^2 - br11_2^2)
      pf11 <- num11 / denom11

      pf <- (pf00 + pf01 + pf10 + pf11) / 4
      relresults[k, 2, 1] <- pf

      pf_w <- pf00 * br00_2 + pf01 * br01_2 + pf10 * br10_2 + pf11 * br11_2
      relresults[k, 11, 1] <- pf_w

      #########################################################################
      # Longitudinal information gain metric

      post00_2 <- post00
      post01_2 <- post01
      post10_2 <- post10
      post11_2 <- post11

      post00_2[post00_2 == 0] <- .0001
      post01_2[post01_2 == 0] <- .0001
      post10_2[post10_2 == 0] <- .0001
      post11_2[post11_2 == 0] <- .0001

      post00_2[post00_2 > .99999] <- .9999
      post01_2[post01_2 > .99999] <- .9999
      post10_2[post10_2 > .99999] <- .9999
      post11_2[post11_2 > .99999] <- .9999

      H100 <- -br00_2 * log(br00_2) - (1 - br00_2) * log(1 - br00_2)
      H101 <- -br01_2 * log(br01_2) - (1 - br01_2) * log(1 - br01_2)
      H110 <- -br10_2 * log(br10_2) - (1 - br10_2) * log(1 - br10_2)
      H111 <- -br11_2 * log(br11_2) - (1 - br11_2) * log(1 - br11_2)

      H200 <- -mean(post00_2 * log(post00_2) + (1 - post00_2) * log(1 - post00_2))
      H201 <- -mean(post01_2 * log(post01_2) + (1 - post01_2) * log(1 - post01_2))
      H210 <- -mean(post10_2 * log(post10_2) + (1 - post10_2) * log(1 - post10_2))
      H211 <- -mean(post11_2 * log(post11_2) + (1 - post11_2) * log(1 - post11_2))

      ig <- ((1 - exp(-2 * (H100 - H200))) + (1 - exp(-2 * (H101 - H201))) + (1 - exp(-2 * (H110 - H210))) + (1 - exp(-2 * (H111 - H211)))) / 4
      relresults[k, 3, 1] <- ig

      ig_w <- (1 - exp(-2 * (H100 - H200))) * br00 + (1 - exp(-2 * (H101 - H201))) * br01 + (1 - exp(-2 * (H110 - H210))) * br10 + (1 - exp(-2 * (H111 - H211))) * br11
      relresults[k, 12, 1] <- ig_w

      #########################################################################
      # Madison, 2019 metric

      # create 4x4 matrix of re-test probabilities
      cell11 <- post00_2^2
      cell22 <- post01_2^2
      cell33 <- post10_2^2
      cell44 <- post11_2^2
      cell12 <- post00_2 * post01_2
      cell13 <- post00_2 * post10_2
      cell14 <- post00_2 * post11_2
      cell23 <- post01_2 * post10_2
      cell24 <- post01_2 * post11_2
      cell34 <- post10_2 * post11_2

      fourmat <- matrix(NA, nrow = 4, ncol = 4)

      fourmat[1, 1] <- mean(cell11)
      fourmat[2, 2] <- mean(cell22)
      fourmat[3, 3] <- mean(cell33)
      fourmat[4, 4] <- mean(cell44)
      fourmat[1, 2] <- fourmat[2, 1] <- mean(cell12)
      fourmat[1, 3] <- fourmat[3, 1] <- mean(cell13)
      fourmat[1, 4] <- fourmat[4, 1] <- mean(cell14)
      fourmat[2, 3] <- fourmat[3, 2] <- mean(cell23)
      fourmat[2, 4] <- fourmat[4, 2] <- mean(cell24)
      fourmat[3, 4] <- fourmat[4, 3] <- mean(cell34)

      relresults[k, 4, 1] <- suppressWarnings(polycor::polychor(fourmat))
    } # end attribute loop

    if (length(attribute.names) == num.atts) {
      rnames <- attribute.names
    } else {
      rnames <- c(paste("Attribute", 1:num.atts, sep = " "))
    }

    cnames <- c(
      "pt bis", "pa forms", "info gain", "polychor", "ave max tr",
      "P(t>.6)", "P(t>.7)", "P(t>.8)", "P(t>.9)", "wt pt bis", "wt pa forms", "wt info gain"
    )
    mnames <- c(paste("T1 to", paste("T", num.time.points, sep = ""), sep = " "))

    if (length(attribute.names) == num.atts) {
      tp_mnames <- paste(attribute.names,
                         c(paste(": T1 to", paste("T", num.time.points, sep = ""), sep = " ")),
                         sep = ""
      )
    } else {
      tp_mnames <- paste(c(paste("Attribute", 1:num.atts, sep = " ")),
                         c(paste(": T1 to", paste("T", num.time.points, sep = ""), sep = " ")),
                         sep = ""
      )
    }
    tp_rnames <- c(1:sum(N))
    tp_cnames <- c("00", "01", "10", "11")

    transclass2 <- as.matrix(noquote(transclass[, 1, ]))

    tc_rnames <- c(1:sum(N))
    if (length(attribute.names) == num.atts) {
      tc_cnames <- paste(attribute.names,
                         c(paste(": T1 to", paste("T", num.time.points, sep = ""), sep = " ")),
                         sep = ""
      )
    } else {
      tc_cnames <- paste(c(paste("Attribute", 1:num.atts, sep = " ")),
                         c(paste(": T1 to", paste("T", num.time.points, sep = ""), sep = " ")),
                         sep = ""
      )
    }
  } # transition option = 1 loop

  ####################################################################################
  ####################################################################################
  # transition.option = 2, first to each
  else if (transition.option == 2) {
    A <- num.atts # number of attributes

    # extract posteriors and profile proportions
    posteriors <- model$posterior
    if (model$G == 1) {
      est_baserates <- data.frame(model$attr.prob)
      transposts <- array(NA, dim = c(N, 4, A * (num.time.points - 1)))
      transclass <- array(NA, dim = c(N, 4, A * (num.time.points - 1)))
    } else {
      est_baserates <- data.frame(rowSums(model$N * model$attr.prob) / sum(model$N))
      transposts <- array(NA, dim = c(sum(N), 4, A * (num.time.points - 1)))
      transclass <- array(NA, dim = c(sum(N), 4, A * (num.time.points - 1)))
    }

    # enumerate all profiles, reverse order from mplus
    classes <- data.frame(model$attribute.patt.splitted)

    # store results for each metric and each attribute
    relresults <- array(NA, dim = c(A, 12, num.time.points - 1))
    baserates <- matrix(NA, nrow = 1, ncol = 4)


    # loop over each transition comparison, c is last dimension of array (row-column-matrix)
    for (c in 1:(num.time.points - 1)) {
      # loop over each attribute
      for (k in 1:A) {
        # rows for the transitions
        tr00 <- which(classes[, k] == 0 & classes[, k + c * A] == 0)
        tr01 <- which(classes[, k] == 0 & classes[, k + c * A] == 1)
        tr10 <- which(classes[, k] == 1 & classes[, k + c * A] == 0)
        tr11 <- which(classes[, k] == 1 & classes[, k + c * A] == 1)

        # base rates for the transitions
        br00 <- sum(est_baserates[tr00, 1])
        br01 <- sum(est_baserates[tr01, 1])
        br10 <- sum(est_baserates[tr10, 1])
        br11 <- sum(est_baserates[tr11, 1])

        baserates[1, 1] <- br00
        baserates[1, 2] <- br01
        baserates[1, 3] <- br10
        baserates[1, 4] <- br11

        # posteriors for the transitions
        if (num.atts == 1 & num.time.points == 2) {
          post00 <- posteriors[, tr00]
          post01 <- posteriors[, tr01]
          post10 <- posteriors[, tr10]
          post11 <- posteriors[, tr11]
        } else {
          post00 <- rowSums(posteriors[, tr00])
          post01 <- rowSums(posteriors[, tr01])
          post10 <- rowSums(posteriors[, tr10])
          post11 <- rowSums(posteriors[, tr11])
        }

        #########################################################################
        # average most likely transition (Madison, 2019)
        alltransposts <- data.frame(cbind(post00, post01, post10, post11))
        transposts[, , k + (c - 1) * num.atts] <- round(as.matrix(alltransposts), 3)
        transclass[, , k + (c - 1) * num.atts] <- apply(alltransposts, 1, which.max)
        transclass[, , k + (c - 1) * num.atts][transclass[, , k + (c - 1) * num.atts] == 1] <- "    00"
        transclass[, , k + (c - 1) * num.atts][transclass[, , k + (c - 1) * num.atts] == 2] <- "    01"
        transclass[, , k + (c - 1) * num.atts][transclass[, , k + (c - 1) * num.atts] == 3] <- "    10"
        transclass[, , k + (c - 1) * num.atts][transclass[, , k + (c - 1) * num.atts] == 4] <- "    11"

        pmax <- matrix(NA, nrow = nrow(alltransposts))
        for (j in 1:nrow(alltransposts)) {
          pmax[j] <- max(alltransposts[j, ])
        }
        relresults[k, 5, c] <- mean(pmax)
        relresults[k, 6, c] <- mean(pmax > .6)
        relresults[k, 7, c] <- mean(pmax > .7)
        relresults[k, 8, c] <- mean(pmax > .8)
        relresults[k, 9, c] <- mean(pmax > .9)

        #########################################################################
        # Longitudinal point biserial metric

        pb00 <- mean(post00^2 - br00^2) / (br00 * (1 - br00))
        pb01 <- mean(post01^2 - br01^2) / (br01 * (1 - br01))
        pb10 <- mean(post10^2 - br10^2) / (br10 * (1 - br10))
        pb11 <- mean(post11^2 - br11^2) / (br11 * (1 - br11))

        if (is.na(pb00) | pb00 == Inf) {
          pb <- (pb01 + pb10 + pb11) / 3
          relresults[k, 1, c] <- pb
          pb_w <- pb01 * br01 + pb10 * br10 + pb11 * br11
          relresults[k, 10, c] <- pb_w
        } else if (is.na(pb01) | pb01 == Inf) {
          pb <- (pb00 + pb10 + pb11) / 3
          relresults[k, 1, c] <- pb
          pb_w <- pb00 * br00 + pb10 * br10 + pb11 * br11
          relresults[k, 10, c] <- pb_w
        } else if (is.na(pb10) | pb10 == Inf) {
          pb <- (pb00 + pb01 + pb11) / 3
          relresults[k, 1, c] <- pb
          pb_w <- pb00 * br00 + pb01 * br01 + pb11 * br11
          relresults[k, 10, c] <- pb_w
        } else if (is.na(pb11) | pb11 == Inf) {
          pb <- (pb00 + pb01 + pb10) / 3
          relresults[k, 1, c] <- pb
          pb_w <- pb00 * br00 + pb01 * br01 + pb10 * br10
          relresults[k, 10, c] <- pb_w
        } else {
          pb <- (pb00 + pb01 + pb10 + pb11) / 4
          relresults[k, 1, c] <- pb
          pb_w <- pb00 * br00 + pb01 * br01 + pb10 * br10 + pb11 * br11
          relresults[k, 10, c] <- pb_w
        }

        if (is.na(pb00)) {
          pb00 <- 0
        }
        if (is.na(pb01)) {
          pb01 <- 0
        }
        if (is.na(pb10)) {
          pb10 <- 0
        }
        if (is.na(pb11)) {
          pb11 <- 0
        }

        #########################################################################
        # Longitudinal parallel forms metric

        est_baserates2 <- est_baserates
        est_baserates2[est_baserates2 == 0] <- 0.00001

        # updated base rates for the transitions
        br00_2 <- sum(est_baserates2[tr00, 1])
        br01_2 <- sum(est_baserates2[tr01, 1])
        br10_2 <- sum(est_baserates2[tr10, 1])
        br11_2 <- sum(est_baserates2[tr11, 1])

        num00 <- NULL
        num01 <- NULL
        num10 <- NULL
        num11 <- NULL

        for (i in 1:nrow(classes)) {
          num00_comp <- (1 / est_baserates2[i, 1]) * (mean(post00 * posteriors[, i]))^2
          num00 <- c(num00, num00_comp)

          num01_comp <- (1 / est_baserates2[i, 1]) * (mean(post01 * posteriors[, i]))^2
          num01 <- c(num01, num01_comp)

          num10_comp <- (1 / est_baserates2[i, 1]) * (mean(post10 * posteriors[, i]))^2
          num10 <- c(num10, num10_comp)

          num11_comp <- (1 / est_baserates2[i, 1]) * (mean(post11 * posteriors[, i]))^2
          num11 <- c(num11, num11_comp)
        }

        num00 <- sum(num00) - br00_2^2
        denom00 <- mean((post00)^2 - br00_2^2)
        pf00 <- num00 / denom00

        num01 <- sum(num01) - br01_2^2
        denom01 <- mean((post01)^2 - br01_2^2)
        pf01 <- num01 / denom01

        num10 <- sum(num10) - br10_2^2
        denom10 <- mean((post10)^2 - br10_2^2)
        pf10 <- num10 / denom10

        num11 <- sum(num11) - br11_2^2
        denom11 <- mean((post11)^2 - br11_2^2)
        pf11 <- num11 / denom11

        pf <- (pf00 + pf01 + pf10 + pf11) / 4
        relresults[k, 2, c] <- pf

        pf_w <- pf00 * br00_2 + pf01 * br01_2 + pf10 * br10_2 + pf11 * br11_2
        relresults[k, 11, c] <- pf_w

        #########################################################################
        # Longitudinal information gain metric

        post00_2 <- post00
        post01_2 <- post01
        post10_2 <- post10
        post11_2 <- post11

        post00_2[post00_2 == 0] <- .0001
        post01_2[post01_2 == 0] <- .0001
        post10_2[post10_2 == 0] <- .0001
        post11_2[post11_2 == 0] <- .0001

        post00_2[post00_2 > .99999] <- .9999
        post01_2[post01_2 > .99999] <- .9999
        post10_2[post10_2 > .99999] <- .9999
        post11_2[post11_2 > .99999] <- .9999

        H100 <- -br00_2 * log(br00_2) - (1 - br00_2) * log(1 - br00_2)
        H101 <- -br01_2 * log(br01_2) - (1 - br01_2) * log(1 - br01_2)
        H110 <- -br10_2 * log(br10_2) - (1 - br10_2) * log(1 - br10_2)
        H111 <- -br11_2 * log(br11_2) - (1 - br11_2) * log(1 - br11_2)

        H200 <- -mean(post00_2 * log(post00_2) + (1 - post00_2) * log(1 - post00_2))
        H201 <- -mean(post01_2 * log(post01_2) + (1 - post01_2) * log(1 - post01_2))
        H210 <- -mean(post10_2 * log(post10_2) + (1 - post10_2) * log(1 - post10_2))
        H211 <- -mean(post11_2 * log(post11_2) + (1 - post11_2) * log(1 - post11_2))

        ig <- ((1 - exp(-2 * (H100 - H200))) + (1 - exp(-2 * (H101 - H201))) + (1 - exp(-2 * (H110 - H210))) + (1 - exp(-2 * (H111 - H211)))) / 4
        relresults[k, 3, c] <- ig

        ig_w <- (1 - exp(-2 * (H100 - H200))) * br00 + (1 - exp(-2 * (H101 - H201))) * br01 + (1 - exp(-2 * (H110 - H210))) * br10 + (1 - exp(-2 * (H111 - H211))) * br11
        relresults[k, 12, c] <- ig_w

        #########################################################################
        # Madison, 2019 metric

        # create 4x4 matrix of re-test probabilities
        cell11 <- post00_2^2
        cell22 <- post01_2^2
        cell33 <- post10_2^2
        cell44 <- post11_2^2
        cell12 <- post00_2 * post01_2
        cell13 <- post00_2 * post10_2
        cell14 <- post00_2 * post11_2
        cell23 <- post01_2 * post10_2
        cell24 <- post01_2 * post11_2
        cell34 <- post10_2 * post11_2

        fourmat <- matrix(NA, nrow = 4, ncol = 4)

        fourmat[1, 1] <- mean(cell11)
        fourmat[2, 2] <- mean(cell22)
        fourmat[3, 3] <- mean(cell33)
        fourmat[4, 4] <- mean(cell44)
        fourmat[1, 2] <- fourmat[2, 1] <- mean(cell12)
        fourmat[1, 3] <- fourmat[3, 1] <- mean(cell13)
        fourmat[1, 4] <- fourmat[4, 1] <- mean(cell14)
        fourmat[2, 3] <- fourmat[3, 2] <- mean(cell23)
        fourmat[2, 4] <- fourmat[4, 2] <- mean(cell24)
        fourmat[3, 4] <- fourmat[4, 3] <- mean(cell34)

        relresults[k, 4, c] <- suppressWarnings(polycor::polychor(fourmat))
      } # end attribute loop
    } # end transition comparison loop

    if (length(attribute.names) == num.atts) {
      rnames <- attribute.names
    } else {
      rnames <- c(paste("Attribute", 1:num.atts, sep = " "))
    }
    cnames <- c(
      "pt bis", "pa forms", "info gain", "polychor", "ave max tr",
      "P(t>.6)", "P(t>.7)", "P(t>.8)", "P(t>.9)", "wt pt bis", "wt pa forms", "wt info gain"
    )
    mnames <- outer(c("T1 to "), c(paste("T", 2:num.time.points, sep = "")), FUN = "paste0")
    dim(mnames) <- NULL

    tp_rnames <- c(1:sum(N))
    tp_cnames <- c("00", "01", "10", "11")


    transclass <- noquote(transclass)
    transclass2 <- matrix(NA, nrow = sum(N), ncol = A * (num.time.points - 1))
    for (m in 1:(A * (num.time.points - 1))) {
      transclass2[, m] <- transclass[, 1, m]
    }
    transclass2 <- noquote(transclass2)

    tc_rnames <- c(1:sum(N))
    if (length(attribute.names) == num.atts) {
      s1 <- attribute.names
    } else {
      s1 <- c(paste("Attribute", 1:num.atts, sep = " "))
    }
    s2 <- outer(c(": T1 to "), c(paste("T", 2:num.time.points, sep = "")), FUN = "paste0")

    tc_cnames <- apply(expand.grid(s1, s2), 1, function(x) paste0(x, collapse = ""))
    tp_mnames <- tc_cnames
  } # end transition option = 2 if statement

  ####################################################################################
  ####################################################################################
  # transition.option = 3, successive
  else {
    A <- num.atts # number of attributes

    # extract posteriors and profile proportions
    posteriors <- model$posterior
    if (model$G == 1) {
      est_baserates <- data.frame(model$attr.prob)
      transposts <- array(NA, dim = c(N, 4, A * (num.time.points - 1)))
      transclass <- array(NA, dim = c(N, 4, A * (num.time.points - 1)))
    } else {
      est_baserates <- data.frame(rowSums(model$N * model$attr.prob) / sum(model$N))
      transposts <- array(NA, dim = c(sum(N), 4, A * (num.time.points - 1)))
      transclass <- array(NA, dim = c(sum(N), 4, A * (num.time.points - 1)))
    }
    # enumerate all profiles, reverse order from mplus
    classes <- data.frame(model$attribute.patt.splitted)

    # store results for each metric and each attribute
    relresults <- array(NA, dim = c(A, 12, num.time.points - 1))
    baserates <- matrix(NA, nrow = 1, ncol = 4)


    mnames <- c()
    tp_mnames <- c()
    tc_cnames <- c()
    # loop over each transition comparison, c is last dimension of array (row-column-matrix)
    for (c in 1:(num.time.points - 1)) {
      temp.name <- paste("T", c, " to T", c + 1, sep = "")
      mnames <- c(mnames, temp.name) # Combines matrix names into list

      # loop over each attribute
      for (k in 1:A) {
        temp.name2 <- paste("Attribute ", k, ": T", c, " to T", c + 1, sep = "")
        tp_mnames <- c(tp_mnames, temp.name2)
        tc_cnames <- c(tc_cnames, temp.name2)

        # rows for the transitions
        tr00 <- which(classes[, k + (c - 1) * A] == 0 & classes[, k + c * A] == 0)
        tr01 <- which(classes[, k + (c - 1) * A] == 0 & classes[, k + c * A] == 1)
        tr10 <- which(classes[, k + (c - 1) * A] == 1 & classes[, k + c * A] == 0)
        tr11 <- which(classes[, k + (c - 1) * A] == 1 & classes[, k + c * A] == 1)

        # base rates for the transitions
        br00 <- sum(est_baserates[tr00, 1])
        br01 <- sum(est_baserates[tr01, 1])
        br10 <- sum(est_baserates[tr10, 1])
        br11 <- sum(est_baserates[tr11, 1])

        baserates[1, 1] <- br00
        baserates[1, 2] <- br01
        baserates[1, 3] <- br10
        baserates[1, 4] <- br11

        # posteriors for the transitions
        if (num.atts == 1 & num.time.points == 2) {
          post00 <- posteriors[, tr00]
          post01 <- posteriors[, tr01]
          post10 <- posteriors[, tr10]
          post11 <- posteriors[, tr11]
        } else {
          post00 <- rowSums(posteriors[, tr00])
          post01 <- rowSums(posteriors[, tr01])
          post10 <- rowSums(posteriors[, tr10])
          post11 <- rowSums(posteriors[, tr11])
        }

        #########################################################################
        # average most likely transition (Madison, 2019)
        alltransposts <- data.frame(cbind(post00, post01, post10, post11))
        transposts[, , k + (c - 1) * num.atts] <- round(as.matrix(alltransposts), 3)
        transclass[, , k + (c - 1) * num.atts] <- apply(alltransposts, 1, which.max)
        transclass[, , k + (c - 1) * num.atts][transclass[, , k + (c - 1) * num.atts] == 1] <- "    00"
        transclass[, , k + (c - 1) * num.atts][transclass[, , k + (c - 1) * num.atts] == 2] <- "    01"
        transclass[, , k + (c - 1) * num.atts][transclass[, , k + (c - 1) * num.atts] == 3] <- "    10"
        transclass[, , k + (c - 1) * num.atts][transclass[, , k + (c - 1) * num.atts] == 4] <- "    11"

        pmax <- matrix(NA, nrow = nrow(alltransposts))
        for (j in 1:nrow(alltransposts)) {
          pmax[j] <- max(alltransposts[j, ])
        }
        relresults[k, 5, c] <- mean(pmax)
        relresults[k, 6, c] <- mean(pmax > .6)
        relresults[k, 7, c] <- mean(pmax > .7)
        relresults[k, 8, c] <- mean(pmax > .8)
        relresults[k, 9, c] <- mean(pmax > .9)

        #########################################################################
        # Longitudinal point biserial metric

        pb00 <- mean(post00^2 - br00^2) / (br00 * (1 - br00))
        pb01 <- mean(post01^2 - br01^2) / (br01 * (1 - br01))
        pb10 <- mean(post10^2 - br10^2) / (br10 * (1 - br10))
        pb11 <- mean(post11^2 - br11^2) / (br11 * (1 - br11))

        if (is.na(pb00) | pb00 == Inf) {
          pb <- (pb01 + pb10 + pb11) / 3
          relresults[k, 1, c] <- pb
          pb_w <- pb01 * br01 + pb10 * br10 + pb11 * br11
          relresults[k, 10, c] <- pb_w
        } else if (is.na(pb01) | pb01 == Inf) {
          pb <- (pb00 + pb10 + pb11) / 3
          relresults[k, 1, c] <- pb
          pb_w <- pb00 * br00 + pb10 * br10 + pb11 * br11
          relresults[k, 10, c] <- pb_w
        } else if (is.na(pb10) | pb10 == Inf) {
          pb <- (pb00 + pb01 + pb11) / 3
          relresults[k, 1, c] <- pb
          pb_w <- pb00 * br00 + pb01 * br01 + pb11 * br11
          relresults[k, 10, c] <- pb_w
        } else if (is.na(pb11) | pb11 == Inf) {
          pb <- (pb00 + pb01 + pb10) / 3
          relresults[k, 1, c] <- pb
          pb_w <- pb00 * br00 + pb01 * br01 + pb10 * br10
          relresults[k, 10, c] <- pb_w
        } else {
          pb <- (pb00 + pb01 + pb10 + pb11) / 4
          relresults[k, 1, c] <- pb
          pb_w <- pb00 * br00 + pb01 * br01 + pb10 * br10 + pb11 * br11
          relresults[k, 10, c] <- pb_w
        }

        if (is.na(pb00)) {
          pb00 <- 0
        }
        if (is.na(pb01)) {
          pb01 <- 0
        }
        if (is.na(pb10)) {
          pb10 <- 0
        }
        if (is.na(pb11)) {
          pb11 <- 0
        }

        #########################################################################
        # Longitudinal parallel forms metric

        est_baserates2 <- est_baserates
        est_baserates2[est_baserates2 == 0] <- 0.00001

        # updated base rates for the transitions
        br00_2 <- sum(est_baserates2[tr00, 1])
        br01_2 <- sum(est_baserates2[tr01, 1])
        br10_2 <- sum(est_baserates2[tr10, 1])
        br11_2 <- sum(est_baserates2[tr11, 1])

        num00 <- NULL
        num01 <- NULL
        num10 <- NULL
        num11 <- NULL

        for (i in 1:nrow(classes)) {
          num00_comp <- (1 / est_baserates2[i, 1]) * (mean(post00 * posteriors[, i]))^2
          num00 <- c(num00, num00_comp)

          num01_comp <- (1 / est_baserates2[i, 1]) * (mean(post01 * posteriors[, i]))^2
          num01 <- c(num01, num01_comp)

          num10_comp <- (1 / est_baserates2[i, 1]) * (mean(post10 * posteriors[, i]))^2
          num10 <- c(num10, num10_comp)

          num11_comp <- (1 / est_baserates2[i, 1]) * (mean(post11 * posteriors[, i]))^2
          num11 <- c(num11, num11_comp)
        }

        num00 <- sum(num00) - br00_2^2
        denom00 <- mean((post00)^2 - br00_2^2)
        pf00 <- num00 / denom00

        num01 <- sum(num01) - br01_2^2
        denom01 <- mean((post01)^2 - br01_2^2)
        pf01 <- num01 / denom01

        num10 <- sum(num10) - br10_2^2
        denom10 <- mean((post10)^2 - br10_2^2)
        pf10 <- num10 / denom10

        num11 <- sum(num11) - br11_2^2
        denom11 <- mean((post11)^2 - br11_2^2)
        pf11 <- num11 / denom11

        pf <- (pf00 + pf01 + pf10 + pf11) / 4
        relresults[k, 2, c] <- pf

        pf_w <- pf00 * br00_2 + pf01 * br01_2 + pf10 * br10_2 + pf11 * br11_2
        relresults[k, 11, c] <- pf_w

        #########################################################################
        # Longitudinal information gain metric

        post00_2 <- post00
        post01_2 <- post01
        post10_2 <- post10
        post11_2 <- post11

        post00_2[post00_2 == 0] <- .0001
        post01_2[post01_2 == 0] <- .0001
        post10_2[post10_2 == 0] <- .0001
        post11_2[post11_2 == 0] <- .0001

        post00_2[post00_2 > .99999] <- .9999
        post01_2[post01_2 > .99999] <- .9999
        post10_2[post10_2 > .99999] <- .9999
        post11_2[post11_2 > .99999] <- .9999

        H100 <- -br00_2 * log(br00_2) - (1 - br00_2) * log(1 - br00_2)
        H101 <- -br01_2 * log(br01_2) - (1 - br01_2) * log(1 - br01_2)
        H110 <- -br10_2 * log(br10_2) - (1 - br10_2) * log(1 - br10_2)
        H111 <- -br11_2 * log(br11_2) - (1 - br11_2) * log(1 - br11_2)

        H200 <- -mean(post00_2 * log(post00_2) + (1 - post00_2) * log(1 - post00_2))
        H201 <- -mean(post01_2 * log(post01_2) + (1 - post01_2) * log(1 - post01_2))
        H210 <- -mean(post10_2 * log(post10_2) + (1 - post10_2) * log(1 - post10_2))
        H211 <- -mean(post11_2 * log(post11_2) + (1 - post11_2) * log(1 - post11_2))

        ig <- ((1 - exp(-2 * (H100 - H200))) + (1 - exp(-2 * (H101 - H201))) + (1 - exp(-2 * (H110 - H210))) + (1 - exp(-2 * (H111 - H211)))) / 4
        relresults[k, 3, c] <- ig

        ig_w <- (1 - exp(-2 * (H100 - H200))) * br00 + (1 - exp(-2 * (H101 - H201))) * br01 + (1 - exp(-2 * (H110 - H210))) * br10 + (1 - exp(-2 * (H111 - H211))) * br11
        relresults[k, 12, c] <- ig_w

        #########################################################################
        # Madison, 2019 metric

        # create 4x4 matrix of re-test probabilities
        cell11 <- post00_2^2
        cell22 <- post01_2^2
        cell33 <- post10_2^2
        cell44 <- post11_2^2
        cell12 <- post00_2 * post01_2
        cell13 <- post00_2 * post10_2
        cell14 <- post00_2 * post11_2
        cell23 <- post01_2 * post10_2
        cell24 <- post01_2 * post11_2
        cell34 <- post10_2 * post11_2

        fourmat <- matrix(NA, nrow = 4, ncol = 4)

        fourmat[1, 1] <- mean(cell11)
        fourmat[2, 2] <- mean(cell22)
        fourmat[3, 3] <- mean(cell33)
        fourmat[4, 4] <- mean(cell44)
        fourmat[1, 2] <- fourmat[2, 1] <- mean(cell12)
        fourmat[1, 3] <- fourmat[3, 1] <- mean(cell13)
        fourmat[1, 4] <- fourmat[4, 1] <- mean(cell14)
        fourmat[2, 3] <- fourmat[3, 2] <- mean(cell23)
        fourmat[2, 4] <- fourmat[4, 2] <- mean(cell24)
        fourmat[3, 4] <- fourmat[4, 3] <- mean(cell34)

        relresults[k, 4, c] <- suppressWarnings(polycor::polychor(fourmat))
      } # end attribute loop
    } # end transition comparison loop

    if (length(attribute.names) == num.atts) {
      rnames <- attribute.names
    } else {
      rnames <- c(paste("Attribute", 1:num.atts, sep = " "))
    }
    cnames <- c(
      "pt bis", "pa forms", "info gain", "polychor", "ave max tr",
      "P(t>.6)", "P(t>.7)", "P(t>.8)", "P(t>.9)", "wt pt bis", "wt pa forms", "wt info gain"
    )

    tp_rnames <- c(1:sum(N))
    tp_cnames <- c("00", "01", "10", "11")

    transclass <- noquote(transclass)
    transclass2 <- matrix(NA, nrow = sum(N), ncol = A * (num.time.points - 1))
    for (m in 1:(A * (num.time.points - 1))) {
      transclass2[, m] <- transclass[, 1, m]
    }
    transclass2 <- noquote(transclass2)

    tc_rnames <- c(1:sum(N))
  } # end transition option = 3 if statement

  ####################################################################################
  ####################################################################################

  dimnames(relresults) <- list(rnames, cnames, mnames)

  dimnames(transposts) <- list(tp_rnames, tp_cnames, tp_mnames)
  dimnames(transclass2) <- list(tc_rnames, tc_cnames)

  relresults <- relresults[, -c(2, 11), ]

  newList <- list(
    "metrics" = round(relresults, 3), "transposts" = transposts,
    "mostlikelytransitions" = transclass2
  )

  return(newList)
}
