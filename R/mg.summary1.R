#' Utility function in \pkg{TDCM}
#'
#' @param model gdina object from mg.tdcm function
#' @param num.atts number of attributes
#' @param num.groups number of groups
#' @param num.time.points number of time points
#' @param attribute.names optional vector to specify attribute names
#' @param group.names optional vector to specify group names
#' @keywords internal
#' @noRd
mg.summary1 <- function(model, num.atts, num.groups, num.time.points, attribute.names, group.names) {
  transition.option <- 1
  A <- num.atts

  # create matrices to store growth and trans probality matrices for each group
  all.growth <- array(NA, dim = c(num.atts, num.time.points, num.groups))
  all.trans <- array(NA, dim = c(2, 2, num.atts, num.groups))

  for (g in 1:num.groups) {
    growth <- matrix(NA, num.atts, num.time.points)

    cnames.growth <- c()
    for (t in 1:num.time.points) {
      temp.growth.c.names <- c(paste0("T", t, "[1]"))
      cnames.growth <- append(cnames.growth, temp.growth.c.names)
    }

    rnames.growth <- c()
    for (i in 1:num.atts) {
      if (length(attribute.names) == num.atts) {
        temp.g.rname <- attribute.names[i]
      } else {
        temp.g.rname <- paste("Attribute", i, sep = " ")
      } # Creates temporary row name per iteration
      rnames.growth <- c(rnames.growth, temp.g.rname) # Combines row names into list
    }
    matrix.names.growth <- c()

    trans <- array(NA, c(2, 2, num.atts)) # Creates an array of empty 2x2 matrices, one for each attribute
    trans.cnames <- c(paste("T", num.time.points, " [0]", sep = ""), paste("T", num.time.points, " [1]", sep = ""))
    trans.rnames <- c("T1 [0]", "T1 [1]")
    matrix.names <- c()

    for (t in 2:num.time.points) { # open time points loop

      for (j in 1:num.atts) { # open attribute loop


        temp.ind00 <- which(model$attribute.patt.splitted[, j] == 0 & model$attribute.patt.splitted[, j + ((t - 1) * num.atts)] == 0) # First step is to create an index
        temp.ind01 <- which(model$attribute.patt.splitted[, j] == 0 & model$attribute.patt.splitted[, j + ((t - 1) * num.atts)] == 1) # using $attibute.patt.splitted
        temp.ind10 <- which(model$attribute.patt.splitted[, j] == 1 & model$attribute.patt.splitted[, j + ((t - 1) * num.atts)] == 0)
        temp.ind11 <- which(model$attribute.patt.splitted[, j] == 1 & model$attribute.patt.splitted[, j + ((t - 1) * num.atts)] == 1)

        temp.sum00 <- sum(model$attribute.patt[temp.ind00, g]) # Finds the actual sums of each pattern
        temp.sum01 <- sum(model$attribute.patt[temp.ind01, g])
        temp.sum10 <- sum(model$attribute.patt[temp.ind10, g])
        temp.sum11 <- sum(model$attribute.patt[temp.ind11, g])


        temp.growth <- temp.sum01 + temp.sum11 # Provides timepoint.x proficiency proportion
        growth[j, 1] <- round(temp.sum10 + temp.sum11, 3) # Provides baserates
        growth[j, t] <- round(temp.growth, 3) ####
      }
    }

    for (j in 1:num.atts) {
      temp.ind00 <- which(model$attribute.patt.splitted[, j] == 0 & model$attribute.patt.splitted[, j + ((num.time.points - 1) * num.atts)] == 0) # First step is to create an index
      temp.ind01 <- which(model$attribute.patt.splitted[, j] == 0 & model$attribute.patt.splitted[, j + ((num.time.points - 1) * num.atts)] == 1) # using $attibute.patt.splitted
      temp.ind10 <- which(model$attribute.patt.splitted[, j] == 1 & model$attribute.patt.splitted[, j + ((num.time.points - 1) * num.atts)] == 0)
      temp.ind11 <- which(model$attribute.patt.splitted[, j] == 1 & model$attribute.patt.splitted[, j + ((num.time.points - 1) * num.atts)] == 1)

      temp.sum00 <- sum(model$attribute.patt[temp.ind00, g]) # Finds the actual sums of each pattern
      temp.sum01 <- sum(model$attribute.patt[temp.ind01, g])
      temp.sum10 <- sum(model$attribute.patt[temp.ind10, g])
      temp.sum11 <- sum(model$attribute.patt[temp.ind11, g])

      base <- c(temp.sum00, temp.sum01, temp.sum10, temp.sum11)
      r1 <- sum(temp.sum00, temp.sum01)
      r2 <- sum(temp.sum10, temp.sum11)
      props <- base / c(r1, r1, r2, r2)

      temp.trans <- matrix(props, nrow = 2, ncol = 2, byrow = TRUE) # Creates a temporary 2x2 to slot into the array
      trans[, , j] <- round(temp.trans, 3) # Replaces the empty matrix in array slot J with the temporary matrix

      temp.name <- paste("Attribute", j, sep = " ") # Creates temporary name for each matrix
      matrix.names <- c(matrix.names, temp.name) # Combines matrix names into list
    }

    dimnames(all.trans) <- list(trans.rnames, trans.cnames, matrix.names)

    all.growth[, , g] <- growth
    all.trans[, , , g] <- trans
  } # end group loop

  # give growth and transition matrices names
  if (length(attribute.names) == num.atts) {
    rownames(all.growth) <- attribute.names
  } else {
    rownames(all.growth) <- paste("Attribute", 1:num.atts, sep = " ")
  }
  colnames(all.growth) <- cnames.growth

  if (length(group.names) == num.groups) {
    dimnames(all.growth)[[3]] <- group.names
  } else {
    dimnames(all.growth)[[3]] <- paste("Group", 1:num.groups, sep = " ")
  }

  if (length(attribute.names) == num.atts) {
    dimnames(all.trans)[[3]] <- paste(attribute.names, ": Time 1 to Time ", num.time.points, sep = "")
  } else {
    dimnames(all.trans)[[3]] <- paste("Attribute ", 1:num.atts, ": Time 1 to Time ", num.time.points, sep = "")
  }

  if (length(group.names) == num.groups) {
    dimnames(all.trans)[[4]] <- group.names
  } else {
    dimnames(all.trans)[[4]] <- paste("Group", 1:num.groups, sep = " ")
  }

  # compute transition reliability
  complexity <- num.atts * num.time.points * num.groups
  if (complexity < 20) {
    esttime <- round(stats::runif(1, 40, 60), 0)
  } else {
    esttime <- round(stats::runif(1, 30, 40), 0)
  }

  if (model$progress == TRUE) {
    print(paste("Summarizing results, progress = ", esttime, "%...", sep = ""), quote = FALSE)
  }
  if (length(attribute.names) == A) {
    rel <- tdcm.rel(model, num.atts, num.time.points,
                    transition.option = transition.option, attribute.names = attribute.names
    )
  } else {
    rel <- tdcm.rel(model, num.atts, num.time.points, transition.option = transition.option)
  }

  newlist1 <- list("trans" = all.trans, "growth" = all.growth, "rel" = rel)
  return(newlist1)
}
