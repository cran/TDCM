#' Utility function in \pkg{TDCM}
#'
#' @param model gdina object from tdcm function
#' @param num.atts number of attributes
#' @param num.time.points number of time points
#' @param attribute.names optional argument to specify attribute names
#' @keywords internal
#' @noRd
summary.option1 <- function(model, num.atts, num.time.points, attribute.names) {
  A <- num.atts
  transition.option <- 1
  growth <- matrix(NA, num.atts, num.time.points)

  cnames.growth <- c()
  for (t in 1:num.time.points) {
    temp.growth.c.names <- c(paste0("T", t, "[1]"))
    cnames.growth <- append(cnames.growth, temp.growth.c.names)
  }

  rnames.growth <- c()
  for (i in 1:num.atts) {
    temp.g.rname <- paste("Attribute", i, sep = " ")
    rnames.growth <- c(rnames.growth, temp.g.rname)
  }
  matrix.names.growth <- c()

  trans <- array(NA, c(2, 2, num.atts)) # Creates an array of empty 2x2 matrices, one for each attribute
  trans.cnames <- c(paste("T", num.time.points, " [0]", sep = ""), paste("T", num.time.points, " [1]", sep = ""))
  trans.rnames <- c("T1 [0]", "T1 [1]")

  matrix.names.trans <- c()

  for (t in 2:num.time.points) {
    for (j in 1:num.atts) {
      temp.ind00 <- which(model$attribute.patt.splitted[, j] == 0 & model$attribute.patt.splitted[, j + ((t - 1) * num.atts)] == 0) # First step is to create an index
      temp.ind01 <- which(model$attribute.patt.splitted[, j] == 0 & model$attribute.patt.splitted[, j + ((t - 1) * num.atts)] == 1) # using $attibute.patt.splitted
      temp.ind10 <- which(model$attribute.patt.splitted[, j] == 1 & model$attribute.patt.splitted[, j + ((t - 1) * num.atts)] == 0)
      temp.ind11 <- which(model$attribute.patt.splitted[, j] == 1 & model$attribute.patt.splitted[, j + ((t - 1) * num.atts)] == 1)

      temp.sum00 <- sum(model$attribute.patt[temp.ind00, 1]) # Finds the actual sums of each pattern
      temp.sum01 <- sum(model$attribute.patt[temp.ind01, 1])
      temp.sum10 <- sum(model$attribute.patt[temp.ind10, 1])
      temp.sum11 <- sum(model$attribute.patt[temp.ind11, 1])

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

    temp.sum00 <- sum(model$attribute.patt[temp.ind00, 1])
    temp.sum01 <- sum(model$attribute.patt[temp.ind01, 1])
    temp.sum10 <- sum(model$attribute.patt[temp.ind10, 1])
    temp.sum11 <- sum(model$attribute.patt[temp.ind11, 1])

    base <- c(temp.sum00, temp.sum01, temp.sum10, temp.sum11)
    r1 <- sum(temp.sum00, temp.sum01)
    r2 <- sum(temp.sum10, temp.sum11)
    props <- base / c(r1, r1, r2, r2)

    temp.trans <- matrix(props, nrow = 2, ncol = 2, byrow = TRUE) # Creates a temporary 2x2 to slot into the array
    trans[, , j] <- round(temp.trans, 3) # Replaces the empty matrix in array slot J with the temporary matrix

    if (length(attribute.names) == A) {
      temp.name <- paste(attribute.names[j], ": Time 1 to Time ", num.time.points, sep = "")
    } else {
      temp.name <- paste("Attribute ", j, ": Time 1 to Time ", num.time.points, sep = "")
    } # Creates temporary name for each matrix
    matrix.names.trans <- c(matrix.names.trans, temp.name) # Combines matrix names into list
  }

  #Growth effect sizes
  growth.effects <- matrix(NA, num.atts, 5)
  growth.effects[, 1] <- growth[, 1] #T1 proportion
  growth.effects[, 2] <- growth[, num.time.points] #Tk proportion

  #growth = difference
  diff <- growth[, num.time.points] - growth[, 1]
  growth.effects[, 3] <-  diff

  #odds ratio
  or1 <- growth[, num.time.points] / (1 - growth[, num.time.points])
  or2 <- growth[, 1] / (1 -growth[, 1])
  growth.effects[, 4] <- round(or1 / or2, 2)

  #Cohen's h: arcsin difference in proportions
  ar1 <- 2 * asin(sqrt(growth[, num.time.points]))
  ar2 <- 2 * asin(sqrt(growth[, 1]))
  growth.effects[, 5] <- round(ar1 - ar2, 2)
  colnames(growth.effects) <- c("T1[1]", paste("T", num.time.points,"[1]", sep=""),
                                 "Growth", "Odds Ratio", "Cohen`s h")

  if (length(attribute.names) == A) {
    dimnames(growth) <- list(attribute.names, cnames.growth)
    rownames(growth.effects) <- attribute.names
  } else {
    dimnames(growth) <- list(rnames.growth, cnames.growth)
    rownames(growth.effects) <- rnames.growth

  }

  dimnames(trans) <- list(trans.rnames, trans.cnames, matrix.names.trans)

  # compute transition reliability
  if (model$progress == TRUE) {
    print("Summarizing results...", quote = FALSE)
  }
  if (length(attribute.names) == A) {
    rel <- tdcm.rel(model, num.atts, num.time.points,
                    transition.option = transition.option, attribute.names = attribute.names
    )
  } else {
    rel <- tdcm.rel(model, num.atts, num.time.points, transition.option = transition.option)
  }





  newlist1 <- list("trans" = trans, "growth" = growth,
                   "growth.effects" = growth.effects,  "rel" = rel)
  return(newlist1)
}
