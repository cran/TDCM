#' Utility function in \pkg{TDCM}
#'
#' @param model model from \bold{tdcm} estimation
#' @param num.atts number of attributes
#' @param num.items number of items at each time point
#' @param num.time.points number of time points
#' @keywords internal
#' @noRd
summary.param <- function(model, num.atts, num.items, num.time.points) {

  temp.param.names1 <- NULL
  param.names1 <- NULL
  lscoefs0 <- NULL

  invariance <- model$invariance

  param.names <- c() ### Empty vector for names
  for (k in 1:num.atts) {
    temp.p.names <- as.data.frame(gtools::combinations(n = num.atts, r = k)) # Creates all possible combinations using k numbers
    for (i in 1:nrow(temp.p.names)) {
      param.names <- c(param.names, paste0("Attr", toString(temp.p.names[i, ])))
    }
  }
  p.names <- gsub(", ", "-Attr", param.names) # Fixes name scheme to match model$coef
  p.names <- append("", p.names)


  ind.vals <- seq(from = 0, to = length(p.names))
  for (i in 1:length(p.names)) {
    p.name.ind <- which(model$coef[, 8] == p.names[i]) # Creates index of all rows that match attribute combo
    p.coefs <- as.matrix(model$coef[p.name.ind, c(2, 8, 6)]) # Stores the columns wanted from model$coefs that match indexed rows
    p.coefs <- noquote(p.coefs)
    assign(paste0("lscoefs", ind.vals[i]), p.coefs) # Creates the names for each sub-dataframe
  }
  lscoef.og <- model$coef # Creates a default with no adjustments
  ls.p.names <- ls(pattern = "lscoefs") # creates a list of all the lscoefs subset NAMES we made
  ls.p.names <- gtools::mixedsort(ls.p.names)
  ls.p.names <- append("lscoef.og", ls.p.names)
  param <- mget(ls.p.names[1:length(p.names)])

  order.check <- c() # This block finds the highest order
  for (i in 1:num.items) {
    order.check <- c(order.check, sum(model$q.matrix[i, ])) # Sums each row of q-matrix
  }
  h.order <- max(order.check) # Returns highest order term

  for (i in 1:num.time.points) {
    assign(paste0("t.param.names", i), c())
  }

  for (t in 1:num.time.points) {
    assign(paste0("param.names", t), c())
  }

  for (k in 1:h.order) { # Set of loops finds combinations for the highest order
    for (t in 1:num.time.points) {
      assign(paste0("temp.param.names", t), as.data.frame(gtools::combinations(num.atts, r = k) + ((t - 1) * num.atts)))
    }
    for (i in 1:nrow(temp.param.names1)) {
      for (t in 1:num.time.points) {
        assign(paste0("param.names", t), c(get(paste0("param.names", t)), paste0(toString(get(paste0("temp.param.names", t))[i, ]))))
      }
    }
  }

  s.parm <- matrix(NA, nrow = num.items, ncol = length(param.names1) + 1) # Creates appropriate empty matrix
  num.count <- nchar(gsub("\\D", "", param.names1)) # Counts the number of integers
  parm.c.names <- c(paste0("\U03BB", 0)) # Sets first column name to lambda0
  for (i in 1:(length(param.names1))) { # This loop goes from 2 to number of parameter names
    for (a in 1:num.atts) {
      if (num.count[i] == a) {
        temp.name <- paste0("\U03BB", a, ",", gsub(", ", "", param.names1[i])) # Creates name lambda1,1 or lambda1,2, or lambda3,123
        parm.c.names <- c(parm.c.names, temp.name)
      } else {}
    }
  }
  colnames(s.parm) <- parm.c.names # Fill in the column names
  rownames(s.parm) <- lscoefs0[, 1] # Fills in the row names using the intercepts which each item should have
  s.parm[, 1] <- lscoefs0[, 3] # Fills in first column with intercepts

  for (t in 1:num.time.points) {
    assign(paste0("new.p.names", t), c())
  }
  for (i in 1:length(param.names1)) {
    for (t in 1:num.time.points) {
      assign(paste0("new.p.names", t), c(get(paste0("new.p.names", t)), paste0("Attr", gsub(", ", "-Attr", get(paste0("param.names", t))[i]))))
      assign(paste0("new.p.names", t), noquote(get(paste0("new.p.names", t))))
    }
  }


  mod <- model$coef

  for (t in 1:num.time.points) {
    assign(paste0("partial.index", t), c())
  }

  for (i in 2:ncol(s.parm)) {
    index <- c()
    for (t in 1:num.time.points) {
      assign(paste0("partial.index", t), as.matrix(mod[which(mod[, 8] == get(paste0("new.p.names", t))[i - 1]), ]))
      index <- rbind(index, get(paste0("partial.index", t)))
    }
    if (nrow(index) != 0) {
      s.parm[index[, 2], i] <- index[, 6]
    }
  }


  class(s.parm) <- "numeric" # Coerces everything to numeric
  s.parm <- round(s.parm, 3) # Rounds
  row.names(s.parm) <- paste("Item ", 1:num.items, sep = "")

  param <- s.parm # Renames for final list return
  param[is.na(param)] <- "  -- " # replace NA
  param <- noquote(param)
  if (invariance == TRUE) {
    param <- param[1:(nrow(param) / num.time.points), ] #
  } else {}

  return(param)
}
