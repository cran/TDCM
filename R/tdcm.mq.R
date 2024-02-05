#' @title Utility function to estimate TDCM with multiple Q-matrices.
#'
#' @param data item response data
#' @param q.matrix specified qumatrices
#' @param num.time.points number of time points
#' @param invariance invariance assumption (T or F)
#' @param rule specific DCM to estimate
#' @param num.q.matrix number of Q-matrices
#' @param num.items number of items for each Q-matrix
#' @param anchor anchor items specified as pairs in a vector
#' @keywords internal
#' @noRd
tdcm.mq <- function(data, q.matrix, num.time.points, invariance = TRUE, rule = "GDINA",
                    num.q.matrix = 1, num.items = c(), anchor = c()) {

  # Initial Data Sorting
  n.items <- ncol(data) # Total Items
  items <- num.items
  N <- nrow(data) # Number of Examinees
  n.att <- ncol(q.matrix) # Number of Attributes

  qnew <- matrix(0, ncol = n.att * num.time.points, nrow = n.items)
  qnew[(1:items[1]), (1:n.att)] <- as.matrix(q.matrix[(1:items[1]), (1:n.att)])
  for (z in 2:num.time.points) { # stack the Q-matrices
    qnew[(sum(items[2:z]) + 1):(sum(items[1:z])), ((1:n.att) + (z - 1) * 2)] <- as.matrix(q.matrix[(sum(items[2:z]) + 1):(sum(items[1:z])), (1:n.att)])
  } # for

  if (length(anchor) == 0) {
    # if no anchor items are specified
    tdcm <- CDM::gdina(data, qnew, linkfct = "logit", method = "ML", rule = rule, progress = FALSE)
  } else {
    # if anchor items are specified, make them equal in the design matrix
    tdcm.1 <- tdcm.base(data, qnew, rule)
    c0 <- tdcm.1$coef
    c.0 <- nrow(c0)
    delta.designmatrix <- diag(nrow = c.0, ncol = c.0)

    for (a in 1:(length(anchor) / 2)) { # collect indices of anchored item parms
      row.anchor <- as.numeric(rownames(c0[c0$itemno %in% anchor[(2 * a - 1):(2 * a)], ]))
      nrows <- length(row.anchor) / 2
      for (b in 1:nrows) {
        delta.designmatrix[row.anchor[nrows + b], ] <- delta.designmatrix[row.anchor[b], ]
      } # end loop to build design matrix
    } # for

    # remove constrained item columns
    rem.cols <- as.numeric(rownames(c0[c0$itemno %in% anchor[c(FALSE, TRUE)], ]))
    delta.designmatrix <- delta.designmatrix[, -rem.cols]
    tdcm <- CDM::gdina(
      data,
      qnew,
      linkfct = "logit",
      method = "ML",
      progress = FALSE,
      delta.designmatrix = delta.designmatrix,
      rule = rule
    ) # tdcm
  } # end else for when anchor items are specified

  tdcm$invariance <- FALSE

  return(tdcm)

} # if
