#' @title Estimating item influence measures.
#'
#' @description Function to estimate estimate item influence measures. Code adapted from (Jurich & Madison, 2023).
#' This function is not available for longitudinal DCMs.
#'
#' @param model a previously calibrated model; an object of class \code{gdina}.
#'
#' @param data a required \eqn{N \times I} matrix. Binary item responses are in the columns.
#'
#' @param fullcorrelation optional logical argument indicating a full or reduced response-classification correlation matrix.
#'
#' @param progress An optional logical indicating whether the function should print the progress of estimation.
#'
#' @details For DCMs, item influence quantifies how much an item impacts classifications. Given an estimated DCM and item response data, this function estimates five item influence measures, including
#' item pull, item override, proportion of attribute information, response-classification correlation (corr1), and response-posterior correlation (corr2).
#' @note
#' Currently, this function currently only runs on DCMs estimated at a single time point. It will not run properly for TDCM objects.
#'
#' @return A list containing several item influence measures.
#'
#' @examples
#' \donttest{
#' ## Item influence illustration
#' #load data (simulated based on Jurich and Bradshaw (2014))
#' qmatrix <- CDM::data.sda6$q.matrix
#' responses <- CDM::data.sda6$data
#'
#' #Estimate the full LCDM
#' model1 <- CDM::gdina(responses, qmatrix, linkfct = "logit", method = "ML")
#'
#' #Estimate item influence measures
#' influence <- TDCM::item.influence(model1, responses)
#'
#' #Summarize influence statistics
#' influence$Pull #item pull
#' influence$Override #item override
#' influence$Information #proportion of attribute information
#' influence$Correlation1 #correlation of responses and classifications
#' influence$Correlation2 #correlation of responses and posterior probabilities
#'
#' }
#'
#' @references
#'
#' Jurich, D. & Madison, M. J. (2023). Measuring item influence for diagnostic classification models. \emph{Educational Assessment}.
#'
#' @export
item.influence <- function(model, data, fullcorrelation = FALSE, progress = TRUE) {

  if (progress == TRUE) {
    print("Calclating item influence measures. Progress = 0%.", quote = FALSE)
  }

  #attribute profile and mastery proportions
  profileprops <- model$attribute.patt
  masteryprops <- model$skill.patt

  #Individual attribute profiles
  postprobs <- as.matrix(model$pattern[, 6:ncol(model$pattern)])
  originalclass <- (postprobs > .5) * 1

  #### now let's compute the influence measures
  I <- length(unique(model$coef$item))
  J <- length(model$NAttr)
  Q <- model$q.matrix

  #1) item pull: pull0 = pull towards non-mastery, pull1 = pull towards mastery
  pull0 <- matrix(NA, nrow = I, ncol = J) #conditional prob of non-mastery, given incorrect response
  pull1 <- matrix(NA, nrow = I, ncol = J) #conditional prob of mastery, given correct response


  if (ncol(Q) == 1) {
    for (i in 1:I) {
      pull0[i, ] <- (1 - mean(originalclass[which(data[, i] == 0), ])) * Q[i, ]
      pull1[i, ] <- mean(originalclass[which(data[, i] == 1), ]) * Q[i, ]
    }
  }
  else {
    for (i in 1:I) {
      pull0[i, ] <- (1 - colMeans(originalclass[which(data[, i] == 0), ])) * Q[i, ]
      pull1[i, ] <- colMeans(originalclass[which(data[, i] == 1), ]) * Q[i, ]
    }
  }


  #Remove missing q-matrix entires, just in case one of the pulls equals exactly zero
  pull0[!sapply(as.data.frame(Q), as.logical)] <- NA
  pull1[!sapply(as.data.frame(Q), as.logical)] <- NA

  #Add Column Names
  colnames(pull0) <- sprintf("Pull0_Att%d", seq(1:J))
  colnames(pull1) <- sprintf("Pull1_Att%d", seq(1:J))

  #Simple Structure Only
  if (all(rowSums(Q) == 1)) {
    pull0 <- as.numeric(round(t(pull0), 2)[!is.na(t(pull0))])
    pull1 <- as.numeric(round(t(pull1), 2)[!is.na(t(pull1))])
    Pull <- data.frame(Item = 1:I, Attribute = max.col(Q), pull0, pull1)
  }

  #Complex Structure
  if (any(rowSums(Q) != 1)) {
    pull0 <- round(pull0, 2)
    pull1 <- round(pull1, 2)
    Pull <- data.frame(Item = 1:I, pull0, pull1)
  }

  Results <- list()

  Results[[1]] <- Pull

  #2) item override: estimate model omitting item, compare classifications, how many switch?
  override <- matrix(NA, nrow = I, ncol = J)
  for (i in 1:I) {

    #remove item i from responses and q-matrix
    q1 <- as.matrix(Q[-i, ])
    responses1 <- as.matrix(data[, -i])

    #estimate lcdm with item removed
    m2 <- CDM::gdina(responses1, q1, linkfct = "logit", progress = FALSE)

    #get the new classificationss
    newclass <- as.matrix(m2$pattern[, 6:ncol(model$pattern)] > .5) * 1

    #now compare to original classifications, how many switch
    override[i, ] <- (1 - colMeans(originalclass == newclass)) * Q[i, ]

    if (i %% 4 == 0 & progress == TRUE) {
      print(paste("Calclating item influence measures. Progress = ", round(i / I * 100, 0), "%.", sep = ""), quote = FALSE)
    }
    if (i == I & progress == TRUE) {
      print("Routine finished. Check results.", quote = FALSE)
    }
  }

  #Some cleanup
  override[!sapply(as.data.frame(Q), as.logical)] <- NA
  override <- round(override, 2)
  colnames(override) <- sprintf("Att%d", seq(1:J))

  #Simple Structure Only
  if (all(rowSums(Q) == 1)) {
    override <- as.numeric(t(override)[!is.na(t(override))])
    override <- data.frame(Item = 1:I, Attribute = max.col(Q), override)
  }

  #Complex Structure
  if (any(rowSums(Q) != 1)) {
    override <- data.frame(Item = 1:I, override)
  }

  Results[[2]] <- override

  #3) item proportion of attribute information
  propinfo <- matrix(NA, nrow = I, ncol = J)

  #attribute information
  info <- as.matrix(CDM::cdi.kli(model)$summary[, 3:ncol(CDM::cdi.kli(model)$summary)])
  #row 1 is total attribute information
  #rows 2 - 6 are each individual item's attribute information

  #Complex Item Loop
  for (i in 1:I) {
    for (j in 1:J) {

      #item proportion of attribute information
      propinfo[i, j] <- info[i + 1, j] / info[1, j]
    }
  }

  #Clean up
  propinfo[!sapply(as.data.frame(Q), as.logical)] <- NA
  propinfo <- round(propinfo, 2)
  colnames(propinfo) <- sprintf("Att%d", seq(1:J))

  #Simple Structure Only
  if (all(rowSums(Q) == 1)) {
    propinfo <- as.numeric(t(propinfo)[!is.na(t(propinfo))])
    propinfo <- data.frame(Item = 1:I, Attribute = max.col(Q), propinfo)
  }

  #Complex Structure
  if (any(rowSums(Q) != 1)) {
    propinfo <- data.frame(Item = 1:I, propinfo)
  }

  Results[[3]] <- propinfo

  #Correlation Descriptive of Responses with class
  Correl1 <- round(stats::cor(data, originalclass), 2)
  rownames(Correl1) <- NULL
  colnames(Correl1) <- sprintf("Att%d", seq(1:J))

  #Correlation Descriptive of Responses with posterior
  Correl2 <- round(stats::cor(data, postprobs), 2)
  rownames(Correl2) <- NULL
  colnames(Correl2) <- sprintf("Att%d", seq(1:J))

  #Only looking at items associated with attribute
  if (fullcorrelation == FALSE) {
    Correl1[!Q] <- NA
    Correl2[!Q] <- NA
  }

  #Simple Structure Only
  if (all(rowSums(Q) == 1 & fullcorrelation == FALSE)) {
    Correl1 <- as.numeric(t(Correl1)[!is.na(t(Correl1))])
    Correl1 <- data.frame(Item = 1:I, Attribute = max.col(Q), Correl1)

    Correl2 <- as.numeric(t(Correl2)[!is.na(t(Correl2))])
    Correl2 <- data.frame(Item = 1:I, Attribute = max.col(Q), Correl2)
  }

  #Complex Structure
  if (any(rowSums(Q) != 1 | fullcorrelation == TRUE)) {
    Correl1 <- data.frame(Item = 1:I, Correl1)
    Correl2 <- data.frame(Item = 1:I, Correl2)
  }


  Results[[4]] <- Correl1
  Results[[5]] <- Correl2

  names(Results) <- c("Pull", "Override", "Information", "Correlation1", "Correlation2")

  return(Results)

}
