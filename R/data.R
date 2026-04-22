#' Several data sets for the \pkg{TDCM} package.
#'
#' @name data.tdcm
#' @rdname data.tdcm
#' @order 0
#' @examples
#' \donttest{
#'
#' #############################
#' ## Example 1: T = 2, A = 4 ##
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#' model <- TDCM::tdcm(data, q.matrix, num.time.points = 2)
#' results <- TDCM::tdcm.summary(model)
#' results$item.parameters
#' results$growth.effects
#'
#' #############################
#' ## Example 3: T = 3, A = 2 ##
#' data <- data.tdcm03$data
#' q1 <- data.tdcm03$q.matrix.1
#' q2 <- data.tdcm03$q.matrix.2
#' q3 <- data.tdcm03$q.matrix.3
#' q <- data.tdcm03$q.matrix.stacked
#'
#' #TDCM with anchor items constrained
#' m1 <- tdcm(data, q, num.time.points = 3, num.q.matrix = 3,
#' anchor = c(1,11,1,21,14,24), num.items = c(10,10,10))
#'
#' #TDCM without anchor items
#' m2 <- tdcm(data, q, num.time.points = 3, num.q.matrix = 3, num.items = c(10,10,10))
#'
#' #Compare models to assess measurement invariance
#' tdcm.compare(m1, m2)
#'
#' #Summarize model 1
#' r1 = tdcm.summary(m1, transition.option = 3)
#' r1$item.parameters
#' r1$growth
#' r1$growth.effects
#'
#' }
NULL

#' @rdname data.tdcm
#' @order 1
#' @docType data
#' @keywords data
#' @format `data.tdcm01` is a simulated dataset with two time points, four
#' attributes, twenty items, one group of size 1000, and a single Q-matrix. The
#' format is a list of two:
#'
#' - `data`: a data frame of binary item responses
#' - `q.matrix`: a data frame specifying the Q-matrix
"data.tdcm01"

#' @rdname data.tdcm
#' @order 2
#' @docType data
#' @keywords data
#' @format `data.tdcm02` is a simulated dataset with three time points, two
#' attributes, ten items, one group of size 2500, and a single Q-matrix. The
#' format is a list of two:
#'
#' - `data`: a data frame of binary item responses
#' - `q.matrix`: a data frame specifying the Q-matrix
"data.tdcm02"

#' @rdname data.tdcm
#' @order 3
#' @docType data
#' @keywords data
#' @format `data.tdcm03` is a simulated dataset with three time points, two
#' attributes, one group of size 1500, and three different ten-item Q-matrices
#' for each time point. Anchor items are specified as items 1/11/21 and items
#' 14/24. The format is a list of five:
#'
#' - `data`: a data frame of binary item responses
#' - `q.matrix.1`: a data frame specifying the Q-matrix for the first time point
#' - `q.matrix.2`: a data frame specifying the Q-matrix for the second time point
#' - `q.matrix.3`: a data frame specifying the Q-matrix for the third time point
#' - `q.matrix.stacked`: data frame specifying the combined Q-matrix for all time
#'   points
"data.tdcm03"

#' @rdname data.tdcm
#' @order 4
#' @docType data
#' @keywords data
#' @format `data.tdcm04` is a simulated dataset with two time points, four
#' attributes, twenty items, two group of size 800 and 900, respectively, and a
#' single Q-matrix. The format is a list of three:
#'
#' - `data`: a data frame of binary item responses
#' - `q.matrix`: a data frame specifying the Q-matrix
#' - `groups`: a vector specifying the examinee group membership
"data.tdcm04"

#' @rdname data.tdcm
#' @order 5
#' @docType data
#' @keywords data
#' @format `data.tdcm05` is a simulated dataset with one time point, four
#' attributes, and twenty items. For use with the 1-PLCDM. The format is a list
#' of two:
#' - `data`: a data frame of binary item responses
#' - `q.matrix`: a data frame specifying the Q-matrix
"data.tdcm05"
