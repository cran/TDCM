#' Several data sets for the 'tdcm' package.
#' @name data.tdcm
#' @rdname data.tdcm
#' @order 0
#' @examples
#' \donttest{
#' ## Example 1: T = 2, A = 4
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#' model <- TDCM::tdcm(data, q.matrix, num.time.points = 2)
#' }
NULL

#' @rdname data.tdcm
#' @order 1
#' @docType data
#' @keywords data
#' @format `data.tdcm01` is simulated sample data that has two time points, four
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
#' @format `data.tdcm02` is simulated data that has three time points, two
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
#' @format `data.tdcm03` is simulated data that has three time points, two
#' attributes, one group of size 1500, and three different ten-item Q-matrices
#' for each time point. Anchor items are specified as items 1/1/21 and items
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
#' @format `data.tdcm04` is simulated data that has two time points, four
#' attributes, twenty items, two group of size 800 and 900, respectively, and a
#' single Q-matrix. The format is a list of three:
#'
#' - `data`: a data frame of binary item responses
#' - `q.matrix`: a data frame specifying the Q-matrix
#' - `groups`: a vector specifying the examinee group memberships
"data.tdcm04"

#' @rdname data.tdcm
#' @order 5
#' @docType data
#' @keywords data
#' @format `data.tdcm05` is simulated data that has two has one time point, four
#' attributes, and twenty items. For use with the 1-PLCDM. The format is a list
#' of two:
#' - `data`: a data frame of binary item responses
#' - `q.matrix`: a data frame specifying the Q-matrix
"data.tdcm05"
