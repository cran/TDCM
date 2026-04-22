#' Plotting TDCM Results
#'
#' `tdcm.plot()` visualizes the results from TDCM analyses.
#'
#' @param results results from \code{\link{tdcm.summary}} or \code{\link{mg.tdcm.summary}}
#'
#' @param attribute.names an optional vector of attribute names to include in plots.
#'
#' @param group.names an optional vector of group names to include in plots.
#'
#' @param type an option to specify the type of plot in single group cases; "both" is default and
#' will produce a line plot and a bar chart; "line" will produce a line plot;
#' and "bar" will produce a bar chart.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' \donttest{
#' ## Example 1: T = 2, A = 4
#' data(data.tdcm01, package = "TDCM")
#' dat1 = data.tdcm01$data
#' qmat1 = data.tdcm01$q.matrix
#'
#' #estimate TDCM with invariance assumed and full LCDM
#' m1 = TDCM::tdcm(dat1, qmat1, num.time.points = 2, invariance = TRUE, rule = "LCDM")
#'
#' #summarize results with tdcm.summary function
#' results1 = TDCM::tdcm.summary(m1)
#'
#' #plot results
#' TDCM::tdcm.plot(results1, attribute.names = c("Addition", "Subtraction",
#' "Multiplication", "Division"))
#' }
#'
#' @export
tdcm.plot <- function(results, attribute.names = c(), group.names = c(), type = "both") {

  if (results$numgroups == 1) {#single group TDCM plots
    #pull number of attributes and time points
    num.atts <- nrow(results$growth)
    numtime <- ncol(results$growth)

    #pull growth estimates
    growth <- results$growth

    #pull att names if specified
    if (length(attribute.names) == num.atts) {
      row.names(growth) <- attribute.names}

    if(type == "both"){

      #line plot of growth
      plot(1:numtime, growth[1, ], type = "b", lwd = 3, pch = 16, xlab = "", ylab = "",
           las = 1, xaxt = "n", yaxt = "n", ylim = c(max(round(min(growth) - .10, 1), 0), min(round(max(growth) + .10, 1), 1)), col = 2)
      if (num.atts >= 2) {
        for (i in 2:num.atts){
          graphics::lines(1:numtime, growth[i, ], pch = 16, type = "b", lwd = 3, col = i + 1)
        }
      }

      if (length(attribute.names) == num.atts) {
        graphics::legend("topleft", attribute.names, col = 1 + 1:num.atts, pch = 16, box.lwd = 0, box.col = 1, bty = "n")
      } else {
        graphics::legend("topleft", c(paste("Attribute", 1:num.atts)), col = 1 + 1:num.atts, pch = 16, box.lwd = 0, box.col = 1, bty = "n", cex = 1.1)
      }

      graphics::mtext("Time Point", side = 1, line = 3, cex = 1.3, las = 1)
      graphics::mtext("Proficiency Proportion", side = 2, las = 3, line = 3, cex = 1.3)
      graphics::axis(1, at = c(1:numtime))
      graphics::axis(2, at = c(seq(max(round(min(growth) - .10, 1), 0), min(round(max(growth) + .10, 1), 1), by = .10)), las = 1, cex = 1.1)

      #bar chart of growth
      graphics::barplot(t(growth), xlab = "", ylab = "", las = 1, beside = TRUE, col = 1 + 1:numtime, ylim = c(0, 1))
      graphics::legend("topleft", c(paste("Time", 1:numtime)), col = 1 + 1:numtime, pch = 15, box.lwd = 1, box.col = 1, bty = "n", cex = 1.1)

      graphics::mtext("Proficiency Proportion", side = 2, las = 3, line = 3, cex = 1.3)
      graphics::mtext("Attribute", side = 1, line = 3, cex = 1.3, las = 1)

      for (i in 1:num.atts){
        graphics::text(x = (numtime * (i - 1) + (i - 1) * 1) + 1:numtime + .5, y = growth[i, ] / 2, labels = round(t(growth[i, ]), 2), cex = 1.1)
      }

      print("**Check the plots window for line and bar plots of growth proportions.", quote = FALSE)

    }

    if(type == "line"){

      #line plot of growth
      plot(1:numtime, growth[1, ], type = "b", lwd = 3, pch = 16, xlab = "", ylab = "",
           las = 1, xaxt = "n", yaxt = "n", ylim = c(max(round(min(growth) - .10, 1), 0), min(round(max(growth) + .10, 1), 1)), col = 2)
      if (num.atts >= 2) {
        for (i in 2:num.atts){
          graphics::lines(1:numtime, growth[i, ], pch = 16, type = "b", lwd = 3, col = i + 1)
        }
      }

      if (length(attribute.names) == num.atts) {
        graphics::legend("topleft", attribute.names, col = 1 + 1:num.atts, pch = 16, box.lwd = 0, box.col = 1, bty = "n")
      } else {
        graphics::legend("topleft", c(paste("Attribute", 1:num.atts)), col = 1 + 1:num.atts, pch = 16, box.lwd = 0, box.col = 1, bty = "n", cex = 1.1)
      }

      graphics::mtext("Time Point", side = 1, line = 3, cex = 1.3, las = 1)
      graphics::mtext("Proficiency Proportion", side = 2, las = 3, line = 3, cex = 1.3)
      graphics::axis(1, at = c(1:numtime))
      graphics::axis(2, at = c(seq(max(round(min(growth) - .10, 1), 0), min(round(max(growth) + .10, 1), 1), by = .10)), las = 1, cex = 1.1)

      print("**Check the plots window for a line plot of growth proportions.", quote = FALSE)

    }

    if(type == "bar"){

      #bar chart of growth
      graphics::barplot(t(growth), xlab = "", ylab = "", las = 1, beside = TRUE, col = 1 + 1:numtime, ylim = c(0, 1))
      graphics::legend("topleft", c(paste("Time", 1:numtime)), col = 1 + 1:numtime, pch = 15, box.lwd = 1, box.col = 1, bty = "n", cex = 1.1)

      graphics::mtext("Proficiency Proportion", side = 2, las = 3, line = 3, cex = 1.3)
      graphics::mtext("Attribute", side = 1, line = 3, cex = 1.3, las = 1)

      for (i in 1:num.atts){
        graphics::text(x = (numtime * (i - 1) + (i - 1) * 1) + 1:numtime + .5, y = growth[i, ] / 2, labels = round(t(growth[i, ]), 2), cex = 1.1)
      }

      print("**Check the plots window for a bar chart of growth proportions.", quote = FALSE)

    }

  }#end single group TDCM plot function


  else {
    tdcm.plot.mg(results = results, attribute.names = attribute.names, group.names = group.names)
  }#end multiple group TDCM plot function


}
