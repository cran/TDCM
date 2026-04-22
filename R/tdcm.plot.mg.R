#' Utility function in \pkg{TDCM}
#'
#' @param results results from \bold{mg.tdcm.summary}
#' @param attribute.names optional vector of attribute names to include in plots.
#' @param group.names optional vector of group names to include in plots.
#' @keywords internal
#' @noRd
tdcm.plot.mg <- function(results, attribute.names = c(), group.names = c()) {
  # mgTDCM plots

  # pull number of groups, attributesn and time points
  numgroups <- results$numgroups
  num.atts <- nrow(results$growth)
  numtime <- ncol(results$growth)

  # pull growth estimates
  growth <- results$growth

  # pull att names if specified
  if (length(attribute.names) == num.atts) {
    row.names(growth) <- attribute.names
  }

  # pull group names if specified
  if (length(group.names) == numgroups) {
    dimnames(growth)[[3]] <- group.names
  }

  # line plots, one for each attribute
  for (i in 1:num.atts) {
    plot(1:numtime, growth[i, , 1],
         type = "b", lwd = 3, pch = 16, xlab = "", ylab = "",
         las = 1, xaxt = "n", yaxt = "n", ylim = c(max(round(min(growth) - .10, 1), 0), min(round(max(growth) + .10, 1), 1)), col = 2)

    if (length(attribute.names) == 0) {
      graphics::title(paste("Attribute ", i, sep = ""))
    } else {
      graphics::title(paste(attribute.names[i]))
    }

    for (j in 2:numgroups) {
      graphics::lines(1:numtime, growth[i, , j], pch = 16, type = "b", lwd = 3, col = 1 + j)
    } # end group loop

    # add details to line plot
    graphics::mtext("Time Point", side = 1, line = 3, cex = 1.3, las = 1)
    graphics::mtext("Proficiency Proportion", side = 2, las = 3, line = 3, cex = 1.3)

    graphics::axis(1, at = c(1:numtime))
    graphics::axis(2, at = c(seq(max(round(min(growth) - .10, 1), 0), min(round(max(growth) + .10, 1), 1), by = .10)), las = 1, cex = 1.1)

    if (length(group.names) == numgroups) {
      graphics::legend("topleft", group.names, col = 1 + 1:numgroups, pch = 16, box.lwd = 0, box.col = 1, bty = "n")
    } else {
      graphics::legend("topleft", c(paste("Group", 1:numgroups)), col = 1 + 1:num.atts, pch = 16, box.lwd = 0, box.col = 1, bty = "n", cex = 1.1)
    }
  } # end attribute loop

  # bar chart of growth
  for (i in 1:num.atts) {
    graphics::barplot(growth[i, , ],
            xlab = "", ylab = "", las = 1, beside = TRUE, col = 1 + 1:numtime,
            ylim = c(0, round(min(max(growth[i, , ]) + .2, 1), 2)))
    graphics::legend("topleft", c(paste("Time", 1:numtime)), col = 1 + 1:numtime, pch = 15, box.lwd = 1, box.col = 1, bty = "n", cex = 1.1)

    graphics::mtext("Proficiency Proportion", side = 2, las = 3, line = 3, cex = 1.3)
    graphics::mtext("Group", side = 1, line = 3, cex = 1.3, las = 1)

    if (length(attribute.names) == 0) {
      graphics::title(paste("Attribute ", i, sep = ""))
    } else {
      graphics::title(paste(attribute.names[i]))
    }

    for (j in 1:numgroups) {
      graphics::text(x = (numtime * (j - 1) + (j - 1) * 1) + 1:numtime + .5, y = growth[i, , j] / 2, labels = round(growth[i, , j], 2), cex = 1.1)
    }
  } # end attribute loop

  print("**Check the plots window for line and bar plots for group growth proportions.", quote = FALSE)

}
