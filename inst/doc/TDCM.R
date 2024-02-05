## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  dpi = 300,
  out.width = "100%"
) # set output options

## ----eval = TRUE--------------------------------------------------------------
standards <- paste0("4.MD.", 1:4)
standards

## ----eval = TRUE--------------------------------------------------------------
# Load the TDCM package and sample dataset
library(TDCM)
data(data.tdcm01, package = "TDCM")

# Get item responses from sample data.
data <- data.tdcm01$data
head(data)

# Get Q-matrix from sample data and rename the attributes to match the standard.
q.matrix <- data.tdcm01$q.matrix
colnames(q.matrix) <- standards
q.matrix

## ----eval = TRUE--------------------------------------------------------------
# Calibrate TDCM with measurement invariance assumed, full LCDM
model1 <- tdcm(data, q.matrix, num.time.points = 2)

## ----eval = TRUE--------------------------------------------------------------
# Summarize the results
results1 <- tdcm.summary(model1, num.time.points = 2, attribute.names = standards)

## ----eval = TRUE--------------------------------------------------------------
item.parameters <- results1$item.parameters
item.parameters

## ----eval = TRUE--------------------------------------------------------------
growth <- results1$growth
growth

## ----include = FALSE----------------------------------------------------------
growth.change <- growth[, 2] - growth[, 1]
growth.similar <- paste0(round(mean(growth.change[1:3]) * 100, digits = 2), "%")
growth.outlier <- paste0(round(growth.change[4] * 100, digits = 2), "%")

## ----eval = TRUE--------------------------------------------------------------
transition.probabilities <- results1$transition.probabilities
transition.probabilities

## ----include = FALSE----------------------------------------------------------
p1 <- transition.probabilities[,, 1]
p1 <- round(p1[[1,2]] * 100, digits = 2)
p1 <- paste0(p1, "%")

## ----eval = TRUE--------------------------------------------------------------
transition.posteriors <- results1$transition.posteriors
head(transition.posteriors)

## ----include = FALSE----------------------------------------------------------
maxp01 <- max(head(transition.posteriors[,,1][,2]))

## ----eval = TRUE--------------------------------------------------------------
results1$reliability

## ----eval = TRUE--------------------------------------------------------------
# Estimate TDCM with measurement invariance not assumed.
model2 <- tdcm(data, q.matrix, num.time.points = 2, invariance = FALSE)

# Compare Model 1 (longitudinal invariance assumed) to Model 2 (invariance not assumed).
tdcm.compare(model1, model2)

## ----eval = TRUE--------------------------------------------------------------
# calibrate TDCM with measurement invariance assumed, DINA measurement model
model3 <- tdcm(data, q.matrix, num.time.points = 2, rule = "DINA")

#calibrate TDCM with measurement invariance assumed, ACDM measurement model
model4 <- tdcm(data, q.matrix, num.time.points = 2, rule = "ACDM")

#compare Model 1 (full LCDM) to Model 3 (DINA)
tdcm.compare(model1, model3)

#compare Model 1 (full LCDM) to Model 4 (ACDM)
tdcm.compare(model1, model4)

## ----eval = TRUE--------------------------------------------------------------
results1$model.fit$Global.Fit.Stats
results1$model.fit$Global.Fit.Tests
results1$model.fit$Global.Fit.Stats2
results1$model.fit$Item.RMSEA
results1$model.fit$Mean.Item.RMSEA

## ----eval = FALSE-------------------------------------------------------------
#  # plot results (check plot viewer for line plot and bar chart)
#  tdcm.plot(results1, attribute.names = standards)

## ----eval = TRUE--------------------------------------------------------------
#load the TDCM library
library(TDCM)

#read data, Q-matrix, and group labels
dat4 <- data.tdcm04$data
qmat4 <- data.tdcm04$q.matrix
groups <- data.tdcm04$groups
head(dat4)

## ----eval = TRUE--------------------------------------------------------------

#calibrate mgTDCM with item and group invariance assumed, full LCDM
mg1 <- mg.tdcm(data = dat4, q.matrix = qmat4, num.time.points = 2, rule = "GDINA", groups = groups, group.invariance = TRUE, item.invariance = TRUE)


## ----eval = TRUE--------------------------------------------------------------

#calibrate mgTDCM with item invariance assumed, full LCDM
mg2 <- mg.tdcm(data = dat4, q.matrix = qmat4, num.time.points = 2, groups = groups, group.invariance = FALSE, item.invariance = TRUE)

#calibrate mgTDCM with group invariance assumed, full LCDM
mg3 <- mg.tdcm(data = dat4, q.matrix = qmat4, num.time.points = 2, groups = groups, group.invariance = TRUE, item.invariance = FALSE)

#calibrate mgTDCM with no invariance assumed, full LCDM
mg4 <- mg.tdcm(data = dat4, q.matrix = qmat4, num.time.points = 2, groups = groups, group.invariance = FALSE, item.invariance = FALSE)

#compare Model 1 (group/item invariance) to Model 2 (no group invariance)
tdcm.compare(mg1, mg2)

#compare Model 1 (group/item invariance) to Model 3 (no item invariance)
tdcm.compare(mg1, mg3)

#compare Model 1 (group/item invariance) to Model 4 (no invariance)
tdcm.compare(model1, model4)

## ----eval = TRUE--------------------------------------------------------------

#summarize results
resultsmg1 <- mg.tdcm.summary(mg1, num.time.points = 2, attribute.names = c("4.MD.1", "4.MD.2", "4.MD.3", "4.MD.4"), group.names = c("Control", "Treatment"))
resultsmg1$item.parameters
resultsmg1$growth
resultsmg1$transition.probabilities
head(resultsmg1$transition.posteriors)
resultsmg1$reliability

## ----eval = TRUE--------------------------------------------------------------

#plot results (check plot viewer for line plots and bar charts)
tdcm.plot(resultsmg1, attribute.names = c("4.MD.1", "4.MD.2", "4.MD.3", "4.MD.4"), 
          group.names = c("Control", "Treatment"))


