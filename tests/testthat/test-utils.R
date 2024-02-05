################################################################################
# tdcm_emit
################################################################################

test_that("tdcm_emit produces the correct output", {
  expect_message(tdcm_emit("hello"), "[tdcm] INFO: hello", fixed = TRUE)
  expect_message(tdcm_emit("hello", func = base::message), "[tdcm] INFO: hello", fixed = TRUE)
  expect_warning(tdcm_emit("hello", func = base::warning), "[tdcm] INFO: hello", fixed = TRUE)
  expect_error(tdcm_emit("hello", func = base::stop), "[tdcm] INFO: hello", fixed = TRUE)
})

################################################################################
# tdcm_warn
################################################################################

test_that("tdcm_warn produces the correct output", {
  expect_warning(tdcm_warn("hello"), "[tdcm] WARN: hello", fixed = TRUE)
})

################################################################################
# tdcm_stop
################################################################################

test_that("tdcm_stop produces the correct output", {
  expect_error(tdcm_stop("hello"), "[tdcm] STOP: hello", fixed = TRUE)
})
