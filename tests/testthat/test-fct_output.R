test_that("exampledata", {
  data <- exampleData()
  expect_equal(length(data), 4)
})
