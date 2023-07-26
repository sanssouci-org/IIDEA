test_that("exampledata", {
  data <- exampleData(type = "microarrays")
  expect_equal(length(data), 4)
})
