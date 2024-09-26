set.seed(20211123)

data.file <- system.file("test_data", "simulated_NUP107_data.csv", package="LOMAR")
locs <- locs_from_csv(data.file)

test_that("Reading from csv file works", {
  # Check that dimensions are as expected
  expect_equal(1+nrow(locs), length(readLines(data.file)))
  expect_equal(ncol(locs), length(names(locs)))
})

test_that("Getting point sets from localizations works", {
  PS <- point_sets_from_locs(locs = locs, keep.locprec = FALSE)
  # Check we get the expected number of point sets
  expect_equal(length(PS), max(locs$site))
  # Check column names
  expect_equal(colnames(PS[[1]]), c("x", "y", "z"))
})
