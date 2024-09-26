set.seed(20211123)

data.file <- system.file("test_data", "simulated_NUP107_data.csv", package="LOMAR")
locs <- locs_from_csv(data.file)
PS <- point_sets_from_locs(locs = locs, keep.locprec = FALSE)

test_that("Persistence diagrams can be computed", {
  Diag <- get_persistence_diagrams(PS, maxdimension = 1, maxscale = 100)
  expect_equal(length(Diag), length(PS))
  expect_equal(length(which(Diag[[3]][, "dimension"] == 1)), 8)
})

test_that("Sliced Wasserstein distance matrix can be computed", {
  Diag <- get_persistence_diagrams(PS[1:20], maxdimension = 1, maxscale = 100)
  ncpu <- parallel::detectCores()
  if(ncpu>=2) {
    ncpu <- 2
  } else {
    ncpu <- 1
  }
  sWD <- get_kernel_matrix(Diag, method = "sWd", M = 10, return.dist = TRUE, ncpu = ncpu)
  dt <- determinant(sWD)
  expect_equal(as.numeric(floor(dt$modulus)), 106)
})