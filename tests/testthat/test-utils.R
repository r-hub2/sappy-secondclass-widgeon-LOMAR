set.seed(20211123)

data.file <- system.file("test_data", "simulated_NUP107_data.csv", package="LOMAR")
locs <- locs_from_csv(data.file)
PS <- point_sets_from_locs(locs = locs, keep.locprec = TRUE)

test_that("Cropping point set works", {
  pts <- crop_point_set(PS[[2]], c(60, 60, 120))
  expect_equal(nrow(pts), 38)
})

test_that("Localization precision can be converted to covariance matrix", {
  PS <- lapply(PS, function(p) {
    p <- cbind(p, locprecz = p[, "locprec"])
  })
  C <- locprec2cov(PS, scale = FALSE)
  expect_equal(length(C), length(PS))
  expect_equal(length(dim(C[[1]])), 3)
})

test_that("Image can be generated from localizations by histogram method", {
  I <- points2img(points = as.data.frame(PS[[2]]), voxel.size = c(10, 10, 30), method = 'histogram')
  expect_equal(dim(I), c(30, 32, 29))
})

test_that("Image can be generated from localizations by photon method", {
  pts <- as.data.frame(cbind(PS[[2]], split(locs$phot, locs$site)[[2]]))
  colnames(pts)[ncol(pts)] <- "phot"
  ncpu <- parallel::detectCores()
  if(ncpu>=2) {
    ncpu <- 2
  } else {
    ncpu <- 1
  }
  suppressWarnings(I <- points2img(points = pts, voxel.size = c(10, 10, 30), method = 'photon', ncpu = ncpu))
  expect_equal(dim(I), c(30, 32, 29))
})

test_that("Circle Hough transform works", {
  I <- points2img(points = as.data.frame(PS[[2]]), voxel.size = c(10, 10, 30), method = 'histogram')
  img <- apply(I, c(1,2), sum) # 2D projection along z
  df <- circle_hough_transform(img, rmin = 5, rmax = 5, threshold = 0.4)
  expect_equal(df[, c("x", "y", "r")], data.frame(x = 15, y = 12, r = 5))
})

test_that("coloc_index works", {
  ci <- coloc_index(P1 = standardize_coordinates(PS[[1]])$X, P2 = standardize_coordinates(PS[[2]])$X)
  ci <- sum(unlist(lapply(ci, mean)))
  expect_equal(ci, 0.650, tolerance = 0.001)
})
