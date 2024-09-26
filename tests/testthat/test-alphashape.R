set.seed(20230328)

X <- read.csv(system.file("test_data", "parasaurolophusA.txt", package="LOMAR"), sep = "\t")
colnames(X) <- c("x", "y", "z")

test_that("Computing 3D shape features works", {
  as <- get_shape(X, alpha = 25)
  features <- shape_features_3d(as)
  # Remove features not on the same tolerance scale for testing
  features <- features[c("major.axis", "minor.axis", "least.axis", "elongation", "flatness", "sphericity")] 
  expect_equal(features, 
               c("major.axis" = 72.026, "minor.axis" = 34.592, "least.axis" = 22.542,
                 "elongation" = 0.480, "flatness" = 0.313, "sphericity" = 0.025),
               tolerance = 1e-2)
})
