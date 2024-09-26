set.seed(20211123)

X <- read.csv(system.file("test_data", "parasaurolophusA.txt", package="LOMAR"), sep = "\t")
Y <- read.csv(system.file("test_data", "parasaurolophusB.txt", package="LOMAR"), sep = "\t")
Z <- read.csv(system.file("test_data", "parasaurolophusC.txt", package="LOMAR"), sep = "\t")

PS <- list(X, Y, Z)
C <- list()
for(i in 1:3) {
  cv <- diag(0.1, ncol(PS[[i]])) + jitter(0.01, amount = 0.01)
  cv <- replicate(nrow(PS[[i]]), cv)
  C[[i]] <- cv
}

test_that("Registration with cpd works", {
  trsf <- cpd(X, Y)
  expect_true(trsf$converged)
  expect_lt(trsf$sigma, 1e-5)
})

test_that("Registration with icp works", {
  trsf <- icp(X, Y, iterations = 20)
  expect_true(trsf$conv)
})

test_that("Joint registration of multiple point clouds works", {
  trsf <- jrmpc(PS, C = C, K = 100, maxIter = 20, tol = 0.01, model.selection = TRUE)
  expect_lte(trsf$iter, 20)
  expect_lte(trsf$conv, 0.025)
})

# test_that("Registration with wgmmreg works", {
#   trsf <- wgmmreg(PS[[1]], PS[[3]], C[[1]], C[[3]], maxIter = 20, tol = 1e-3)
#   expect_lte(trsf$c, 3e-3)
#   expect_true(trsf$converged)
# })
