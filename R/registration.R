## D: dimension of the point sets
## M, N: number of points in each point sets
## X: N x D matrix with the first point set
## Y: M x D matrix with the second point set (= the GMM centroids)
## R: D x D rotation matrix
## t: 1 x D translation vector
## s: scaling factor
## P: Posterior probabilities of the GMM centroids given the data point
## w: weight of the uniform distribution

#' apply_transformation
#' 
#' Apply rotation and translation to a point set
#'
#' @param X a point set as an N x D matrix
#' @param R D x D rotation matrix
#' @param t 1 x D translation vector
#' @param s scaling factor
#' @return transformed point set as a N x D matrix
#' @export

apply_transformation <- function(X, R, t, s) {
  return(t(apply(X, 1, function(p) s * R %*% as.matrix(p) + t)))
}

#' tr
#' 
#' Compute the trace of a matrix
#' 
#' @param x matrix
#' @return trace of the matrix
#' @export

tr <- function(x) {
  return(sum(diag(x)))
}

#' standardize_coordinates
#' 
#' Transform coordinates to have mean 0 and standard deviation 1
#' 
#' @param X point set as N x D matrix
#' @return a list of 
#'      X: standardized matrix, 
#'      mu: vector of means, 
#'      sigma: standard deviation
#' @export

standardize_coordinates <- function(X) {
  mu <- colMeans(X)
  Xs <- t(apply(X, 1, function(x) {x - mu}))
  sigma <- sqrt(sum(Xs * Xs))
  Xs <- Xs/sigma
  return(list(X = Xs, mu = mu, sigma = sigma))
}

#' restore_coordinates
#' 
#' Restore coordinates from mean 0 and standard deviation 1
#' to their original distribution
#' 
#' @param X standardized point set as N x D matrix
#' @param mu 1 x D vector of means
#' @param sigma standard deviation
#' @return N X D matrix of unstandardized coordinates
#' @export

restore_coordinates <- function(X, mu, sigma) {
  X <- X * sigma
  X <- t(apply(X, 1, function(x) {x + mu}))
  return(X)
}

#' cpd
#' 
#' Affine and rigid registration of two point sets using the coherent point drift algorithm. 
#' See: Myronenko A., Song X. (2010): "Point-Set Registration: Coherent Point Drift",
#' IEEE Trans. on Pattern Analysis and Machine Intelligence, vol. 32, issue 12, pp. 2262-2275.
#' 
#' @param X reference point set, a N x D matrix
#' @param Y point set to transform, a M x D matrix, 
#' @param w noise weight in the range [0, 1)
#' @param weights a M x N matrix of point correspondence weights
#' @param scale logical (default: FALSE), whether to use scaling
#' @param maxIter maximum number of iterations to perform (default: 100)
#' @param subsample if set, use this randomly selected fraction of the points
#' @param tol tolerance for determining convergence
#' @return a list of
#' \itemize{
#'    \item Y: transformed point set, 
#'    \item R: rotation matrix, 
#'    \item t: translation vector, 
#'    \item s: scaling factor, 
#'    \item P: matrix of correspondence probabilities between the two point sets,
#'    \item sigma: final variance,
#'    \item iter: number of iterations performed,
#'    \item converged: boolean, whether the algorithm has converged.
#' }      
#' @examples
#' data.file1 <- system.file("test_data", "parasaurolophusA.txt", package = "LOMAR",
#'  mustWork = TRUE)
#' PS1 <- read.csv(data.file1, sep = '\t', header = FALSE)
#' data.file2 <- system.file("test_data", "parasaurolophusB.txt", package = "LOMAR",
#'  mustWork = TRUE)
#' PS2 <- read.csv(data.file2, sep = '\t', header = FALSE)
#' transformation <- cpd(PS1, PS2, maxIter = 10, tol = 1e-3)
#' \dontrun{
#' # Visualize registration outcome
#' library(rgl)
#' plot3d(PS1, col = "blue")
#' points3d(PS2, col = "green")
#' points3d(transformation[['Y']], col = "magenta")
#' }
#' @export

cpd <- function(X, Y, w = 0, weights = NULL, scale = FALSE, maxIter = 100, subsample = NULL, tol = 1e-4) {
  if (dim(Y)[2]!=dim(X)[2]) {
    stop("ERROR: Point sets must have same dimension")
  }
  if (w >= 1 || w < 0) {
    stop("ERROR: w must be in the interval [0,1)")
  }
  if (!is.null(subsample) && (subsample >= 1 || subsample<=0)) {
    stop("ERROR: subsample must be a fraction between 0 and 1")
  }
  D <- dim(X)[2]

  X <- as.matrix(X)
  Xs <- standardize_coordinates(X)
  X <- Xs[["X"]]
  Y <- as.matrix(Y)
  Ys <- standardize_coordinates(Y)
  Y <- Ys[["X"]]
  
  Y0 <- NULL
  if(!is.null(subsample)) {
    sample.X <- sample(nrow(X), floor(subsample * nrow(X)))
    sample.Y <- sample(nrow(Y), floor(subsample * nrow(Y)))
    X <- X[sample.X,]
    Y0 <- Y
    Y <- Y[sample.Y,]
    if(!is.null(weights)) {
      weights <- weights[sample.Y, sample.X]
    }
  }
  M <- dim(Y)[1]
  N <- dim(X)[1]
  
  R <- diag(D)
  t <- rep(0, D)
  s <- 1
  Yt <- apply_transformation(Y, R, t, s)
  sigma <- sum(proxy::dist(X,Y)^2)/(D*M*N)
  P <- matrix(NA, nrow = M, ncol = N)
  sigma.old <- Inf
  
  steps <- 0
  converged <- FALSE
  while(steps<maxIter) {
    steps <- steps + 1
  
    ## E step. Compute P.
    noise <- (2*pi*sigma)^(0.5*D)*(w*M)/((1-w)*N)
    dist <- proxy::dist(Yt,X)^2
    if(!is.null(weights)) {
      dist <- dist / (weights+1e-9)
    }
    P <- exp(-dist/(2*sigma)) + 1e-9
    P <- apply(P, 2, function(x) {x/(sum(x) + noise)})
    
    ## M step. Compute R, t, sigma
    Np <- sum(P)
    muX <-  colSums(X * colSums(P))/Np 
    muY <-  colSums(Y * rowSums(P))/Np
    Xhat <- t(apply(X, 1, function(x) {x - muX}))
    Yhat <- t(apply(Y, 1, function(x) {x - muY}))
    
    A <- t(Xhat) %*% t(P) %*% Yhat
    svd.A <- svd(A)
    C <- diag(D)
    C[length(C)] <- det(svd.A$u) * det(svd.A$v)
    R <- svd.A$u %*% C %*% t(svd.A$v)
    if(scale) {
      s <- sum(A*R)/tr(t(Yhat) %*% (diag(rowSums(P)) %*% Yhat))
    }      
    t <- muX - s * R %*% muY

    Yt <- apply_transformation(Y, R, t, s)
    
    sigma.old <- sigma
    if(scale) {
      sigma <- (tr(Xhat %*% t(Xhat) %*% diag(colSums(P))) - s * sum(A*R))/(Np*D)
    } else {
      sigma <- (tr(Xhat %*% t(Xhat) %*% diag(colSums(P))) +
                  tr(t(Yhat) %*% (diag(rowSums(P)) %*% Yhat)) - 2 * sum(A*R))/(Np*D)
    }
    conv <- abs(log(sigma.old) - log(sigma))
    if(conv < tol) {
      converged <- TRUE
      break
    }
  }
  if(!is.null(Y0)) { # We've worked with a downsampled point set
    Yt <- apply_transformation(Y0, R, t, s)
  }
  Yt <- restore_coordinates(Yt, Xs[["mu"]], Ys[["sigma"]])
  return(list(Y = Yt, R = R, t = t*Ys[["sigma"]], s = s, P = P, sigma = sigma, iter = steps, converged = converged))
}

#' icp
#' 
#' Rigid registration of two point sets using the iterative closest point algorithm. 
#' 
#' @param X reference point set, a N x D matrix
#' @param Y point set to transform, a M x D matrix, 
#' @param weights vector of length nrow(Y) containing weights for each point in Y. Not implemented.
#' @param scale logical (default: FALSE), whether to use scaling.
#' @param subsample if set, use this randomly selected fraction of the points
#' @param iterations number of iterations to perform (default: 100)
#' @param tol tolerance for determining convergence
#' @return a list of
#' \itemize{
#'    \item Y: transformed point set, a M x D matrix,
#'    \item R: rotation matrix, 
#'    \item t: translation vector, 
#'    \item s: scaling factor, 
#'    \item iter: number of iterations performed,
#'    \item conv: boolean, whether the algorithm has converged.
#' }
#' @examples
#' data.file1 <- system.file("test_data", "parasaurolophusA.txt", package = "LOMAR",
#'  mustWork = TRUE)
#' PS1 <- read.csv(data.file1, sep = '\t', header = FALSE)
#' data.file2 <- system.file("test_data", "parasaurolophusB.txt", package = "LOMAR",
#'  mustWork = TRUE)
#' PS2 <- read.csv(data.file2, sep = '\t', header = FALSE)
#' transformation <- icp(PS1, PS2, iterations = 10, tol = 1e-3)
#' \dontrun{
#' # Visualize registration outcome
#' library(rgl)
#' plot3d(PS1, col = "blue")
#' points3d(PS2, col = "green")
#' points3d(transformation[['Y']], col = "magenta")
#' }
#' @export

icp <- function(X, Y, weights = NULL, iterations = 100, subsample = NULL, scale = FALSE, tol = 1e-3) {
  if (dim(Y)[2]!=dim(X)[2]) {
    stop("ERROR: Point sets must have same dimension")
  }
  if (!is.null(subsample) && (subsample >= 1 || subsample<=0)) {
    stop("ERROR: subsample must be a fraction between 0 and 1")
  }
  D <- dim(X)[2]
  
  X <- as.matrix(X)
  Xs <- standardize_coordinates(X)
  X <- Xs[["X"]]
  Y <- as.matrix(Y)
  Ys <- standardize_coordinates(Y)
  Y <- Ys[["X"]]
  
  Y0 <- NULL
  if(!is.null(subsample)) {
    sample.X <- sample(nrow(X), floor(subsample * nrow(X)))
    sample.Y <- sample(nrow(Y), floor(subsample * nrow(Y)))
    X <- X[sample.X,]
    Y0 <- Y
    Y <- Y[sample.Y,]
    if(!is.null(weights)) {
      weights <- weights[sample.Y, sample.X]
    }
  }
  M <- dim(Y)[1]
  N <- dim(X)[1]
  
  R <- diag(D)
  t <- rep(0, D)
  s <- 1
  Yt <- Y
  
  iter <- 0
  converged <- FALSE
  err <- 1e154
  R.final <- diag(D)
  t.final <- rep(0,D)
  s.final <- 1
  while(iter<iterations) {
    iter <- iter + 1
    
    # Find closest points in X
    nn <- RANN::nn2(X, Yt, k = 1, searchtype = "standard")
    Xc <- X[nn[["nn.idx"]], ]
    
    old.err <- err
    err <- sum((Xc-Yt)^2)
    if(abs(old.err-err)/err<tol) {
      converged <- TRUE
      break
    }
 
    # Compute R and t
    muX <-  colMeans(Xc) 
    muY <-  colMeans(Yt)
    Xhat <- as.matrix(Xc-muX)
    Yhat <- as.matrix(Yt-muY)
    
    A <- t(Xhat) %*% Yhat
    svd.A <- svd(A)
    C <- diag(D)
    C[length(C)] <- det(svd.A$u) * det(svd.A$v)
    R <- svd.A$u %*% C %*% t(svd.A$v)
    if(scale) {
      s <- sum(A*R)/sum(t(Yhat) %*% Yhat)
      s.final <- s.final * s
    }  
    t <- muX - s * R %*% muY
    
    Yt <- apply_transformation(Yt, R, t, s)
    R.final <- R %*% R.final
    t.final <- t.final + t
  }
  if(!is.null(Y0)) { # We've worked with a downsampled point set
    Yt <- apply_transformation(Y0, R.final, t.final, s.final)
  }
  Yt <- restore_coordinates(Yt, Xs[["mu"]], Ys[["sigma"]])
  return(list(Y = Yt, R = R.final, t = t.final*Ys[["sigma"]], s = s, iter = iter, conv = converged))
}

#' wgmmreg
#' 
#' Rigid registration of two point sets by minimizing the Wasserstein distance between GMMs
#' 
#' @param X reference point set, a N x D matrix
#' @param Y point set to transform, a M x D matrix,
#' @param CX array of covariance matrices for each point in X
#' @param CY array of covariance matrices for each point in Y
#' @param wx (optional) vector of mixture weights for X.
#' @param wy (optional) vector of mixture weights for Y.
#' @param maxIter maximum number of iterations to perform (default: 200)
#' @param subsample if set, use this randomly selected fraction of the points
#' @param tol tolerance for determining convergence (default: 1e-8)
#' @return a list of 
#'  \itemize{
#'     \item Y: transformed point set, 
#'     \item R: rotation matrix, 
#'     \item t: translation vector,
#'     \item c: final value of the cost function,
#'     \item converged: logical, whether the algorithm converged.
#'  }
#' @examples
#' data.file1 <- system.file("test_data", "parasaurolophusA.txt", package = "LOMAR",
#'  mustWork = TRUE)
#' PS1 <- read.csv(data.file1, sep = '\t', header = FALSE)
#' data.file2 <- system.file("test_data", "parasaurolophusB.txt", package = "LOMAR",
#'  mustWork = TRUE)
#' C1 <- diag(0.1, ncol(PS1)) + jitter(0.01, amount = 0.01)
#' C1 <- replicate(nrow(PS1),C1)
#' PS2 <- read.csv(data.file2, sep = '\t', header = FALSE)
#' C2 <- diag(0.1, ncol(PS2)) + jitter(0.01, amount = 0.01)
#' C2 <- replicate(nrow(PS2),C2)
#' transformation <- wgmmreg(PS1, PS2, C1, C2, subsample = 0.1, maxIter = 30, tol = 1e-4)
#' \dontrun{
#' # Visualize registration outcome
#' library(rgl)
#' plot3d(PS1, col = "blue")
#' points3d(PS2, col = "green")
#' points3d(transformation[['Y']], col = "magenta")
#' }
#' @export

wgmmreg <- function(X, Y, CX, CY, wx = NULL, wy = NULL, maxIter = 200, subsample = NULL, tol = 1e-8) {
  
  if (dim(Y)[2]!=dim(X)[2]) {
    stop("ERROR: Point sets must have same dimension")
  }
  if (!is.null(subsample) && (subsample >= 1 || subsample<=0)) {
    stop("ERROR: subsample must be a fraction between 0 and 1")
  }
  D <- dim(X)[2]
  
  X <- as.matrix(X)
  Xs <- standardize_coordinates(X)
  X <- Xs[["X"]]
  Y <- as.matrix(Y)
  Ys <- standardize_coordinates(Y)
  Y <- Ys[["X"]]
  
  Y0 <- NULL
  if(!is.null(subsample)) {
    Y0 <- Y
    sample.X <- sample(nrow(X), floor(subsample * nrow(X)))
    sample.Y <- sample(nrow(Y), floor(subsample * nrow(Y)))
    X <- X[sample.X,]
    CX <- CX[,,sample.X]
    Y <- Y[sample.Y,]
    CY <- CY[,,sample.Y]
  }
  D2 <- D*D
  if (D == 2) {
    Tr <- c(0,0,0) # translation vector + rotation angle
  } else {
    Tr <- c(0,0,0,0,0,0,1) # translation vector + quaternion
  }
  
  # Pre-compute sqrtm(sqrtm(CXi) %*% CYj %*% sqrtm(CXi)) for all covariance matrices
  k1 <- nrow(X)
  k2 <- nrow(Y)
  S <- array(0, dim = c(k1, k2, D, D))
  for(i in 1:k1) {
    for(j in 1:k2) {
      S1 <- CX[,,i]
      S2 <- CY[,,j]
      sqS1 <- pracma::sqrtm(S1)$B
      S[i,j,,] <- pracma::sqrtm(sqS1 %*% S2 %*% sqS1)$B
    }
  }
  if(is.null(wx)) {
    wx <- rep(1, k1)/k1
  }
  if(is.null(wy)) {
    wy <- rep(1, k2)/k2
  }
  # Registration by minimizing Wasserstein distance between GMMs
  converged = FALSE
  opt <- stats::optim(par = Tr, fn = costWd, gr = gradientWd,
                      method = "BFGS",
                      control = list(maxit = maxIter, reltol = tol),
                      X = X, Y = Y, CX = CX, CY = CY, S = S, w1 = wx, w2 = wy)
  if(opt$convergence != 0) {
    warning("Optimisation did not converge.\n")
    if(opt$convergence == 1) {
      warning("Maximum number of optimisation iterations reached.\n")
    } else {
      warning("There's been an error during optimisation.\n")
    }
  } else {
    converged = TRUE
  }
  cost <- opt$value
  # Extract rotation matrix and translation vector
  if (D == 2) {
    r <- opt$par[3] # rotation angle in radians
    R <- matrix(c(cos(r), -sin(r), sin(r), cos(r)), ncol = 2, nrow = 2, byrow = TRUE)
    t <- Tr[1:2] 
  } else {
    R <- q2r(opt$par[4:7]) # Convert quaternion to rotation matrix
    t <- opt$par[1:3]
  }
  if(!is.null(Y0)) {
    # Restore full Y after subsampling
    Y <- Y0
  }
  Yt <- apply_transformation(Y, R, t, 1)
  Yt <- restore_coordinates(Yt, Xs[["mu"]], Ys[["sigma"]])
  return(list(Y = Yt, R = R, t = t*Ys[["sigma"]], c = cost, converged = converged))
}

#' costWd
#' 
#' Objective function to minimize when using GMMs
#'
#' @param Tr Transformation vector as translation vector + rotation (angle in 2d, quaternion in 3d))
#' @param X matrix of means of first GMM (i.e. reference point set)
#' @param Y matrix of means of second GMM (i.e. moving point set)
#' @param CX array of covariance matrices of first GMM such that X[i,] has covariance matrix CX[,,i]
#' @param CY array of covariance matrices of second GMM such that Y[i,] has covariance matrix CY[,,i]
#' @param w1 (optional) vector of mixture weights of first GMM.
#' @param w2 (optional) vector of mixture weights of second GMM.
#' @param S (optional) array of pre-computed sqrtm(sqrtm(CX[,,i]) \%*\% CY[,,j] \%*\% sqrtm(CX[,,i]))
#' @return cost value
#' @export

costWd <- function(Tr, X, Y, CX, CY, w1 = NULL, w2 = NULL, S = NULL) {
  d <- dim(X)[2]
  d2 <- d*d
  if (d == 2) {
    r <- Tr[3] # rotation angle in radians
    R <- matrix(c(cos(r), -sin(r), sin(r), cos(r)), ncol = 2, nrow = 2, byrow = TRUE)
    t <- Tr[1:2] 
  } else {
    R <- q2r(Tr[4:7]) # Convert quaternion to rotation matrix
    t <- Tr[1:3]
  }
  Yt <- apply_transformation(Y, R, t, 1)
  c <- GMM_Wd(X, Yt, CX, CY, w1 = w1, w2 = w2, S = S)
  return(c[['d']])
}

#' gradientWd
#' 
#' Gradient of the objective function with respect to rotation and translation parameters
#'
#' @param Tr Transformation vector as translation vector + rotation (angle in 2d, quaternion in 3d))
#' @param X matrix of means of first GMM (i.e. reference point set)
#' @param Y matrix of means of second GMM (i.e. moving point set)
#' @param CX array of covariance matrices of first GMM such that X[i,] has covariance matrix C1[,,i]
#' @param CY array of covariance matrices of second GMM such that Y[i,] has covariance matrix C2[,,i]
#' @param w1 (optional) vector of mixture weights of first GMM.
#' @param w2 (optional) vector of mixture weights of second GMM.
#' @param S (optional) array of pre-computed sqrtm(sqrtm(CX[,,i]) \%*\% CY[,,j] \%*\% sqrtm(CX[,,i]))
#' @return gradient vector
#' @export

gradientWd <- function(Tr, X, Y, CX, CY, w1 = NULL, w2 = NULL, S = NULL) {
  k1 <- nrow(X)
  k2 <- nrow(Y)
  if(is.null(w1)) {
    w1 <- rep(1, k1)/k1
  }
  if(is.null(w2)) {
    w2 <- rep(1, k2)/k2
  }
  d <- dim(X)[2]
  d2 <- d*d
  if (d == 2) {
    r <- Tr[3] # rotation angle in radians
    R <- matrix(c(cos(r), -sin(r), sin(r), cos(r)), ncol = 2, nrow = 2, byrow = TRUE)
    dRdq <- matrix(c(-sin(r), -cos(r), cos(r), -sin(r)), ncol = 2, nrow = 2, byrow = TRUE) # dR/dr
    t <- Tr[1:2] 
  } else {
    R <- q2r(Tr[4:7]) # Convert quaternion to rotation matrix
    dRdq <- q2dr(Tr[4:7]) # dR/dq
    t <- Tr[1:3]
  }
  Yt <- apply_transformation(Y, R, t, 1)
 
  ot <- GMM_Wd(X, Yt, CX, CY, w1 = w1, w2 = w2, S = S)[['ot']]
  
  gradient <- rep(0, length(Tr))
  for(i in 1:k1) {
    for(j in 1:k2) {
      dfdR <- -2 * ot[i,j] * (matrix(X[i,], ncol = 1) %*% matrix(Y[j,], nrow = 1))
      if(d == 2) {
        dfdq <- dfdR * dRdq
        gradient[3] <- gradient[3] + sum(dfdq)
      } else {
        gradient[4] <- gradient[4] + sum(dfdR * dRdq[[1]])
        gradient[5] <- gradient[5] + sum(dfdR * dRdq[[2]])
        gradient[6] <- gradient[6] + sum(dfdR * dRdq[[3]])
        gradient[7] <- gradient[7] + sum(dfdR * dRdq[[4]])
      }
      dfdt <- ot[i,j] * as.matrix(Yt[j,] - X[i,])
      gradient[1:d] <- gradient[1:d] + dfdt
    }
  }
  return(gradient)
}

#' Gaussian_Wd
#' 
#' Compute 2-Wasserstein distance between two Gaussian distributions
#' 
#' @param m1 mean of first distribution
#' @param m2 mean of second distribution
#' @param S1 variance of first distribution
#' @param S2 variance of second distribution
#' @param S (optional) matrix of pre-computed sqrtm(sqrtm(S1) \%*\% S2 \%*\% sqrtm(S1))
#' @return distance value
#' @export

Gaussian_Wd <- function(m1, m2, S1, S2, S = NULL) {
  if(is.null(S)) {
    sqS1 <- pracma::sqrtm(S1)$B # faster than expm::sqrtm
    S <- pracma::sqrtm(sqS1 %*% S2 %*% sqS1)$B
  }
  d <- sum((m1 - m2)^2) + tr(S1+S2-2*S)
  return(d)
}

#' GMM_Wd
#' 
#' Compute 2-Wasserstein distance between two Gaussian mixture models
#' See: Delon J, Desolneux A. (2019) A Wasserstein-type distance in the space of Gaussian Mixture Models. hal-02178204v2
#' 
#' @param m1 matrix of means of first GMM
#' @param m2 matrix of means of second GMM
#' @param S1 array of covariance matrices of first GMM such that m1[i,] has covariance matrix S1[,,i]
#' @param S2 array of covariance matrices of second GMM such that m2[i,] has covariance matrix S2[,,i]
#' @param w1 (optional) vector of mixture weights of first GMM.
#' @param w2 (optional) vector of mixture weights of second GMM.
#' @param S (optional) array of pre-computed sqrtm(sqrtm(S1[,,i]) \%*\% S2[,,j] \%*\% sqrtm(S1[,,i]))
#' @return list of distance value d and optimal transport matrix ot
#' @export

GMM_Wd <- function(m1, m2, S1, S2, w1 = NULL, w2 = NULL, S = NULL) {
  k1 <- nrow(m1)
  k2 <- nrow(m2)
  if(is.null(w1)) {
    w1 <- rep(1, k1)/k1
  }
  if(is.null(w2)) {
    w2 <- rep(1, k2)/k2
  }
  # Matrix of pairwise distances between all Gaussians from the two GMMs
  M <- matrix(0, nrow = k1, ncol = k2)
  for(i in 1:k1) {
    M[i,] <- sapply(1:k2, FUN = function(j) {Gaussian_Wd(m1[i,], m2[j,], S1[,,i], S2[,,j], S = S[i,j,,])})
  }
  # Compute optimal transport plan
  ot <- transport::transport(w1, w2, M, method = "networkflow")
  # Reshape data frame into a matrix
  ot <- reshape2::acast(ot, from~to, value.var="mass")
  ot[is.na(ot)] <- 0
  d <- sum(ot * M)
  return(list(d = d, ot = ot))
}

#' Get derivative of 3D rotation matrix from quaternion
#'
#' @param q quaternion
#' @return derivative of rotation matrix

q2dr <- function(q) {
  
  a <- q[4]
  b <- q[1]
  c <- q[2]
  d <- q[3]
  
  a2 <- a * a
  b2 <- b * b
  c2 <- c * c
  d2 <- d * d
  
  z <- a2 + b2 + c2 + d2
  z2 <- z * z
  
  R <- q2r(q)
  
  dR1 <- matrix(0, nrow = 3, ncol = 3)
  dR2 <- matrix(0, nrow = 3, ncol = 3)
  dR3 <- matrix(0, nrow = 3, ncol = 3)
  dR4 <- matrix(0, nrow = 3, ncol = 3)
  
  # Diagonal terms
  dR1[1,1] <- +4 * b * (c2 + d2) / z2;
  dR2[1,1] <- -4 * c * (b2 + a2) / z2;
  dR3[1,1] <- -4 * d * (b2 + a2) / z2;
  dR4[1,1] <- +4 * a * (c2 + d2) / z2;
  
  dR1[2,2] <- -4 * b * (c2 + a2) / z2;
  dR2[2,2] <- +4 * c * (b2 + d2) / z2;
  dR3[2,2] <- -4 * d * (c2 + a2) / z2;
  dR4[2,2] <- +4 * a * (b2 + d2) / z2;
  
  dR1[3,3] <- -4 * b * (d2 + a2) / z2;
  dR2[3,3] <- -4 * c * (a2 + d2) / z2;
  dR3[3,3] <- +4 * d * (b2 + c2) / z2;
  dR4[3,3] <- +4 * a * (b2 + c2) / z2;
  
  # Off-diagonal terms
  dR1[1,2] <- +2 * c / z - 2 * b * R[1,2] / z2;
  dR2[1,2] <- +2 * b / z - 2 * c * R[1,2] / z2;
  dR3[1,2] <- -2 * a / z - 2 * d * R[1,2] / z2;
  dR4[1,2] <- -2 * d / z - 2 * a * R[1,2] / z2;
  
  dR1[1,3] <- +2 * d / z - 2 * b * R[1,3] / z2;
  dR2[1,3] <- +2 * a / z - 2 * c * R[1,3] / z2;
  dR3[1,3] <- +2 * b / z - 2 * d * R[1,3] / z2;
  dR4[1,3] <- +2 * c / z - 2 * a * R[1,3] / z2;
  
  dR1[2,1] <- +2 * c / z - 2 * b * R[2,1] / z2;
  dR2[2,1] <- +2 * b / z - 2 * c * R[2,1] / z2;
  dR3[2,1] <- +2 * a / z - 2 * d * R[2,1] / z2;
  dR4[2,1] <- +2 * d / z - 2 * a * R[2,1] / z2;
  
  dR1[2,3] <- -2 * a / z - 2 * b * R[2,3] / z2;
  dR2[2,3] <- +2 * d / z - 2 * c * R[2,3] / z2;
  dR3[2,3] <- +2 * c / z - 2 * d * R[2,3] / z2;
  dR4[2,3] <- -2 * b / z - 2 * a * R[2,3] / z2;
  
  dR1[3,1] <- +2 * d / z - 2 * b * R[3,1] / z2;
  dR2[3,1] <- -2 * a / z - 2 * c * R[3,1] / z2;
  dR3[3,1] <- +2 * b / z - 2 * d * R[3,1] / z2;
  dR4[3,1] <- -2 * c / z - 2 * a * R[3,1] / z2;
  
  dR1[3,2] <- +2 * a / z - 2 * b * R[3,2] / z2;
  dR2[3,2] <- +2 * d / z - 2 * c * R[3,2] / z2;
  dR3[3,2] <- +2 * c / z - 2 * d * R[3,2] / z2;
  dR4[3,2] <- +2 * b / z - 2 * a * R[3,2] / z2;
  
  return(list(dR1,dR2,dR3,dR4))
}

#' Convert quaternion to rotation matrix
#' http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
#'
#' @param q quaternion
#' @return rotation matrix

q2r <- function(q) {
  
  a <- q[4];
  b <- q[1];
  c <- q[2];
  d <- q[3];
  
  a2 <- a * a;
  b2 <- b * b;
  c2 <- c * c;
  d2 <- d * d;
  
  bc <- b * c;
  cd <- c * d;
  bd <- b * d;
  ab <- a * b;
  ac <- a * c;
  ad <- a * d;
  
  R <- matrix(0, nrow = 3, ncol = 3)
  
  # Diagonal terms
  R[1,1] <- a2 + b2 - c2 - d2;
  R[2,2] <- a2 - b2 + c2 - d2;
  R[3,3] <- a2 - b2 - c2 + d2;
  
  # Off-diagonal terms
  R[1,2] <- 2 * (bc - ad);
  R[1,3] <- 2 * (bd + ac);
  R[2,1] <- 2 * (bc + ad);
  R[2,3] <- 2 * (cd - ab);
  R[3,1] <- 2 * (bd - ac);
  R[3,2] <- 2 * (cd + ab);
  
  # Normalise (make rotation matrix orthogonal)
  z <- a2 + b2 + c2 + d2;
  R <- R / z;
  
  return(R)
}

#' jrmpc
#' 
#' Joint registration of multiple point sets
#' See: G.  D.  Evangelidis,  D.  Kounades-Bastian,  R.  Horaud,  andE. Z. Psarakis. 
#' A generative model for the joint registration of multiple point sets. 
#' In European Conference on Computer Vision, pages 109–122. Springer, 2014
#' 
#' @param V list of point sets as N x D matrices
#' @param C (optional) list of arrays of covariance matrices with C[[j]][,,i] the covariance matrix associated with point i of set j.
#' @param K (optional) number of components of the GMM, defaults to the average number of points in a set.
#' @param g (optional) proportion of noisy points, defaults to 1/K. If set, priors will be initialized uniformly.
#' @param initialPriors (optional) vector of length K of prior probabilities. Defaults to uniform distribution using g.
#'                      If set, will determine g so it is an error to specify g with initialPriors.
#' @param maxIter maximum number of iterations to perform (default: 100).
#' @param fixedVarIter number of iterations before starting variance updates
#' @param updatePriors logical, whether to update priors at each iteration (default: TRUE).
#' @param tol tolerance for determining convergence (default: 1e-2).
#' @param initializeBy (optional) how to initialize the GMM means. Defaults to distributing the means on the surface of the sphere enclosing all (centred) sets.
#'                     Currently supported values are:
#'                     \itemize{
#'                      \item 'sampling': sample from the data,
#'                      \item a K x D matrix of points
#'}
#' @param model.selection whether to perform model selection (default: FALSE). 
#'                        If set to TRUE, GMM components with no support in the data are deleted.
#' @param model.selection.threshold value below which we consider a GMM component has no support, set to 1/K if not explicitly given                      
#' @param rotation.only if set to TRUE, no translation is performed (default: FALSE)
#' @return a list of 
#' \itemize{
#'     \item Y: list of transformed point sets as N x d matrices,
#'     \item R: list of d x d rotation matrices, one for each point set in V,
#'     \item t: list of translation vectors, one for each point set in V,
#'     \item M: centres of the GMM,
#'     \item S: variances of the GMM.
#'     \item a: list of posterior probabilities as N x K matrices
#'     \item iter: number of iterations
#'     \item conv: error value used to evaluate convergence relative to tol
#'     \item z: support scores of the GMM components
#' }
#' @examples
#' X <- read.csv(system.file("test_data", "parasaurolophusA.txt", package="LOMAR",
#'  mustWork = TRUE), sep = "\t")
#' Y <- read.csv(system.file("test_data", "parasaurolophusB.txt", package="LOMAR",
#' mustWork = TRUE), sep = "\t")
#' Z <- read.csv(system.file("test_data", "parasaurolophusC.txt", package="LOMAR",
#' mustWork = TRUE), sep = "\t")
#' PS <- list(X, Y, Z)
#' C <- list()
#' for(i in 1:3) {
#'  cv <- diag(0.1, ncol(PS[[i]])) + jitter(0.01, amount = 0.01)
#'  cv <- replicate(nrow(PS[[i]]), cv)
#'  C[[i]] <- cv
#'}
#' transformation <- jrmpc(PS, C = C, K = 100, maxIter = 20, tol = 0.01, 
#' model.selection = TRUE)
#' \dontrun{
#' # Visualize registration outcome
#' library(rgl)
#' colours <- c("blue", "green", "magenta")
#' Yt <- transformation[['Y']]
#' plot3d(Yt[[1]], col = colours[1])
#' for(i in 2:length(Yt)) {
#'  points3d(Yt[[i]], col = colours[i])
#' }
#' # Visualize GMM centres highlighting those with high variance
#' GMM <- as.data.frame(cbind(transformation[['M']], transformation[['S']]))
#' colnames(GMM) <- c("x", "y", "z", "S")
#' colours <- rep("blue", nrow(GMM))
#' # Find high variance components
#' threshold <- quantile(transformation[['S']], 0.75)
#' high.var.idx <- which(transformation[['S']]>threshold)
#' colours[high.var.idx] <- "red"
#' plot3d(GMM[, c("x", "y", "z")], col = colours, type = 's', size = 2, box = FALSE, xlab = '', 
#'       ylab = '', zlab = '', xlim = c(-0.15,0.15), ylim = c(-0.15, 0.15), 
#'       zlim = c(-0.15, 0.15))
#' }
#' @export

jrmpc <-function(V, C = NULL, K = NULL, g = NULL, initialPriors = NULL, updatePriors = TRUE, 
                 maxIter = 100, fixedVarIter = 0, tol = 1e-2, initializeBy = NULL,
                 model.selection = FALSE, model.selection.threshold = NULL, rotation.only = FALSE) {

  M <- length(V); # Number of point sets
  d <- ncol(V[[1]]) # Dimensions of the point sets
  
  if(!is.null(initialPriors) && !is.null(g)) {
    stop("ERROR: only one of initialPriors and g should be specified.")
  }
  
  if(is.null(C)) {
    C <- list()
    for(j in 1:M) {
      C[[j]] <- replicate(nrow(V[[j]]),  diag(1, d))
    }
  }

  ####################
  ## Initialization ##
  ####################
  
  ## Centre on origin
  Vs <- list()
  for(j in 1:M) {
    Vs[[j]] <- standardize_coordinates(V[[j]])
    V[[j]] <- Vs[[j]][["X"]]
  }
  # Radius of enclosing sphere
  radius <- max(unlist(lapply(V, function(v) {as.vector(proxy::dist(v))})))/2
  ## Number of componemts in the GMM
  if(is.null(K)) {
    K <- round(mean(sapply(V, nrow))) # Mean number of points in a set
    if(model.selection) { 
      # Make K large enough
      K <- 2 * K
    }
  } else if(!is.null(initializeBy) && !is.null(nrow(initializeBy))) {
    K <- nrow(initializeBy)
  }
  ## Initial GMM means
  X <- matrix(NA, nrow = K, ncol = d)
  if(!is.null(initializeBy) && !is.null(nrow(initializeBy))) {
    # Use given matrix of centres
    if(nrow(initializeBy) != K || ncol(initializeBy) != d) {
      stop("ERROR: Incorrect dimensions for initializeBy matrix of centres.")
    }
    X <- initializeBy
  } else if(!is.null(initializeBy) && initializeBy == 'sampling') {
    ## Sample from the data
    idx <- sample(1:M, K, replace = TRUE)
    for(i in 1:K) {
      k <- sample(1:nrow(V[[idx[i]]]), 1)
      X[i,] <- V[[idx[i]]][k,]
    }
  } else {
    ## Distribute on the enclosing sphere
    for(i in 1:K) {
      X[i,] <- stats::rnorm(d)
      l <- sqrt(sum(X[i,]^2))
      X[i,] <- (X[i,]/l) * radius
    }
  }
  if(model.selection && is.null(model.selection.threshold)) {
    model.selection.threshold <- 1/K
  }
  
  ## Initial rotations and translations
  R <- lapply(1:M, function(i) { diag(d) })
  if(rotation.only) {
    t <- rep(list(rep(0, d)), M)
  } else {
    t <- lapply(V, function(v) {colMeans(X) - colMeans(v)})
  }
  
  TV <- list()
  for(j in 1:M) {
    TV[[j]] <- t(apply(V[[j]], 1, function(p) {R[[j]] %*% as.matrix(p) + t[[j]]}))
  }

  ## Initial GMM variances S as median distance between X and all the points in V
  dist <- unlist(lapply(TV, function(y) {as.vector(proxy::dist(y, X)^2)}))
  S <- rep(stats::median(dist), K)
  
  if(is.null(initialPriors)) {
    p <- rep(1/(K+1),K)
  } else {
    p <- initialPriors
    g <- (1-sum(p))/sum(p)
  }
  
  if(is.null(g)) {
    g <- (1-sum(p))/sum(p)
  } else {
    p <- rep(1/(K*(g+1)),K)
  }
  
  N <- sum(sapply(V, nrow))
  
  ## Parameter h should be proportional to the volume that
  ## encompasses all the point sets.
  ## Set to it the volume of the enclosing sphere.
  if(d == 2) {
    h = pi*radius^2
  } else {
    h <- (4/3)*pi/(radius^3)
  }
  
  u <- g/(h*(g+1)) # noise component
  
  alpha <- list() # posterior probabilities
  lambda <- list()
  W <- list() # Mixture weights
  for(j in 1:M) {
    n <- nrow(TV[[j]])
    alpha[[j]] <- matrix(0, nrow = n, ncol = K)
    W[[j]] <- matrix(NA, nrow = K, ncol = d)
  }
  
  P <- list()
  e <- matrix(1, nrow = K, ncol = 1)
  
  #########################
  ##     EM algorithm    ##
  #########################
  
  for (iter in 1:maxIter) {
    
    ## E-step: posterior probabilities
    for(j in 1:M) {
      n <- nrow(TV[[j]])
      for(i in 1:n) {
        sumk <- 0
        sqe <- stats::mahalanobis(X, TV[[j]][i,], C[[j]][,,i])
        for(k in 1:K) {
          alpha[[j]][i,k] <- (p[k]*S[k]^(-1.5))*exp(-0.5*sqe[k]/S[k])
          sumk <- sumk + alpha[[j]][i,k]
        }
        alpha[[j]][i,] <- alpha[[j]][i,]/(sumk+u) # Normalization
      }
    }
    
    ## M-step: rigid transformation
    W <- list()
    R <- list()
    t <- list()
    for(j in 1:M) {
      n <- nrow(V[[j]])
      lambda[[j]] <- as.matrix(colSums(alpha[[j]]))
      W[[j]] <- t(alpha[[j]]) %*% V[[j]]
      W[[j]] <- apply(W[[j]],2, function(x) { x/S})
      mW <- as.matrix(colSums(W[[j]]))
      b <- lambda[[j]]/S
      mX <- t(X) %*% b
      sumOfWeights <- drop(t(lambda[[j]]) %*% (1/S))
      P <- (t(X) %*% W[[j]]) - (mX %*% t(mW))/sumOfWeights
      svd.A <- svd(P)
      tmp <- diag(d)
      tmp[length(tmp)] <- det(svd.A$u) * det(svd.A$v)
      R[[j]] <- svd.A$u %*% tmp %*% t(svd.A$v)
      if(rotation.only) {
        t[[j]] <- rep(0,d)
      } else {
        t[[j]] <- (mX - R[[j]] %*% mW)/sumOfWeights
      }
    }
    
    ## M-step: Update GMMs
    
    ## Transform point sets
    for(j in 1:M) {
      TV[[j]] <- t(apply(V[[j]], 1, function(p) { R[[j]] %*% as.matrix(p) + t[[j]]}))
    }
    
    ## Update GMM means
    X.old <- X
    tmp <- do.call(cbind, lambda)
    den <- rowSums(tmp)
    X <- list()
    for(j in 1:M) {
      X[[j]] <- t(alpha[[j]]) %*% TV[[j]]
    }
    X <- Reduce("+", X)
    X <- X/den
    
    ## Update GMM variances S
    if(iter>fixedVarIter) {
      S <- list()
      for(j in 1:M) {
        sqe <- matrix(NA, nrow = nrow(TV[[j]]), ncol = K)
        for(i in 1:nrow(TV[[j]])) {
          sqe[i,] <- stats::mahalanobis(X, TV[[j]][i,], C[[j]][,,i])
        }
        S[[j]] <- alpha[[j]] * sqe
        S[[j]] <- colSums(S[[j]])
      }
      S <- Reduce("+", S)/(3*den) + 1e-9
    }
    
    ## Update priors
    if(updatePriors) {
      p <- den/((g+1)*sum(den))
    }
    
    ## Check for convergence 
    conv <- sum(abs(X.old - X))
    if(conv < tol) {
      converged <- TRUE
      break
    }
    
    ## Compute support for the GMM components
    ## using minimum message length approach as in
    ## M. A. T. Figueiredo and A. K. Jain, “Unsupervised learning of finite
    ## mixture models,” IEEE Transactions on Pattern Analysis and Machine
    ## Intelligence, vol. 24, no. 3, pp. 381–396, 2002.
    z <- 0
    for(j in 1:M) {
      z <- z + lambda[[j]]
    }
    z <- z - d
    z[z<0] <- 0
    z <- z/sum(z)
    if(model.selection) {
      ## Check if we should remove components from the GMM
      comp.to.remove <- which(z<=model.selection.threshold)
      n.to.remove <- length(comp.to.remove)
      if(n.to.remove>0 && n.to.remove<K) { # No support for some components, remove them
        K <- K - n.to.remove
        X <- X[-comp.to.remove,]
        S <- S[-comp.to.remove]
        for(j in 1:M) {
          alpha[[j]] <- alpha[[j]][,-comp.to.remove] 
        }
      }
    }    
    
  }
  for(j in 1:M) {
    TV[[j]] <- restore_coordinates(TV[[j]], rep(0,d), Vs[[j]][["sigma"]])
    t[[j]] <- t[[j]] * Vs[[j]][["sigma"]]
  }
  return(list(Y = TV, R = R, t = t, M = X, S = S, a = alpha, iter = iter, conv = conv, z = z))
}

#' multiple_registration
#' 
#' Registration of multiple point sets using tree-guided progressive registration 
#' followed by iterative refinement.
#'
#' @param PS list of point sets
#' @param registration pairwise registration method to use
#' @param refine.iter Maximum number of refinement iterations (default: 20)
#' @param ... additional arguments to the registration method
#' @return  a list of 
#' \itemize{
#'     \item Y: list of transformed point sets as N x d matrices
#' }
#' @export

multiple_registration <- function(PS, registration, refine.iter = 20, ...) {
  
  registration.name <- deparse(substitute(registration))
  
  params <- list(...)
  # Select extra parameters used by registration function
  params.idx <- which(names(params) %in% names(as.list(args(registration))))
  if(length(params.idx)>0) {
    params <- params[params.idx]
  } else {
    params <- NULL
  }

  if(registration.name == 'wgmmreg' && !methods::hasArg(C)) {
    stop("Provide list of covariance matrices as argument C when using wgmmreg")
  }
  
  if(length(params$C) > 0) {
   C <- params$C 
  } else {
    C <- NULL
  }
  n <- length(PS); # Number of point sets
  
  # Build distance matrix
  D <- matrix(Inf, nrow = n, ncol = n)
  for(i in 1:n) {
    X <- standardize_coordinates(PS[[i]])$X
    for(j in 1:n) {
      if (i == j) { next }
      Y <- standardize_coordinates(PS[[j]])$X
      nn <- RANN::nn2(X, Y, k = 1, searchtype = "standard")
      Xc <- X[nn[["nn.idx"]], ]
      D[i, j] <- mean((Xc - Y)^2)
    }
  }
  # Progressive registration using guide tree
  dend <- stats::hclust(stats::as.dist(D), method = 'complete')
  PSc <- list()
  Cc <- list()
  CX <- NULL
  CY <- NULL
  for(i in 1:nrow(dend$merge)) {
    r <- dend$merge[i,]
    if(r[1]<0) {
      X <- PS[[abs(r[1])]]
      if(length(C) > 0) {
        CX <- C[[abs(r[1])]]
      }
    } else {
      X <- PSc[[r[1]]]
      if(length(C) > 0) {
        CX <- Cc[[r[1]]]
      }
    }
    if(r[2]<0) {
      Y <- PS[[abs(r[2])]]
      if(length(C) > 0) {
        CY <- C[[abs(r[2])]]
      }
    } else {
      Y <- PSc[[r[2]]]
      if(length(C)> 0) {
        CY <- Cc[[r[2]]]
      }
    }
    if(registration.name == 'wgmmreg') {
      Yr <- registration(X = X, Y = Y, 
                         CX = CX, CY = CY,
                         subsample = params$subsample)$Y
    } else {
      Yr <- registration(X, Y, ...)$Y
    }
    PSc[[i]] <- rbind(X, Yr)
    Cc[[i]] <- abind::abind(CX, CY, along = 3)
  }
  PSr <- list()
  for(i in 1:n) {
    if(registration.name == 'wgmmreg') {
      PSr[[i]] <- registration(X = PSc[[(n-1)]], Y = PS[[i]], 
                               CX = Cc[[(n-1)]], CY = C[[i]],
                               subsample = params$subsample)$Y
    } else {
      PSr[[i]] <- registration(X = PSc[[(n-1)]], Y = PS[[i]], ...)$Y
    }
  }
  
  # Refinement
  for(i in 1:refine.iter) {
    for(j in 1:n) {
      # Remove point set j from global registration 
      tmp <- do.call(rbind, PSr[-j])
      n0 <- floor(nrow(tmp)/(n-1))
      tmp <- downsample(tmp, n = n0, k = ceiling(sqrt(n0)))
      # Register point set j to the rest
      if(registration.name == 'wgmmreg') {
        Ctmp <- abind::abind(C[-j], along = 3)
        PSr[[j]] <- registration(X = tmp, Y = PSr[[j]], 
                                 CX = Ctmp, CY = C[[j]],
                                 subsample = params$subsample)$Y
      } else {
        PSr[[j]] <- registration(tmp, PSr[[j]], ...)$Y
      }
    }
  }
  
  return(list(Y = PSr))
}
