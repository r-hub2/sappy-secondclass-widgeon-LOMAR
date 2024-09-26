
#' binning
#'
#' Binning in 1D, 2D or 3D.
#' 
#' Copied from package aws which is no longer in CRAN.
#' Original author:  Joerg Polzehl (polzehl@wias-berlin.de)
#' who adapted code of function binning in package sm.
#'
#' @param x design matrix, dimension n x d with d in 1:3.
#' @param y either a response vector of length n or NULL.
#' @param nbins vector of length d containing number of bins for each dimension, may be set to NULL.
#' @param xrange range for endpoints of bins for each dimension, either matrix of dimension 2 x d or NULL. xrange is increased if the cube defined does not contain all design points.
#' @return a list with elements:
#'   \itemize{
#'     \item x matrix of coordinates of non-empty bin centers
#'     \item x.freq number of observations in nonempty bins
#'     \item midpoints.x1 bin centers in dimension 1
#'     \item midpoints.x2 bin centers in dimension 2
#'     \item midpoints.x3 bin centers in dimension 3
#'     \item breaks.x1 break points dimension 1
#'     \item breaks.x2 break points dimension 2
#'     \item breaks.x3 break points dimension 3
#'     \item table.freq number of observations per bin
#'     \item means means of y in non-empty bins (if y isn't NULL)
#'     \item devs standard deviations of y in non-empty bins (if y isn't NULL)
#'   }
#' @export

binning <- function (x, y, nbins, xrange = NULL) {
  if(any(nbins<2)) stop("binning - need at least 2 bins")
  dx <- dim(x)
  if (is.null(dx))
    d <- 1
  else
    d <- dx[2]
  if (d > 3) {
    warning("Binning only implemented in 1D, 2D and 3D")
    return(NULL)
  }
  if (length(nbins) < d || any(nbins < 2)) {
    warning("Invalid values for nbins")
    return(NULL)
  }
  if (!is.null(y) && length(y) * d != length(x)) {
    warning("Dimensions of design matrix incompatible with length of response vector")
    return(NULL)
  }
  if (is.null(xrange)) {
    xrange <- if (d == 1)
      range(x)
    else
      apply(x, 2, range)
  } else {
    if ((d == 1 &&
         length(xrange) != 2) || (d > 1 && any(dim(xrange) != c(2, d)))) {
      warning("Dimensions of xrange incorrect ")
      return(NULL)
    }
    xrange <-
      if (d == 1)
        range(x, xrange)
    else
      apply(rbind(x, xrange), 2, range)
  }
  xnames <- if (d > 1)
    dimnames(x)[[2]]
  else
    names(x)
  breaks.x1 <- seq(xrange[1], xrange[2], length = nbins[1] + 1)
  if (d > 1)
    breaks.x2 <- seq(xrange[1, 2], xrange[2, 2], length = nbins[2] + 1)
  if (d > 2)
    breaks.x3 <- seq(xrange[1, 3], xrange[2, 3], length = nbins[3] + 1)
  f1 <- cut(if (d == 1)
    x
    else
      x[, 1], breaks = breaks.x1)
  if (d > 1)
    f2 <- cut(x[, 2], breaks = breaks.x2)
  if (d > 2)
    f3 <- cut(x[, 3], breaks = breaks.x3)
  freq <- switch(d, table(f1), table(f1, f2), table(f1, f2, f3))
  dimnames(freq) <- NULL
  midpoints.x1 <-
    (breaks.x1[-1] + breaks.x1[-(nbins[1] + 1)]) / 2
  if (d > 1)
    midpoints.x2 <- (breaks.x2[-1] + breaks.x2[-(nbins[2] + 1)]) / 2
  if (d > 2)
    midpoints.x3 <- (breaks.x3[-1] + breaks.x3[-(nbins[3] + 1)]) / 2
  z1 <- midpoints.x1
  if (d > 1)
    z2 <- midpoints.x2
  if (d > 2)
    z3 <- midpoints.x3
  X <- switch(d, z1,
              cbind(rep(z1, length(z2)),
                    rep(z2, rep(
                      length(z1), length(z2)
                    ))),
              cbind(rep(z1, length(z2) * length(z3)),
                    rep(z2, rep(
                      length(z1) * length(z3), length(z2)
                    )),
                    rep(z3, rep(
                      length(z1) * length(z2), length(z3)
                    ))))
  X.f <- as.vector(freq)
  id <- (X.f > 0)
  if (d > 1)
    X <- X[id,]
  else
    X <- X[id]
  if (d > 1)
    dimnames(X) <- list(NULL, xnames)
  else
    names(X) <- xnames
  X.f <- X.f[id]
  result <- list(
    x = X,
    x.freq = X.f,
    midpoints.x1 = midpoints.x1,
    midpoints.x2 = if (d > 1)
      midpoints.x2
    else
      NULL,
    midpoints.x3 = if (d > 2)
      midpoints.x3
    else
      NULL,
    breaks.x1 = breaks.x1,
    breaks.x2 = if (d > 1)
      breaks.x2
    else
      NULL,
    breaks.x3 = if (d > 2)
      breaks.x3
    else
      NULL,
    table.freq = freq
  )
  if (!is.null(y) && !all(is.na(y))) {
    result$means <- as.numeric(tapply(y, switch(
      d, list(f1),
      list(f1, f2), list(f1, f2, f3)
    ),
    mean))[id]
    result$devs <- as.numeric(tapply(y, switch(
      d, list(f1),
      list(f1, f2), list(f1, f2, f3)
    ),
    function(x)
      sum((x - mean(
        x
      )) ^ 2)))[id]
  }
  result
}

#' crop_point_set
#'
#' Retain points in the set that are within the given distance from the geometric median of the set.
#' Using the geometric median is more robust than using the centre of mass (i.e. mean).
#'
#' @param point.set a point set as a matrix with columns x,y,z.
#' @param size vector of distances from the target region centre along each axis. 
#'             Points are discarded if they are outside the ellipsoid defined by size and centred on 
#'             the given position.
#' @param center (optional) coordinates of the centre of the target region. If not given, 
#'               default to the geometric median of the point set.
#' @return point set as a matrix with columns x,y,z.
#' @export

crop_point_set <- function(point.set, size, center = NULL) {
    xyz.idx <- match(colnames(point.set), c("x", "y", "z"))
    xyz.idx <- xyz.idx[!is.na(xyz.idx)]
    if(is.null(center)) {
      center <- pracma::geo_median(point.set[, xyz.idx])$p
      center <- center[names(center) %in% c('x','y','z')]
    }
    if(length(size) == 1) { size <- rep(size, length(xyz.idx))}
    if(length(xyz.idx)>length(size)) {
        stop("ERROR: point set and crop size dimensionalities don't match.")
    }
    # Center the points
    point.set[,xyz.idx] <- t(apply(point.set[, xyz.idx, drop = FALSE], 1, function(x) {x-center}))
    # Remove points not within the ellipsoid defined by size
    idx.to.remove <- which(apply(point.set[, xyz.idx, drop = FALSE], 1, function(x) {sum(x^2/size^2)>1}))
    if(length(idx.to.remove>0)) {
        point.set <- point.set[-idx.to.remove,, drop = FALSE]
    }
    # Restore original coordinates
    point.set[,xyz.idx] <- t(apply(point.set[, xyz.idx, drop = FALSE], 1, function(x) {x+center}))
    return(point.set)
}

#' ps2ary
#'
#' Convert a list of 3d point sets to a 4d array.
#' Also works for 2d point sets to 3d array conversion.
#'
#' @param point.sets a list of point sets.
#' @param dims vector of dimensions of the axes (x,y in 2d, x,y,z in 3d).
#' @return a 3d or 4d array.
#' @export

ps2ary <- function(point.sets, dims) {
    if(!(length(dims) %in% c(2,3))) {
        stop("ERROR: Dimensions can only be 2d or 3d.")
    }
    if(typeof(point.sets) != "list") {stop("ERROR: Argument point.sets must be a list.")}
    n <- length(point.sets)
    ts <- array(0, dim = c(dims, n))
    # Translate so that we start at coordinate 1 along each axis
    for(i in 1:n) {
        mins <- apply(point.sets[[i]],2, min)
        point.sets[[i]] <- t(apply(point.sets[[i]], 1, function(x) {x - mins}))
    }
    for(i in 1:n) {
        for(j in 1:nrow(point.sets[[i]])) {
            a <- round(point.sets[[i]][j,]) # round in case we have interpolated coordinates (e.g. from registration)
            if(length(dims) == 3) {
                ts[a[1], a[2], a[3], i] <- ts[a[1], a[2], a[3], i] + 1
            } else {
                ts[a[1], a[2], i] <- ts[a[1], a[2], i] + 1
            }
        }
    }
    return(ts)
}

#' ary2ps
#'
#' Convert a 4d array to a list of 3d point sets. 
#' The points are formed by extracting the coordinates of array values strictly above the given cut-off (default 0).
#'
#' @param ary a 4d array with last dimension indexing instances.
#' @param bkg Extract points for array values strictly above this (default = 0)
#' @return a list of point sets.
#' @export

ary2ps <- function(ary, bkg = 0) {
    PS <- list()
    n <- dim(ary)[4]
    ## Form point sets
    for (j in 1:n) {
        PS[[j]] <- as.matrix(which((ary[,,,j]>bkg), arr.ind = TRUE))
    }
    # Rename columns as x,y,z
    PS <- lapply(seq(PS), function(i) { X <- as.matrix(PS[[i]]); colnames(X) <- c("x","y","z"); return(X)})
    return(PS)
}

#' locprec2cov
#'
#' Converts localization precision columns to a list of arrays of covariance matrices
#'
#' @param point.sets a list of n point sets with locprec columns (locprecz column required for 3D data)
#' @param scale logical, whether to scale the localization precision by the variance of the coordinates
#' @return a list of 2x2xn or 3x3xn arrays.
#' @export

locprec2cov <- function(point.sets, scale = FALSE) {
  C <- list()
  d <- 3 # Assume 3D data
  # 3D data require locprecz column
  # If it's not in the first point set, consider we have 2D data
  if(!("locprecz" %in% colnames(point.sets[[1]]))) { d <- 2 }
  for(i in 1:length(point.sets)) {
    ns <- nrow(point.sets[[i]])
    C[[i]] <- array(0, dim = c(d,d,ns))
    if(scale) {
      v <- apply(point.sets[[i]], 2, stats::var)
    }
    for(j in 1:ns) {
      if (d==3) {
        # locprec applies to x and y
        C[[i]][,,j] <- diag(point.sets[[i]][j, c("locprec", "locprec", "locprecz")], 3)
        if(scale) {
          C[[i]][,,j] <- C[[i]][,,j] / v[c("x", "y", "z")]
        }
      } else {
        C[[i]][,,j] <- diag(point.sets[[i]][j, c("locprec", "locprec")], 2)
        if(scale) {
          C[[i]][,,j] <- C[[i]][,,j] / v[c("x", "y")]
        }
      }
    }
  }
  return(C)
}

#' downsample
#'
#' Weighted downsampling of a point set.
#' If point weights are not provided, they are computed to be proportional to the local density
#' around each point.
#'
#' @param point.set a point set
#' @param n integer, sample size. 
#' @param k integer, number of nearest neighbours to consider to estimate local density
#' @param weights a vector of probability weights
#' @return a point set
#' @export

downsample <- function(point.set, n = NULL, k = NULL, weights = NULL) {
  if(is.null(k) && is.null(weights)) {
    stop("One of k or weights should be provided")
  }
  if(is.null(n)) {
    stop("Sample size must be provided")
  }
  if(is.null(weights)) {
    weights <- local_densities(point.set, k)
    if(max(weights) == min(weights)) {
      weights <- rep((1/nrow(point.set)), nrow(point.set))
    } else {
      weights <- (weights - min(weights))/(max(weights) - min(weights))  # Rescale to [0,1]
    }
  }
  replace <- FALSE
  if(nrow(point.set) <= n) {
    message("Sample size is greater than the number of points, drawing with replacement to reach sample size.")
    replace <- TRUE
  }
  sample.idx <- sample(nrow(point.set), n, prob = weights, replace = replace)
  return(point.set[sample.idx,])
}

#' Circle Hough transform
#'
#' Extract coordinates of the centres of circles from a 2D image using the Hough transform
#'
#' @importFrom foreach %dopar%
#' @param pixels input data, either a matrix representing a 2D image or a data frame of signal coordinates with columns x, y.
#'               For images, background is expected to be 0 and signal to have positive values.
#' @param rmin minimum search radius.
#' @param rmax maximum search radius.
#' @param resolution number of steps in the circle transform (default: 360). This represents the maximum number of votes a point can get.
#' @param threshold score threshold between 0 and 1.
#' @param min.separation distance between circle centres below which overlapping circles are considered the same and merged (default to 0.25*rmin)
#' @param ncpu number of threads to use to speed up computation (default: 1)
#' @return a data frame with columns x, y, r and score
#' @examples
#' point.set <- data.frame(x = c(-9.8,-5.2,12.5,2.5,4.5,1.3,-0.2,0.4,9.3,-1.4,0.5,-1.1,-7.7),
#'                         y = c(-4.2,1.5,-0.5,12,-3,-7.2,10.9,6.7,-1.3,10,6.7,-6.2,2.9))
#' circles <- circle_hough_transform(pixels = point.set, rmin = 3, rmax = 6, resolution = 100,
#'                                   threshold = 0.1, ncpu = 1)
#' @export

circle_hough_transform <- function(pixels, rmin, rmax, threshold, resolution = 360, min.separation = rmin/4, ncpu = 1) {
  
  if (threshold<=0 || threshold>1) {
    stop("ERROR: Threshold must be between 0 and 1")
  }
  
  circles.found <- NULL
  d <- vector("numeric", 2) # dimensions of the image or region to consider
  
  if(methods::is(pixels, "matrix") && length(dim(pixels))==2) { 
    ## Extract table of coordinates of non-zero pixels
    coords <- which(pixels!=0, arr.ind = TRUE)
    if(length(coords)==0) {
      stop("No pixel with intensity above 0 was found in the image.")
    }
    pixels <- as.data.frame(cbind(coords, pixels[coords]))
    colnames(pixels) <- c("x","y","value")
  } 
  if (all(c("x","y") %in% colnames(pixels))) {
    range.x <- range(pixels$x)
    range.y <- range(pixels$y)
    d[1] <- ceiling(range.x[2]-range.x[1])
    d[2] <- ceiling(range.y[2]-range.y[1])
  } else {
    stop("ERROR: unrecognized input data")
  }
  
  radii <- seq(from = rmin, to = rmax, by = 1)
  
  cluster <- parallel::makeCluster(ncpu)
  doParallel::registerDoParallel(cluster)
  on.exit(parallel::stopCluster(cluster))
  
  # Hough transform
  A <- ff::ff(0, dim = c(d[1], d[2], length(radii))) # array(0, c(d[1], d[2], length(radii)))
  tmp <- foreach::foreach(k = 1:length(radii), .combine = 'c', .packages = c('ff')) %dopar% {
    r <- radii[k]
    for(i in 1:nrow(pixels)) {
      x <- pixels$x[i]
      y <- pixels$y[i]
      for(theta in seq(from = 0, to = 2*pi, length = resolution)) {
        a <- round(x + r * cos(theta))
        b <- round(y + r * sin(theta))
        if(a > 0 && a <= d[1] && b > 0 && b <= d[2]) {
          A[a, b, k] <- A[a, b, k] + 1
        }
      }
    }
    NULL
  }
  
  ## Find local maxima in accumulator
  ## Iteratively search for the global maximum and
  ## remove the corresponding circle
  ## Repeat until we hit the threshold
  done <- FALSE
  while(!done) {
    value.max <- max(A[])
    score <- value.max/resolution
    if (score <= threshold) {
      done <- TRUE
    } else {
      coords <- which(A[] == value.max, arr.ind = TRUE)
      ## Note: There can be more than one point with max value
      ## Remove the corresponding circles from the accumulator
      for (k in 1:nrow(coords)) {
        r.idx <- coords[k, 3]
        r <- radii[r.idx]
        coords[k, 3] <- r
        x <- coords[k, 1]
        y <- coords[k, 2]
        x0 <- max(c(1, x-r))
        y0 <- max(c(1, y-r))
        x1 <- min(c(d[1], x+r))
        y1 <- min(c(d[2], y+r))
        A[x0:x1, y0:y1,] <- 0
      }
      if(nrow(coords)>1) {
        ## Check if the centres are close to each other
        ## and merge those that are too close
        dst <- as.matrix(stats::dist(coords[,-3]))
        diag(dst) <- NA
        rws <- unique(as.vector(which(dst<min.separation, arr.ind = TRUE)))
        if(length(rws)>1) {
          new.centre <- round(apply(coords[rws,], 2, mean))
          coords <- coords[-rws,]
          coords <- rbind(coords, new.centre)
        }
      }
      circles <- cbind(coords, score)
      circles.found <- rbind(circles.found, circles)
    }
  }
  if (length(circles.found)>0) {
    colnames(circles.found) <- c("x", "y", "r", "score")
  }
  return(as.data.frame(circles.found))
}

#' points2img
#'
#' Convert a data frame of point coordinates into an image.
#' Expected photon count at each voxel is computed as in:
#' F. Huang, S. L. Schwartz, J. M. Byars, and K. A. Lidke, “Simultaneous multiple-emitter fitting for single
#' molecule super-resolution imaging,” Biomed. Opt. Express 2(5), 1377–1393 (2011).
#'
#' @importFrom foreach %dopar%
#' @param points a point set as a data frame of coordinates with columns x,y,z.
#' @param channels vector of channels to consider, must be values present in the input data frame channel column
#' @param voxel.size a numeric vector of length 3 indicating the size of the voxel along x,y and z in the same unit as the coordinates (e.g. nm)
#' @param method how to calculate voxel values. Available methods are:
#'   \itemize{
#'     \item 'histogram': value is the number of points (i.e. emitters) in the voxel
#'     \item 'photon': value is the expected number of photons from the points in the voxel. Input data frame must have columns locprec, locprecz and phot[on].
#'}
#' @param ncpu number of threads to use to speed up computation (default: 1)
#' @return an array of dimensions x,y,z and channels if applicable
#' @examples 
#' point.set <- data.frame(x = c(-9.8,-5.2,12.5,2.5,4.5,1.3,-0.2,0.4,9.3,-1.4,0.5,-1.1,-7.7),
#'                         y = c(-4.2,1.5,-0.5,12,-3,-7.2,10.9,6.7,-1.3,10,6.7,-6.2,2.9),
#'                         z = c(3.4,-3.8,-1.4,1.8,3.5,2.5,2.6,-4.8,-3.8,3.9,4.1,-3.6,-4))
#' img <- points2img(point.set, voxel.size = c(2,2,2), method = 'histogram')
#' @export

points2img <- function(points, voxel.size, method, channels = NULL, ncpu = 1) {
  
  if (!all(c("x","y","z") %in% colnames(points))) {
    stop("ERROR: Columns x,y,z must be present in input data frame.")
  }
  
  if(!is.null(channels)) {
    if(!("channel" %in% colnames(points))) {
      stop("ERROR: channel column must be present in input data frame.")
    } else if(!all(channels %in% unique(points$channel))) {
      stop("ERROR: Unknown channel requested.")
    }
  }  
  range.x <- range(points$x)
  range.y <- range(points$y)
  range.z <- range(points$z)
  image.size <- ceiling(c((range.x[2]-range.x[1])/voxel.size[1], 
                          (range.y[2]-range.y[1])/voxel.size[2], 
                          (range.z[2]-range.z[1])/voxel.size[3]))
  if(any(image.size>.Machine$integer.max)) {
    stop("ERROR: Target image size too big.")
  }
  # Define voxels
  x.bins <- cut(points$x, breaks = image.size[1])
  y.bins <- cut(points$y, breaks = image.size[2])
  z.bins <- cut(points$z, breaks = image.size[3])
  # Assign points to voxels
  points$x.bin <- as.numeric(x.bins)
  points$y.bin <- as.numeric(y.bins)
  points$z.bin <- as.numeric(z.bins)
  # Extract voxels boundaries
  x.bins <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", levels(x.bins)) ),
                  upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(x.bins))))
  y.bins <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", levels(y.bins)) ),
                  upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(y.bins))))
  z.bins <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", levels(z.bins)) ),
                  upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(z.bins))))
  x.bins.breaks <- sort(unique(c(x.bins[,"lower"], x.bins[,"upper"])))
  y.bins.breaks <- sort(unique(c(y.bins[,"lower"], y.bins[,"upper"])))
  z.bins.breaks <- sort(unique(c(z.bins[,"lower"], z.bins[,"upper"])))
  
  # Voxel position is the centre of the corresponding bin
  x.bins.centres <- apply(x.bins, 1, mean)
  y.bins.centres <- apply(y.bins, 1, mean)
  z.bins.centres <- apply(z.bins, 1, mean)
   
  if(length(channels)>1) {
    I <- ff::ff(0, dim = c(image.size, length(channels)))
  } else {
    I <- ff::ff(0, dim = image.size) 
  }
  
  if(method == 'photon') {
    if (!("phot" %in% colnames(points) || "photon" %in% colnames(points))) {
      stop("ERROR: Column 'phot[on]' must be present in input data frame to use method 'photon'.")
    }
    if (!("locprec" %in% colnames(points))) {
      stop("ERROR: Column 'locprec' must be present in input data frame to use method 'photon'.")
    }
    if (!("locprecz" %in% colnames(points))) {
      warning("WARNING: Column 'locprecz' not present in input data frame. Using values from locprec column for z.")
      points$locprecz <- points$locprec
    }
    cluster <- parallel::makeCluster(ncpu)
    doParallel::registerDoParallel(cluster)
    on.exit(parallel::stopCluster(cluster))
    # Calculate photon contribution of each point to each voxel
    # A point contributes to voxels within 3*locprec of its position
    i <- NULL
    tmp <- foreach::foreach(i = 1:nrow(points), .combine = 'c', .packages = c('pracma','ff')) %dopar% {
      # for(i in 1:nrow(points)) {
      p <- points[i,]
      # Find corners of the bounding box (i.e. +/- 3*locprec along each axis)
      range.vx <- findInterval(c(p$x - 3*p$locprec, p$x + 3*p$locprec), x.bins.breaks, all.inside = TRUE)
      range.vy <- findInterval(c(p$y - 3*p$locprec, p$y + 3*p$locprec), y.bins.breaks, all.inside = TRUE)
      range.vz <- findInterval(c(p$z - 3*p$locprecz, p$z + 3*p$locprecz), z.bins.breaks, all.inside = TRUE)
      # Make sure we don't get out of the image
      range.vx[2] <- min(range.vx[2], image.size[1])
      range.vy[2] <- min(range.vy[2], image.size[2])
      range.vz[2] <- min(range.vz[2], image.size[3])
      # Get all voxels in the bounding box
      voxels <- expand.grid(x = range.vx[1]:range.vx[2], y = range.vy[1]:range.vy[2], z = range.vz[1]:range.vz[2])
      for(j in 1:nrow(voxels)) {
        v <- voxels[j,]
        # Get voxel position in the original coordinate system
        vx <- x.bins.centres[v$x]
        vy <- y.bins.centres[v$y]
        vz <- z.bins.centres[v$z]
        # Test if voxel is within ellipsoid
        test <- (p$x-vx)^2/(9*p$locprec^2) + (p$y-vy)^2/(9*p$locprec^2) + (p$z-vz)^2/(9*p$locprecz^2)
        if(test>1) { next }
        sigma <- p$locprec * sqrt(2)
        sigma.z <- p$locprecz * sqrt(2)
        if(length(channels) > 1) {
          c <- match(p$channel, channels)
          if(!is.na(c)) {
            I[v$x,v$y,v$z,c] <- I[v$x,v$y,v$z,c] + p$phot/4 * (pracma::erf((vx-p$x+0.5)/(sigma)) - pracma::erf((vx-p$x-0.5)/(sigma))) * (pracma::erf((vy-p$y+0.5)/(sigma)) - pracma::erf((vy-p$y-0.5)/(sigma))) * (pracma::erf((vz-p$z+0.5)/(sigma.z)) - pracma::erf((vz-p$z-0.5)/(sigma.z)))
          }
        } else {
          I[v$x,v$y,v$z] <- I[v$x,v$y,v$z] + p$phot/4 * (pracma::erf((vx-p$x+0.5)/(sigma)) - pracma::erf((vx-p$x-0.5)/(sigma))) * (pracma::erf((vy-p$y+0.5)/(sigma)) - pracma::erf((vy-p$y-0.5)/(sigma))) * (pracma::erf((vz-p$z+0.5)/(sigma.z)) - pracma::erf((vz-p$z-0.5)/(sigma.z)))
        }
        NULL # Return NULL to save memory
      }
    }
  } else if(method == 'histogram') {
    if(length(channels) > 1) {
      for(c in channels) {
        pts <- points[which(points$channel == c), c('x','y','z')]
        bins <- binning(x = pts, y = NULL, nbins = image.size[1:3])
        I[,,,c] <- bins$table.freq
      }
    } else {
      if(length(channels == 1)) {
        pts <- points[which(points$channel == channels), c('x','y','z')]
        bins <- binning(x = pts, y = NULL, nbins = image.size)
      } else {
        bins <- binning(x = points[,c('x','y','z')], y = NULL, nbins = image.size)
      }
      I <- array(bins$table.freq, dim = dim(bins$table.freq))
    }
  } else {
    stop(paste0("ERROR: Unknown method: ", method))
  }
  return(I[])
}

#' points_from_roi
#'
#' Extract points within given bounding box.
#' Points are translated so that (0,0,0) correspond to the bounding box corner defined by
#' roi['min',c('x','y','z')]
#'
#' @param points a point set as a data frame of coordinates with columns x,y,z.
#' @param roi a data frame with columns x,y,z and rows min and max defining a bounding box
#' @return a data frame with same columns as input
#' @export

points_from_roi <- function(points, roi) {
  idx.to.keep <- which(points$x>roi['min','x'] & points$x<roi['max','x'] 
                       & points$y>roi['min','y'] & points$y<roi['max','y']
                       & points$z>roi['min','z'] & points$z<roi['max','z'])
  if(length(idx.to.keep)>0) {
      return(points[idx.to.keep,, drop = FALSE])
  } else {
      stop("ERROR: ROI is empty.")
  }
}

#' locs2ps
#'
#' Cluster localizations into point sets using DBSCAN
#'
#' @param points a point set as a data frame of coordinates with columns x,y,z.
#' @param eps DBSCAN parameter, size of the epsilon neighbourhood
#' @param minPts DBSCAN parameter, number of minimum points in the eps region
#' @param keep.locprec logical (default: TRUE), whether to preserve the localization precision columns
#' @param keep.channel logical (default: TRUE), whether to preserve channel information column
#' @param cluster.2d logical (default: FALSE), whether to cluster only using x,y (and ignore z)
#' @return a list of matrices with columns x,y,z and eventually locprec[z] and names set to the cluster indices.
#' @export

locs2ps <- function(points, eps, minPts, keep.locprec = TRUE, keep.channel = TRUE, cluster.2d = FALSE) {
  
  is.2d <- FALSE
  if(!("z" %in% colnames(points))) {
    is.2d <- TRUE
  }
  if(keep.locprec && keep.channel) {
    if(is.2d) {
      points <- points[, c("x", "y", colnames(points)[grep("locprec|channel", colnames(points))])]
    } else {
      points <- points[, c("x", "y", "z", colnames(points)[grep("locprec|channel", colnames(points))])]
    }
  } else if(keep.locprec && !keep.channel) {
    if(is.2d) {
      points <- points[, c("x", "y", colnames(points)[grep("locprec", colnames(points))])]
    } else {
      points <- points[, c("x", "y", "z", colnames(points)[grep("locprec", colnames(points))])]
    }
  } else if(!keep.locprec && keep.channel) {
    if(is.2d) {
      points <- points[, c("x", "y", colnames(points)[grep("channel", colnames(points))])]
    } else {
      points <- points[, c("x", "y", "z", colnames(points)[grep("channel", colnames(points))])]
    }
  } else {
    if(is.2d) {
      points <- points[, c("x", "y")]
    } else {
      points <- points[, c("x", "y", "z")]
    }
  }
  if(cluster.2d || is.2d) {
    dbr <- dbscan::dbscan(points[, c("x", "y")], eps = eps, minPts = minPts, borderPoints = FALSE)
  } else {
    dbr <- dbscan::dbscan(points[, c("x", "y", "z")], eps = eps, minPts = minPts, borderPoints = FALSE)
  }
  points$site <- dbr$cluster
  idx.to.remove <- which(points$site==0, arr.ind = TRUE)
  if(length(idx.to.remove)>0) {
    points <- points[-idx.to.remove,, drop = FALSE]
  }
  PS <- split(points, points$site)
  PS <- lapply(PS, function(x) { x["site"] <- NULL; as.matrix(x) })
  return(PS)
}

#' idx2rowcol
#'
#' Convert indices into a dist object to row, column coordinates of the corresponding distance matrix
#'
#' @param idx vector of indices
#' @param n size of the n x n distance matrix
#' @return a matrix with two columns nr and nc
#' @export

idx2rowcol <- function(idx,n) {
  nr <- ceiling(n-(1+sqrt(1+4*(n^2-n-2*idx)))/2)
  nc <- n-(2*n-nr+1)*nr/2+idx+nr
  return(cbind(nr,nc))
}

#' rotz
#'
#' Create a rotation matrix representing a rotation of theta radians about the z-axis
#'
#' @param theta angle in radians
#' @return a 3x3 rotation matrix
#' @export

rotz <- function(theta) {
  c <- cos(theta)
  s <- sin(theta)
  R <- matrix(c(c, -s, 0, s, c, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  return(R)
}

#' roty
#'
#' Create a rotation matrix representing a rotation of theta radians about the y-axis
#'
#' @param theta angle in radians
#' @return a 3x3 rotation matrix
#' @export

roty <- function(theta) {
  c <- cos(theta)
  s <- sin(theta)
  R <- matrix(c(c, 0, s, 0, 1, 0, -s, 0, c), nrow = 3, ncol = 3, byrow = TRUE)
  return(R)
}

#' rotx
#'
#' Create a rotation matrix representing a rotation of theta radians about the x-axis
#'
#' @param theta angle in radians
#' @return a 3x3 rotation matrix
#' @export

rotx <- function(theta) {
  c <- cos(theta)
  s <- sin(theta)
  R <- matrix(c(1, 0, 0, 0, c, -s, 0, s, c), nrow = 3, ncol = 3, byrow = TRUE)
  return(R)
}

#' local_densities
#' 
#' Compute local point density at each point of a point set
#' 
#' Local density is computed as in Ning X, Li F, Tian G, Wang Y (2018) 
#' An efficient outlier removal method for scattered point cloud data. 
#' PLOS ONE 13(8):e0201280. https://doi.org/10.1371/journal.pone.0201280
#'
#' @param X point set, a N x D matrix
#' @param k (optional) number of nearest neighbors used (defaults to all points).
#' @return vector of density value for each point
#' @export
#' 

local_densities <- function(X, k = NULL) {
  
  n <- nrow(X)
  if(is.null(k)) {
    k <- n-1
  }

  knnDist <- dbscan::kNNdist(X, k, all = TRUE)
  d <- apply(knnDist, 1, mean)
  densities <- rep(0, n)
  for(i in 1:n) {
    for(j in 1:k) {
      densities[i] <- densities[i] + exp(-knnDist[i,j]/d[i])
    }
  }
  densities <- densities/k
  
  return(densities)
}

#' find_elbow
#' 
#' Find elbow in a 2D curve represented by a list of ordered values
#' 
#' This function finds the point with maximum distance from the line between
#' the first and last points.
#' Adapted from StackOverflow:
#' http://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
#' 
#' @param values vector of values in decreasing order
#' @return index and value of the selected point
#' @export

find_elbow <- function(values) {
  
  n <- length(values)
  
  # Check to see if values are ordered
  if (is.unsorted(values) && is.unsorted(rev(values))) {
    stop("Values must be in decreasing order") 
  }
  
  coords <- cbind(seq.int(n), as.matrix(values))
  # Vector between first and last point
  line <- coords[n,] - coords[1,]
  # Normalize the vector
  line <- line/sqrt(sum(line^2));
  # Vector between all points and first point
  V1 <- sweep(coords, 2, coords[1,]) # default function in sweep is - (minus)
  # To calculate the distance to the line, we split V1 into two
  # components, V2 that is parallel to the line and V3 that is perpendicular.
  # The distance is the norm of the part that is perpendicular to the line.
  # We find V2 by projecting V1 onto the line by taking the scalar product
  # of V1 with the unit vector on the line (this gives us the length of the
  # projection of V1 onto the line).
  # We get V2 by multiplying this scalar product by the unit vector.
  # The perpendicular vector is V1 - V2
  V2 <- (V1 %*% line) %*% line
  V3 <- V1 - V2
  # Distance to line is the norm of V3
  dist <- sqrt(rowSums(V3^2))
  idx <- which.max(dist)
  
  return(coords[idx,])
}

#' denoise 
#'
#' Point density is estimated using a Gaussian mixture model and points in low
#' density regions are considered as noise and removed. 
#'
#' @param points a data frame with columns x,y,z.
#' @param k integer, number of mixture components for the GMM
#' @param prob probability level in the range [0,1] to identify high density regions
#' @return a point set
#' @export

denoise <- function(points, k = 16, prob = 0.3) {
  model <- mclust::densityMclust(points[, c("x","y","z")], G = k, plot = FALSE)
  # Find and keep high-density regions
  hdr.threshold <- mclust::hdrlevels(model$density, prob = prob)
  noise.idx <- which(model$density<hdr.threshold)
  if(length(noise.idx)>0) {
    return(points[-noise.idx,])
  } else {
    return(points)
  }
}

#' group_events 
#'
#' Localisation events are grouped by recursively clustering mutual nearest neighbours.
#' Neighbours are determined using the Mahalanobis distance to account for anisotropy in 
#' the localisation precision. Since the Mahalanobis distance has approximately a 
#' chi-squared distribution, a distance threshold can be chosen from a chi-squared table 
#' where the number of degrees of freedom is the dimension and alpha can be seen
#' as the probability of missing a localization event generated from the same fluorophore
#' as the event under consideration.
#'
#' @param points a data frame with columns x,y,z.
#' @param locprec localization precision in x,y
#' @param locprecz localization precision along z, defaults to locprec
#' @param p confidence level, see description. Defaults to 0.1
#' @return a list with two elements:
#'   \itemize{
#'     \item points: a point set as data frame with columns x,y,z
#'     \item membership: a vector of integers indicating the cluster to which each input point is allocated.
#'   }
#' @export

group_events <- function(points, locprec = NULL, locprecz = NULL, p = 0.1) {
  if(is.null(locprec)) {
    stop("Value for the locprec parameter required")
  }
  if(is.null(locprecz)) {
    locprecz <- locprec
  }
  points <- points[, c("x", "y", "z"), drop = FALSE]
  if(nrow(points) == 1) {
    pts <- points
    membership <- 1
  } else {
    # Transform coordinates such that the Euclidean distance between the transformed 
    # points is equal to the Mahalanobis distance with locprec^2 variance in the original
    # space. This avoids computing a potentially very large custom distance matrix for 
    # the nearest neighbours search which can be efficiently done by default with the
    # Euclidean distance.
    sigma <- c(locprec, locprec, locprecz)
    pts <- t(apply(points, 1, function(x) {x/sigma}))
    threshold <- stats::qchisq(p = p, df = ncol(points), lower.tail = FALSE)
    pts <- cbind(id = seq(1:nrow(pts)), pts)
    membership <- seq(1:nrow(pts))
    nn <- RANN::nn2(pts[,-1], k = 2, searchtype = 'radius', radius = threshold)
    idx <- nn$nn.idx[,-1]
    iter <- 0
    while(sum(idx) > 0) {
      iter <- iter + 1
      for(i in 1:nrow(pts)) {
        if(idx[i] > 0 && idx[idx[i]] == i) { # Identify mutual nearest neighbours
          to.update <- c(pts[i,1], pts[idx[i],1])
          if(iter>1) {
            # Both point i and its nearest neighbour can now be the centre of mass of multiple points
            to.update <- which(membership == membership[pts[i,1]])
            to.update <- c(to.update, which(membership == membership[pts[idx[i],1]]))
          }
          membership[to.update] <- pts[i,1]
          pts[i, -1] <- colMeans(points[to.update,,drop = FALSE])/sigma
          pts[idx[i],] <- NA
          idx[i] <- 0 # indicate we've already processed this point
        }
      }
      pts <- pts[stats::complete.cases(pts),, drop = FALSE]
      if(nrow(pts) > 1) {
        nn <- RANN::nn2(pts[,-1], k = 2, searchtype = 'radius', radius = threshold)
        idx <- nn$nn.idx[,-1]
      } else {
        idx <- 0
      }
    }
    if(nrow(pts) == 1) {
      pts <- pts[,-1, drop = FALSE] * sigma
    } else {
      pts <- t(apply(pts[,-1, drop = FALSE], 1, function(x) {x*sigma}))
    }
  }
  return(list(points = pts, membership = membership))
}

#' dist_to_line
#'
#' Compute distance between a set of points and a line defined by two points
#'
#' @param pts a data frame or matrix with 3 columns of coordinates
#' @param a vector of coordinates of a point on the line
#' @param b a second point on the line 
#' @return vector of distances
#' @export

dist_to_line <- function(pts, a = NULL, b = NULL) {
  if(is.null(a)) {
    stop("Need at least one point on the line.\n")
  }
  if(is.null(b)) {
    stop("Need a second point on the line.\n")
  }
  ba <- b - a
  l <- sum(ba^2)
  pa <- apply(pts, 1, function(p) { p - a})
  # The cross product gives a vector perpendicular to the line
  q <- apply(pa, 2, function(x) { pracma::cross(x, ba) })
  # The norm of the perpendicular vector is the distance to the line
  if(length(q)>0) {
    dist.to.line <- apply(q, 2, function(y) {sum(y^2)/l})
  } else {
    dist.to.line <- 0
  }
  return(sqrt(dist.to.line))
}

#' coloc_index
#'
#' Compute a co-localization index between two sets of points. 
#' Adapted from:
#' Willems and MacGillavry, A coordinate-based co-localization index to quantify
#' and visualize spatial associations in single-molecule localization microscopy.
#' Sci Rep 12, 4676 (2022). https://doi.org/10.1038/s41598-022-08746-4
#'
#' This can be seen as measuring the similarity between two spatial distributions.
#' Co-clustering in dense structures can give values above 1.
#' 
#' Localization precision is optional but if used then all locprec parameters
#' must be specified.
#' 
#' @param P1 a point set as matrix or data frame with columns x,y,z.
#' @param locprec1 (optional) localization precision in x,y for P1
#' @param locprecz1 (optional) localization precision along z for P1
#' @param P2 a point set as matrix or data frame with columns x,y,z.
#' @param locprec2 (optional) localization precision in x,y for P2
#' @param locprecz2 (optional) localization precision along z for P2
#' @return a list with two elements:
#'   \itemize{
#'     \item vector of co-localization indices for points in P1 relative to P2
#'     \item vector of co-localization indices for points in P2 relative to P1
#'   }
#' @export

coloc_index <- function(P1, locprec1 = NULL, locprecz1 = NULL, P2, locprec2 = NULL, locprecz2 = NULL) {
  
  if(ncol(P1) != ncol(P2)) {
    stop("The two point sets must have the same dimensions (i.e. number of columns)")
  } 
  if(!all(c("x", "y", "z") %in% colnames(P1))) {
    stop("P1 must have columns named x, y and z")
  }
  if(!all(c("x", "y", "z") %in% colnames(P2))) {
    stop("P2 must have columns named x, y and z")
  }
  if(!is.null(locprec1) && !is.null(locprecz1) && !is.null(locprec2) && !is.null(locprecz2)) {
    # To account for anisotropy in z, transform coordinates such that the Euclidean 
    # distance between the transformed points is equal to the Mahalanobis distance
    # with locprec^2 variance in the original space.
    sigma1 <- c(locprec1, locprec1, locprecz1)
    pts1 <- t(apply(P1[, c("x", "y", "z")], 1, function(x) {x/sigma1}))
    sigma2 <- c(locprec2, locprec2, locprecz2)
    pts2 <- t(apply(P2[, c("x", "y", "z")], 1, function(x) {x/sigma2}))
  } else {
    pts1 <- P1[, c("x", "y", "z")]
    pts2 <- P2[, c("x", "y", "z")]
  }
  nnP1 <- RANN::nn2(pts1, k = 2)
  mnnd1 <- mean(nnP1$nn.dists[,2])
  nnP2 <- RANN::nn2(pts2, k = 2)
  mnnd2 <- mean(nnP2$nn.dists[,2])
  # Assume that no point is within mnnd of more than 25% of the other points
  nn1 <- RANN::nn2(pts1, k = floor(nrow(pts1)/4), searchtype = 'radius', radius = mnnd1)
  nn2 <- RANN::nn2(pts1, k = floor(nrow(pts2)/4), searchtype = 'radius', radius = mnnd2)
  # Local density as number of points within mnnd
  ld1 <- apply(nn1$nn.idx, 1, function(x) {length(which(x>0))})
  ld2 <- apply(nn2$nn.idx, 1, function(x) {length(which(x>0))})
  
  # P2 points that are within mnnd2 of P1 points
  nn12 <- RANN::nn2(pts2, pts1, k = floor(nrow(pts2)/4), searchtype = 'radius', radius = mnnd2)
  ld12 <- apply(nn12$nn.idx, 1, function(x) {length(which(x>0))})
  coloc.idx12 <- ld12/mean(ld2)
  # P1 points that are within mnnd1 of P2 points
  nn21 <- RANN::nn2(pts1, pts2, k = floor(nrow(pts1)/4), searchtype = 'radius', radius = mnnd1)
  ld21 <- apply(nn21$nn.idx, 1, function(x) {length(which(x>0))})
  coloc.idx21 <- ld21/mean(ld1)
  
  return(list(coloc.idx12, coloc.idx21))
}
