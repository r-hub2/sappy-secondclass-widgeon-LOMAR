#' get_shape
#' 
#' Get the the alpha-shape of a point set.
#' If not given, the function automatically determines alpha using a downsampled 
#' point set. 
#' As a consequence, alpha and therefore the computed shape can vary slightly
#' between runs.
#' 
#' @param points a data frame with columns x, y, z.
#' @param alpha (optional) positive number
#' @return an alpha-shape object of class ashape3d
#' @export

get_shape <- function(points, alpha = NULL) {
  
  if(is.null(alpha)) {
    # Find a good value for alpha
    k <- floor(sqrt(nrow(points)))
    P1 <- downsample(points[,c("x", "y", "z")], round(nrow(points)*0.3), k = k)
    P1 <- P1[, c("x", "y", "z")]
    S <- as.matrix(stats::dist(P1))
    d <- stats::density(S)
    alpha <- d$x[which.max(d$y)]
  }  
  # Some points configurations cause ashape3d() to crash.
  # Adding noise to the points helps but we can sometimes still
  # end up in a degenerate case so we may need to try a few times.
  as <- NULL
  attempts <- 0
  while(is.null(as) && attempts <= 3) {
    attempts <- attempts + 1
    try({
      as <- alphashape3d::ashape3d(jitter(as.matrix(points[, c("x", "y", "z")], amount = 1)), alpha = alpha, pert = FALSE) # pert = F avoids storing an extra copy of the point set
    })
  }
  return(as)
}

#' get_surface_area
#'
#' Compute the surface area of an alpha-shape by summing the surfaces of the
#' boundary triangles
#'
#' @param as an alpha-shape object of class ashape3d
#' @return a numeric value
#' @export

get_surface_area <- function(as) {
  idx <- which(as$triang[,9] == 2)
  S <- 0
  for(i in idx) {
    A <- as$x[as$triang[i,1],]
    B <- as$x[as$triang[i,2],]
    C <- as$x[as$triang[i,3],]
    S <- S + 0.5*abs((B-A) %*% (C-A))
  }
  return(S)
}

#' shape_features_3d
#'
#' Compute shape features of a 3D alpha-shape object
#' 
#' Features are:
#'   - major.axis, minor.axis and least.axis: Lengths of the axes of the fitted ellipsoid
#'   - elongation: from 0 (line) to 1 (globular)
#'   - flatness: from 0 (flat) to 1 (spherical)
#'   - max.feret.diameter: Maximum Feret diameter
#'   - max.inscribed.radius: Radius of the largest inscribed sphere
#'   - sphericity: from 0 (not spherical) to 1 (perfect sphere)
#'   - concavity: fraction of the convex hull volume not in the object
#'   - volume
#'   - area: area of the surface of the alpha-shape
#'   - centre.x, centre.y and centre.z: coordinates of the geometric median of the points in the alpha-shape
#'
#' @param as an alpha-shape object of class ashape3d
#' @return a named vector of numeric values or NULL if no non-singular vertices
#' @export

shape_features_3d <- function(as) {
  # Remove singular vertices
  idx <- as$vertex[which(as$vertex[,5] < 3),1]
  if(length(idx)==0) { 
    # No non-singular vertices
    return(NULL)
  } else {
  # Axes of fitted ellipsoid
  pca <- stats::prcomp(as$x[idx,])
  axis.lengths <- pca$sdev
  # Elongation ranges from 0 -> line to 1 -> globular
  elongation <- axis.lengths[2]/axis.lengths[1]
  # Flatness goes from 0 -> flat to 1 -> spherical
  flatness <- axis.lengths[3]/axis.lengths[1]
  # Feret diameter
  max.feret.d <- max(proxy::dist(as$x))
  # Radius of maximum inscribed sphere
  centre <- pracma::geo_median(as$x)$p
  insc.r <- min(proxy::dist(t(centre), as$x))
  volume <- alphashape3d::volume_ashape3d(as)
  area <- get_surface_area(as)
  # Sphericity goes from 0 -> not spherical to 1 -> perfect sphere
  sphericity <- (36*pi*(volume)^2)/(area^3)
  # Concavity with 0 -> convex
  as <- alphashape3d::ashape3d(as, alpha = Inf)
  cxh.volume <- alphashape3d::volume_ashape3d(as, indexAlpha = 2)
  concavity <- (cxh.volume - volume) / cxh.volume
  return(c(major.axis = axis.lengths[1], minor.axis = axis.lengths[2], 
           least.axis = axis.lengths[3], elongation = elongation, 
           flatness = flatness, max.feret.diameter = max.feret.d, 
           max.inscribed.radius = insc.r, sphericity = sphericity, 
           concavity = concavity, volume = volume, area = area,
           centre.x = centre[1], centre.y = centre[2], centre.z = centre[3]))
  }
}

#' dist_to_boundary
#'
#' Given a point set and an alpha-shape, get the distance of each point
#' to the closest boundary point of the alpha-shape.
#' Points inside the shape get negative values.
#'
#' @param points a data frame with x,y,z columns
#' @param shape an object of class ashape3d with a single alpha value
#' @return vector of distances (negative values indicate points inside the shape)
#' @export

dist_to_boundary <- function(points, shape) {
  P <- as.matrix(points[,c("x","y","z")])
  is.inside <- alphashape3d::inashape3d(shape, indexAlpha = 1, P)
  is.inside <- ifelse(is.inside, -1, 1)
  d <- NA
  if(length(P)>0 && nrow(P)>0) {
    # Boundary points
    # Column 5 of the shape vertex matrix contains 1 for interior points,
    # 2 for regular points and 3 for singular points.
    # Boundary points are all non-interior points. 
    # We ignore singular points.
    idx <- shape$vertex[which(shape$vertex[,5] == 2),1]
    b.points <- shape$x[idx,c("x", "y", "z")]
    # Get nearest boundary point for each point
    nn <- RANN::nn2(b.points, P, k = 1) 
    d <- nn$nn.dists
  }
  return(d*is.inside)
}

#' scale_alpha_shape
#'
#' Uniformly scale an alpha-shape.
#' Note that this computes the alpha-shape of the scaled point set 
#' associated with the input alpha-shape.
#'
#' @param as an alpha-shape object of class ashape3d
#' @param s scaling factor
#' @return an object of class ashape3d
#' @export

scale_alpha_shape <- function(as, s) {
  Xs <- standardize_coordinates(as$x)
  X <- Xs[['X']] * s
  X <- restore_coordinates(X, mu = Xs[['mu']], sigma = Xs[['sigma']])
  a <- as$alpha * s
  AS <- NULL
  attempts <- 1
  while(is.null(AS) && attempts <= 5) {
    attempts <- attempts + 1
    try({
      AS <- alphashape3d::ashape3d(jitter(X, amount = s), alpha = a, pert = FALSE)
    })
  }
  return(AS)
}
