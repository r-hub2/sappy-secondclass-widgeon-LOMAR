#' locs_from_csv
#'
#' Reads and filters single molecule localization events from a csv file as typically output by the SMAP software.
#' The main columns of interest are the coordinates (x, y, z), point set membership (site) and localization 
#' precision (locprec and locprecz).
#'
#' @param file a csv file with columns x[nm], y[nm], z[nm] and optionally site[numbers], channel, locprec[nm] and locprecz[nm], other columns are ignored.
#' @param roi region of interest, keep points within the specified volume. Must be a data frame with columns x,y,z and rows min and max defining a bounding box.
#' @param channels vector of integers indicating which channel(s) of a multicolour experiment to get data from.
#' @param frame.filter vector of min and max values, filter out points from frames outside the specified range. 
#' @param llrel.filter vector of min and max values, filter out points on log-likelihood (for fitted data).
#' @param locprec.filter filter out points with locprec value greater than the specified number. Points with locprec == 0 are also removed.
#' @param locprecz.filter filter out points with locprecz value greater than the specified number. Points with locprecz == 0 are also removed.
#' @return a data frame with columns x,y,z, optionally site, locprec and locprecz.
#' @examples
#' data.file <- system.file("test_data", "simulated_NUP107_data.csv", package = "LOMAR",
#'  mustWork = TRUE)
#' locs <- locs_from_csv(file = data.file, locprec.filter = 20)
#' @export

locs_from_csv <- function(file = NULL, 
                          roi = NULL,
                          channels = NULL,
                          frame.filter = NULL,
                          llrel.filter = NULL,
                          locprec.filter = 0,
                          locprecz.filter = 0) {
  message("Reading file ", file, "\n")
  if(tools::file_ext(file) == "mat") { #Assume we have SMAP output file in Matlab hdf5 format
    points <- as.data.frame(rhdf5::h5read(file, name = "saveloc/loc"))
  } else {
    points <- data.table::fread(file, header = TRUE, stringsAsFactors = FALSE)
  }
  colnames(points) <- sub("^xnm$", "x", colnames(points))
  colnames(points) <- sub("^ynm$", "y", colnames(points))
  colnames(points) <- sub("^znm$", "z", colnames(points))
  colnames(points) <- sub("^sitenumbers$", "site", colnames(points))
  colnames(points) <- sub("^locprecnm$", "locprec", colnames(points))
  colnames(points) <- sub("^locprecznm$", "locprecz", colnames(points))
  if(!is.null(roi)) {
    points <- points_from_roi(points, roi)
  }
  if(!is.null(channels) && (!("channel" %in% colnames(points))) ) {
    stop(paste0("ERROR: Column channel is missing in file ",file,".\n"))
  }
  idx.to.remove <- c()
  if(!is.null(channels)) {
    idx.to.remove <- c(idx.to.remove, which(!(points$channel %in% channels), arr.ind = TRUE))
  }
  if(!is.null(llrel.filter)) {
    if(!"LLrel" %in% colnames(points)) {
      stop(paste0("ERROR: Column LLrel (indicating log-likelihood) is missing in file ",file,".\n"))
    }
    if(is.null(llrel.filter['max'])) {
      llrel.filter['max'] <- Inf
    }
    if(is.null(llrel.filter['min'])) {
      llrel.filter['min'] <- -Inf
    }
    idx.to.remove <- c(idx.to.remove, which(points$LLrel < llrel.filter['min'] | points$LLrel > llrel.filter['max'], arr.ind = TRUE))
  }
  if(!is.null(frame.filter)) {
    if(!"frame" %in% colnames(points)) {
      stop(paste0("ERROR: Column frame is missing in file ",file,".\n"))
    }
    if(is.null(frame.filter['max'])) {
      frame.filter['max'] <- Inf
    }
    if(is.null(frame.filter['min'])) {
      frame.filter['min'] <- -Inf
    }
    idx.to.remove <- c(idx.to.remove, which(points$frame <= frame.filter['min'] | points$frame >= frame.filter['max'], arr.ind = TRUE))
  }
  if(locprec.filter || locprecz.filter) {
    if(length(grep("locprec",colnames(points)))==0) {
      stop(paste0("ERROR: Can't filter on localization precision: no locprec or locprecz column in file ",file,".\n"))
    }
    ## Remove points with low localization precision (i.e. high locprec value)
    ## Also remove points with locprec of 0 as this is not expected
    if("locprec" %in% colnames(points) && locprec.filter) {
      idx.to.remove <- c(idx.to.remove, which(points$locprec > locprec.filter | points$locprec == 0, arr.ind = TRUE))
    }
    if("locprecz" %in% colnames(points) && locprecz.filter) {
      idx.to.remove <- c(idx.to.remove, which(points$locprecz > locprecz.filter | points$locprecz == 0, arr.ind = TRUE))
    }
  }
  if(length(idx.to.remove)>0) {
    points <- points[-idx.to.remove,, drop = FALSE]
  }
  return(points)
}

#'  point_sets_from_locs
#' 
#'  Extracts list of point sets from a data frame of single molecule localization coordinates.
#'  By default, uses point set membership indicated in the site column.
#'  
#' @param locs, a data frame with columns x[nm], y[nm], z[nm] and optionally site[numbers], locprec[nm] and locprecz[nm], other columns are ignored.
#' @param channels vector of integers indicating which channel(s) of a multicolour experiment to extract point sets from.
#' @param min.cardinality filter out point sets with less than the specified number of points.
#' @param max.cardinality filter out point sets with more than the specified number of points.
#' @param crop.size remove points from a set if they are further away than the specified distance from the center of the set.
#' @param keep.locprec logical (default:TRUE). Whether to keep locprec information for each point.
#' @param sample.size returns this number of randomly selected point sets. Selects the point sets after applying eventual filtering.
#' @param ignore.site logical (default: FALSE), set to TRUE if point set membership is not present or needed.
#' @param cluster.points logical (default: FALSE), whether to cluster the points using DBSCAN (only if ignore.site is also TRUE).
#' @param eps DBSCAN parameter, size of the epsilon neighbourhood
#' @param minPts DBSCAN parameter, number of minimum points in the eps region
#' @return a list of matrices with columns x,y,z, optionally locprec and name set to the value of the site column (if applicable).
#' @examples
#' data.file <- system.file("test_data", "simulated_NUP107_data.csv", package = "LOMAR",
#'  mustWork = TRUE)
#' locs <- locs_from_csv(file = data.file, locprec.filter = 20)
#' point.sets <- point_sets_from_locs(locs, keep.locprec = TRUE, min.cardinality = 15)
#' @export

point_sets_from_locs <- function(locs = NULL,
                                 channels = NULL,
                                 min.cardinality = NULL,
                                 max.cardinality = NULL,
                                 crop.size = NULL,
                                 keep.locprec = TRUE,
                                 sample.size = NULL,
                                 ignore.site = FALSE,
                                 cluster.points = FALSE,
                                 eps = NULL,
                                 minPts = NULL) {
  PS <- list()
  if(keep.locprec && length(grep("locprec", colnames(locs)))==0) {
    stop(paste0("ERROR: Can't keep localization precision: no locprec or locprecz column in file ",file,".\n"))
  }
  if(!is.null(channels) && (!("channel" %in% colnames(locs))) ) {
    stop(paste0("ERROR: Column channel is missing in file ",file,".\n"))
  }
  if(!is.null(channels)) {
    idx.to.remove <- which(!(locs$channel %in% channels), arr.ind = TRUE)
    if(length(idx.to.remove)>0) {
      locs <- locs[-idx.to.remove,, drop = FALSE]
    }
  }
  if(!ignore.site) { # Use site information provided
    # Site == 0 indicates noise/low quality event
    idx.to.remove <- which(locs$site==0, arr.ind = TRUE)
    if(length(idx.to.remove)>0) {
      locs <- locs[-idx.to.remove,, drop = FALSE]
    }
    if(keep.locprec) {
      locs <- locs[, c("x", "y", "z", colnames(locs)[grep("locprec", colnames(locs))], "site")]
    } else {
      locs <- locs[, c("x", "y", "z", "site")]
    }
    PS <- split(locs, locs$site)
    PS <- lapply(PS, function(x) { x["site"] <- NULL; as.matrix(x) })
  } else if(cluster.points) { # Ignore site info if provided and generate new one by clustering
    PS <- locs2ps(locs, eps = eps, minPts = minPts)
  } else { # Ignore site info if provided and don't cluster
    if(keep.locprec) { 
      PS[[1]] <- locs[, c("x", "y", "z", colnames(locs)[grep("locprec", colnames(locs))])]
    } else {
      PS[[1]] <- locs[, c("x", "y", "z")]
    }
  }
  if(!is.null(crop.size)) {
    PS <- lapply(PS, function(x) {crop_point_set(x, crop.size)})
  }
  if(!is.null(min.cardinality)) {
    idx.to.remove <- NULL
    for(j in 1:length(PS)) {
      if(nrow(PS[[j]])<min.cardinality) {
        idx.to.remove <- c(idx.to.remove, j)
      }
    }
    if(length(idx.to.remove)>0) {
      PS[idx.to.remove] <- NULL
    }
  }
  if(!is.null(max.cardinality)) {
    idx.to.remove <- NULL
    for(j in 1:length(PS)) {
      if(nrow(PS[[j]])>max.cardinality) {
        idx.to.remove <- c(idx.to.remove, j)
      }
    }
    if(length(idx.to.remove)>0) {
      PS[idx.to.remove] <- NULL
    }
  }
  if(!is.null(sample.size)) {
    if(sample.size>length(PS)) {
      warning("Selected sample size is greater than number of point sets. Using all sets.")
    } else {
      PS <- PS[sample(length(PS), sample.size)]
    }
  }
  PS <- lapply(PS, as.matrix)
  return(PS)    
}

#' point_sets_from_tiffs
#'
#' Read in single molecule localization events from a series of 3D images in TIFF files where each image file
#' represents a point set.
#'
#' @param image_dir path to a directory containing the TIFF files.
#' @param pattern regular expression, select images whose file path matches the given pattern.
#' @param image.size vector of length 3 containing the size of the images along each dimension, e.g. c(40,40,40).
#' @param sample.size if set, selects this number of images at random. A sample size larger than the available number of samples produces a warning and is ignored.
#' @param sample.first if TRUE, samples are selected before applying any eventual filtering. This is more efficient as it avoids reading all data files.
#' @param min.cardinality if set, filter out all point sets with less than the specified number of points.
#' @param max.cardinality if set, filter out all point sets with more than the specified number of points.
#' @param crop.size vector of length 3 containing the desired reduced size of the images along each dimension, e.g. c(30,30,30).
#' @examples
#' data.dir <- system.file("test_data/img", package = "LOMAR", mustWork = TRUE) 
#' point_sets <- point_sets_from_tiffs(image_dir = data.dir, pattern = "\\.tiff?$",
#'  image.size = c(64, 64, 4), min.cardinality = 10)
#' @return a list with two elements:
#'   \itemize{
#'     \item point.sets: a list of point sets as matrices with columns x,y,z and 
#'     \item file.names: a vector of paths to the TIFF files from which the point sets were extracted.
#'   }
#' @export

point_sets_from_tiffs <- function(image_dir = NULL, 
                                  pattern = NULL,
                                  image.size = NULL,
                                  sample.size = NULL, 
                                  sample.first = FALSE,
                                  min.cardinality = NULL, 
                                  max.cardinality = NULL, 
                                  crop.size = NULL) {
  tnsr <- NULL
  if(is.null(image.size)) {
    stop("ERROR: Image size required.")
  }
  ## Read all TIFF images in the given directory into a 4D array
  image.list <- dir(path = image_dir, pattern = '\\.tiff?$', recursive = TRUE)
  if(!is.null(pattern)) {
    image.list <- grep(pattern = pattern, image.list, value = TRUE)
  }
  if(!is.null(sample.size) && sample.first) {
    if(sample.size>length(image.list)) {
      warning("Sample size larger than number of images.")
    } else {
      random.sample.idx <- sample(1:length(image.list), sample.size, replace = FALSE)
      image.list <- image.list[random.sample.idx]
    }
  }
  image.paths <- c()
  tnsr <- array(as.numeric(NA), dim = c(image.size, length(image.list)))
  for (i in 1:length(image.list)) {
    path <- paste(image_dir,"/", image.list[i], sep="")
    image <- EBImage::readImage(path)
    tnsr[,,,i] <- as.array(image)
    image.paths <- c(image.paths, path)
  }
  n <- dim(tnsr)[4]
  if(!is.null(crop.size)) {
    if(any(crop.size>image.size)) {
      stop("ERROR: Crop size must be smaller than image size along all dimensions.")
    }
    offset <- (image.size - crop.size)/2
    min <- offset + 1;
    max <- image.size - offset
    tnsr <- tnsr[min[1]:max[1], min[2]:max[2], min[3]:max[3],]
  }
  PS <- list()
  bkg <- 0
  ## Form point sets
  for (j in 1:n) {
    PS[[j]] <- as.matrix(which((tnsr[,,,j]>bkg), arr.ind = TRUE))
  }
  if(!is.null(min.cardinality)) {
    idx.to.remove <- NULL
    for(j in 1:length(PS)) {
      # Use length() instead of nrow() because point sets reduced 
      # to one point are not data frames anymore
      if(length(PS[[j]])/3<min.cardinality) {
        idx.to.remove <- c(idx.to.remove, j)
      }
    }
    if(length(idx.to.remove)>0) {
      PS[idx.to.remove] <- NULL
      image.paths <- image.paths[-idx.to.remove]
    }
  }
  if(!is.null(max.cardinality)) {
    idx.to.remove <- NULL
    for(j in 1:length(PS)) {
      # Use length() instead of nrow() because point sets reduced 
      # to one point are not data frames anymore
      if(length(PS[[j]])/3>max.cardinality) {
        idx.to.remove <- c(idx.to.remove, j)
      }
    }
    if(length(idx.to.remove)>0) {
      PS[idx.to.remove] <- NULL
      image.paths <- image.paths[-idx.to.remove]
    }
  }
  if(!is.null(sample.size) && !sample.first) {
    if(sample.size>length(PS)) {
      warning("Sample size larger than number of point sets.")
    } else {
      random.sample.idx <- sample(1:length(PS), sample.size, replace = FALSE)
      PS <- PS[random.sample.idx]
      image.paths <- image.paths[random.sample.idx]
    }
  }
  
  # Rename columns as x,y,z
  PS <- lapply(seq(PS), function(i) { X <- as.matrix(PS[[i]]); colnames(X) <- c("x","y","z"); return(X)})
  return(list(point.sets = PS, file.names = image.paths))
}

#' img2ps
#'
#' Read an image into a point set. 
#' The points are formed by extracting the coordinates of voxel values strictly above the given cut-off (default 0).
#'
#' @param img either a 2d or 3d array or a path to a file containing a 2d or 3d image.
#' @param bkg Extract points for values strictly above this (default = 0).
#' @param crop.size vector (of length 2 or 3) containing the desired reduced size of the images along each dimension, e.g. c(30,30,30).
#' @return a point set as matrix with columns x,y[,z]
#' @examples
#' img.file <- system.file("test_data/img", "alien1_3d.tif", package = "LOMAR",
#'  mustWork = TRUE) 
#' point_set <- img2ps(img = img.file, bkg = 0)
#' @export

img2ps <- function(img = NULL, bkg = 0, crop.size = NULL) {
  if(methods::is(img,"character") && file.exists(img)) {
    img <- as.array(EBImage::readImage(img))
  }
  image.size <- dim(img)
  n <- length(dim(img))
  if(!is.null(crop.size)) {
    if(any(crop.size>image.size)) {
      stop("ERROR: Crop size must be smaller than image size along all dimensions.")
    }
    offset <- (image.size - crop.size)/2
    min <- offset + 1;
    max <- image.size - offset
    if(n == 2) {
      img <- img[min[1]:max[1], min[2]:max[2]]
    } else if(n==3) {
      img <- img[min[1]:max[1], min[2]:max[2], min[3]:max[3]]
    }
  }
  PS <- as.matrix(which(img>bkg, arr.ind = TRUE))
  if(n == 2) {
    colnames(PS) <- c("x","y")
  } else if(n == 3) {
    colnames(PS) <- c("x","y","z")
  } else {
    stop(paste0("Image doesn't have 2 or 3 dimensions."))
  }
  return(PS)
}
