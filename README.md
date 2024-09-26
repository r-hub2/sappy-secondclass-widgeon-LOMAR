# LOMAR

An R package to deal with point sets from single molecule localization microscopy.
This package provides four sets of functionalities:
  - data input: read SMLM data as point sets either from csv files or from TIFF images
  - registration: point sets registration using different algorithms
  - topological data analysis: compute similarity between point sets using persistent homology
  - alpha-shapes: compute alpha-shapes of point sets and derive 3d shape features
  
  [A preprint describing the package is available on bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.30.493957v1).

# Installation

The package is [available on CRAN](https://cran.r-project.org/package=LOMAR).

![Download counts](https://cranlogs.r-pkg.org/badges/grand-total/LOMAR) (from the RStudio mirror)

``` R
install.packages("LOMAR")
```

To install the development version of this package, run (from within R):

``` R
library(devtools)
install_git('https://git.embl.de/heriche/lomar')
```

This package depends on these other packages:
  * data.table
  * FNN
  * foreach
  * parallel
  * doParallel
  * proxy
  * reshape2
  * pracma
  * transport
  * RANN
  * ff
  * dbscan
  * EBImage (from Bioconductor)
  * tools
  * rhdf5
  * mclust
  * alphashape3d

  
