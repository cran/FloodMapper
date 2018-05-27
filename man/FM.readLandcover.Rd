\name{FM.readLandcover}
\alias{FM.readLandcover}

\title{
	Load Land Cover Raster
}

\description{
	After DEM raster is loaded, this function should be executed to load land cover raster.
}

\usage{
FM.readLandcover(lcfile, rerun = FALSE)
}

\arguments{
	\item{lcfile}{
		full file path of the land cover raster.
	}
	\item{rerun}{
		(OPTIONAL) assign TRUE if this is a rerun.
	}
}

\author{
Xander Wang <xiuquan.wang@gmail.com>
}

\examples{
## Load FloodMapper package
library(FloodMapper)

## Refer to Step 3 in the sample code of FM.start():
## Step 3. Load the 1000km land cover raster from the sample dataset
}
