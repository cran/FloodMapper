\name{FM.readDEM}
\alias{FM.readDEM}

\title{
	Load Digital Elevation Model (DEM)
}

\description{
	After a model run is initiated, this function should be first executed to load DEM raster.
}

\usage{
FM.readDEM(demfile, waterdepth = 0, waterdepthfile = "", rerun = FALSE)
}

\arguments{
	\item{demfile}{
		full file path of the DEM raster.
	}
	\item{waterdepth}{
		(OPTIONAL) a number represents initial surface water depth (unit: m) to be applied to all grid cells.
	}
	\item{waterdepthfile}{
		(OPTIONAL) full file path of the initial surface water depth raster (unit: m). If a raster file is specified, the model will use this raster to initialize the surface water depth.
	}
	\item{rerun}{
		OPTIONAL, assign TRUE if this is a rerun.
	}
}

\author{
Xander Wang <xiuquan.wang@gmail.com>
}

\examples{
## Load FloodMapper package
library(FloodMapper)

## Refer to Step 2 in the sample code of FM.start():
## Step 2. Load the 1000km DEM raster for Lafayette Parish, 
##         LA, USA from the sample dataset.
}
