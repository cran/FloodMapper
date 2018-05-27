\name{FM.readDraininlets}
\alias{FM.readDraininlets}

\title{
	Load Raster for Street Drain Inlets (OPTIONAL)
}

\description{
	After soil rasters are loaded, users may use this function to load the raster for street drain inlets in order to consider the effects of urban stormwater systems on mitigating floods.
}

\usage{
FM.readDraininlets(draininletfile, rerun = FALSE)
}

\arguments{
	\item{draininletfile}{
		full file path of the drain inlet raster.
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

## Refer to Step 5 in the sample code of FM.start():
## Step 5. Load the 1000km raster for street drain inlet
}
