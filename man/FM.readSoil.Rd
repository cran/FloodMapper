\name{FM.readSoil}
\alias{FM.readSoil}

\title{
	Load Rasters for Soil Texture, Soil Depth (Groudwater Level), and Soil Moisture Content
}

\description{
	After land cover raster is loaded, this function should be executed to load rasters for soil texture, surface soil depth (groundwater level), and initial soil moisture content.
}

\usage{
FM.readSoil(soilfile, soildepth, soildepthfile = "", soilmoistcontent, 
    soilmoistcontentfile = "", rerun = FALSE)
}

\arguments{
	\item{soilfile}{
		full file path of the soil texture raster.
	}
	\item{soildepth}{
		a number indicates the depth of unsaturated soil (i.e., aquifer or groundwater level) to be applied to all grid cells, unit: m.
	}
	\item{soildepthfile}{
		full file path of the soil depth raster. If a raster file is specified, the model will use this raster to initialize the depth of unsaturated soil.
	}
	\item{soilmoistcontent}{
		a number indicates the initial soil moisture content (range: 0.0 - 1.0) to be applied to all grid cells.
	}
	\item{soilmoistcontentfile}{
		full file path of the soil moisture content raster. If a raster file is specified, the model will use this raster to initialize the soil moisture content.
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

## Refer to Step 4 in the sample code of FM.start():
## Step 4. Load the 1000km rasters for soil texture and 
##         soil depth from the sample dataset, set the 
##         initial soil moisture content to 0.4785.
}

