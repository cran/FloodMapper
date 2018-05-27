\name{FM.readPrecip}
\alias{FM.readPrecip}

\title{
	Load Precipitation Data
}

\description{
	After all rasters (DEM, land cover, and soil) are loaded, users may use this function to load precipitation data. Precipitation data can come from multiple weather stations or gridded remote sensing data.
}

\usage{
FM.readPrecip(datatype = 0, txtfile, startdatetime, enddatetime, 
    interval, idwpower = 2)
}

\arguments{
	\item{datatype}{
		0 = weather stations, 1 = raster data from NWP models, reanalysis datasets, or other remote sensing datasets.
	}
	\item{txtfile}{
		full file path of the text file required for precipitation data. If the data comes from multiple weather stations, the format of this text file should be given as follows:\cr
		#: ------------------------------------------\cr
		#: ID   Lon         Lat         Datafile\cr
		#: 1    -92.0296    30.1800     p1.txt\cr
		#: 2    -92.0260    30.1779     p2.txt\cr
		#: 3    ...         ...         ...\cr
		#: ------------------------------------------\cr
		If the precipitation data comes from rasters, the format of this text file should be given as follows:\cr
		#: ------------------------------------------\cr
		#: ID   Rasterfile\cr
		#: 1    p_2016-06-09_12-00-00.tif\cr
		#: 2    p_2016-06-09_13-00-00.tif\cr
		#: 3    ...\cr
		#: ------------------------------------------
	}
	\item{startdatetime}{
		start datetime of precipitation records, should be formated as "YYYY-MM-DD hh:mm:ss". E.g., "2016-06-09 12:00:00".
	}
	\item{enddatetime}{
		end datetime of precipitation records, should be formated as "YYYY-MM-DD hh:mm:ss". E.g., "2016-06-09 13:00:00".
	}
	\item{interval}{
		interval of precipitation records, unit: seconds. For example, 1 hr = 60 min * 60 sec = 3600 sec.
	}
	\item{idwpower}{
		(OPTIONAL) the power of IDW method, default: 2, suggested options: 1, 2, or 3.
	}
}

\author{
Xander Wang <xiuquan.wang@gmail.com>
}

\examples{
## Load FloodMapper package
library(FloodMapper)

## Refer to Step 6 in the sample code of FM.start():
## Step 6. Load precipitation data from one weather station
}

