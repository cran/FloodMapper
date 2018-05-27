\name{FM.readInflow}
\alias{FM.readInflow}

\title{
	Load Inflow Data (OPTIONAL)
}

\description{
	After precipitatoin data is loaded, users may use this function to load inflow data in case that the selected domain is expected to receive water from its outside. The inflow data can come from monitoring stations or other raster datasets.
}

\usage{
FM.readInflow(datatype = 0, txtfile, startdatetime, enddatetime, interval)
}

\arguments{
	\item{datatype}{
		0 = monitoring stations, 1 = raster data (NOT supported in this version)
	}
	\item{txtfile}{
		full file path of the text file required for inflow data. If the data comes from multiple monitoring stations, the format of this text file should be given as follows:\cr
		#: ------------------------------------------\cr
		#: ID   Lon         Lat         Datafile\cr
		#: 1    -92.0296    30.1800     in1.txt\cr
		#: 2    -92.0260    30.1779     in2.txt\cr
		#: 3    ...         ...         ...\cr
		#: ------------------------------------------\cr
		If the inflow data comes from rasters, the format of this text file should be given as follows:\cr
		#: ------------------------------------------\cr
		#: ID   Rasterfile\cr
		#: 1    in_2016-06-09_12-00-00.tif\cr
		#: 2    in_2016-06-09_13-00-00.tif\cr
		#: 3    ...\cr
		#: ------------------------------------------
	}
	\item{startdatetime}{
		start datetime of inflow records, should be formated as "YYYY-MM-DD hh:mm:ss". E.g., "2016-06-09 12:00:00".
	}
	\item{enddatetime}{
		end datetime of inflow records, should be formated as "YYYY-MM-DD hh:mm:ss". E.g., "2016-06-09 13:00:00".
	}
	\item{interval}{
		interval of inflow records, unit: seconds. E.g., 1 hr = 60 min * 60 sec = 3600 sec.
	}
}

\author{
Xander Wang <xiuquan.wang@gmail.com>
}

\examples{
## Load FloodMapper package
library(FloodMapper)

## Refer to Step 7 in the sample code of FM.start():
## Step 7. Load inflow data from one weather station
}

