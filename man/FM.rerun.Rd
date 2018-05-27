\name{FM.rerun}
\alias{FM.rerun}

\title{
	Rerun An Existing Model
}

\description{
	Once a model run is completed successfully, users may want to adjust some parameters and rerun it. If this is the case, users can use this function.
}

\usage{
FM.rerun(runname, startdatetime, enddatetime, outputinterval, 
    internaltimestep = 30, debug = TRUE, wdbreakpoints = NA,
    workdir = "", animation = FALSE, bgtype = 0,
    aerialraster = "", pdfoutput = FALSE)
}

\arguments{
	\item{runname}{
		refer to FM.init().
	}
	\item{startdatetime}{
		refer to FM.init().
	}
	\item{enddatetime}{
		refer to FM.init().
	}
	\item{outputinterval}{
		refer to FM.init().
	}
	\item{internaltimestep}{
		refer to FM.init().
	}
	\item{debug}{
		refer to FM.init().
	}
	\item{wdbreakpoints}{
		refer to FM.init().
	}
	\item{workdir}{
		refer to FM.init().
	}
	\item{animation}{
		refer to FM.start().
	}
	\item{bgtype}{
		refer to FM.start().
	}
	\item{aerialraster}{
		refer to FM.start().
	}
	\item{pdfoutput}{
		refer to FM.start().
	}
}

\author{
Xander Wang <xiuquan.wang@gmail.com>
}

\examples{
## Load FloodMapper package
library(FloodMapper)

s_return <- NA

## Assume that the model run with a name of "Lafayette_1km" has been completed successfully,
## now you want to rerun this model after adjust some parameter files (under the input folder),
## you can run the following command to rerun the model.
\donttest{
    s_return <- FM.rerun(runname = "Lafayette_1km",
        startdatetime = "2016-08-12 10:00:00", 
        enddatetime = "2016-08-12 11:00:00", outputinterval = 3600, 
        internaltimestep = 300, wdbreakpoints = seq(0, 3000, 150),
        workdir = tempdir(), animation = TRUE, bgtype = 0,
        aerialraster = "", pdfoutput = FALSE)
}

## Check if there are any errors
if (!is.na(s_return)) {
    if (s_return != "") {
        cat("Model rerun is failed because of the following error:\n")
        cat(s_return)
    } else {
        cat("Model rerun is successful!\n")
    }
}
}
