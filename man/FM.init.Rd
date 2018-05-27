\name{FM.init}
\alias{FM.init}

\title{
	Model Initialization
}

\description{
	This function needs to be executed in order to initiate a FloodMapper model run. A unique run name should be specified to identify the model run.
}

\usage{
FM.init(runname, startdatetime, enddatetime, outputinterval,
    internaltimestep = 30, debug = TRUE, rerun = FALSE,
    wdbreakpoints = NA, workdir = "")
}

\arguments{
	\item{runname}{
		a unique name to identify the model run. A folder with this name will be created under your current work directory.
	}
	\item{startdatetime}{
		start datetime of the model run, should be formated as "YYYY-MM-DD hh:mm:ss". E.g., "2016-06-09 12:00:00".
	}
	\item{enddatetime}{
		end datetime of the model run, should be formated as "YYYY-MM-DD hh:mm:ss". E.g., "2016-06-09 13:00:00".
	}
	\item{outputinterval}{
		time step for the model outputs, unit: seconds.
	}
	\item{internaltimestep}{
		(OPTIONAL) internal integral time step, unit: seconds, default: 30 seconds.
	}
	\item{debug}{
		(OPTIONAL) if the debug mode is enabled, a log file (named log.txt) will be generated under the run folder.
	}
	\item{rerun}{
		(OPTIONAL) assign TRUE if this is a rerun.
	}
	\item{wdbreakpoints}{
		(OPTIONAL) a vector defines how to display the flood depth with different colors (unit: mm). This setting will be only applied to the outputted PNG/PDF files.
	}
	\item{workdir}{
		Full path to a user-specified work directory. If not specified, the temporary folder defined by tempdir() will be used.
	}
}

\author{
Xander Wang <xiuquan.wang@gmail.com>
}

\examples{
## Load FloodMapper package
library(FloodMapper)

## Refer to Step 1 in the sample code of FM.start():
## Step 1. Initiate a new model run with a name of "Lafayette_1km"
}
