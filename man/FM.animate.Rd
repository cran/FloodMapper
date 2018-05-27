\name{FM.animate}
\alias{FM.animate}

\title{
	Create An Animation for An Existing Model Run
}

\description{
	Once a model run is completed successfully, users can use this function to create an animation into a *.gif file.
}

\usage{
FM.animate(runname = "", workdir = "")
}

\arguments{
	\item{runname}{
		(OPTIONAL) the name of an existing model run. If it is empty, the current model run will be selected.
	}
	\item{workdir}{
		(OPTIONAL) full path to user's work directory. If it is empty, the work directory of current model run will be selected.
	}
}

\author{
Xander Wang <xiuquan.wang@gmail.com>
}

\examples{
## Load FloodMapper package
library(FloodMapper)

s_return <- NA

## Assume that the current model run with a name of "Lafayette_1km" 
## has been completed successfully, you can run the following
## command to create an animation in *.gif format.
\donttest{
    s_return <- FM.animate("Lafayette_1km", workdir = tempdir())
}

## Check if there are any errors
if (!is.na(s_return)) {
    if (s_return != "") {
        cat("Creating animation is failed because of the following error:\n")
        cat(s_return)
    } else {
        cat("Creating animation is successful!\n")
    }
}
}
