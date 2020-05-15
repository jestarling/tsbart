\name{makeCutpoints}
\alias{makeCutpoints}
\title{makeCutpoints: Create list of cutpoints for a design matrix X.}
\usage{
makeCutpoints(X, gridlen=10000)
}
\arguments{
\item{X}{A matrix; the design matrix for the tsBART model. (Can use the output of \code{\link{makeModelMatrix}}, which expands categorical predictors to columns of indicators when creating design matrix.) Defaults to NULL.}

\item{gridlen}{Length of grid for continuous predictors.  Defaults to 10,000.}
}
\value{
\item{cutpoints}{A list of vectors, where each is a cutpoint for each covariate (column of design matrix X).}
}
\description{
This function creates a list of cutpoints for use with tsbart().
}
\examples{
\dontrun{
## Create list of cutpoints.
cuts <- makeCutpoints(X, gridlen=1000)
}}
\author{
Jennifer E. Starling
}
\keyword{tsbart}

