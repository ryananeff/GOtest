\name{shorten.MSigDB.terms}
\alias{shorten.MSigDB.terms}
\title{Shorten MSigDB gene set terms}
\description{
  Simplify/shorten MSigDB gene set terms to make them more space economic when plotting.
}
\usage{
shorten.MSigDB.terms(x, maxlen=32, Tolower=TRUE, unchange=FALSE, parms=list())
}
\arguments{
 \item{x}{a vector of text strings of MSigDB gene set terms.}
 \item{maxlen}{maximum lenth allowed.}
 \item{Tolower}{logic, whether to convert upper case letters to lower case.}
 \item{unchange}{logic, whether to return \code{x} unchanged.}
 \item{parms}{a named list specifying the paramter values of \code{maxlen}, \code{Tolower}, and \code{unchange}. Useful when calling inside a function.}
}
\details{
See function source code for details.
}
\seealso{\code{\link{plotGseaEnrTable}}}
\examples{
shorten.MSigDB.terms('GO_NEGATIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE', 
	parms = list(maxlen=40, Tolower=TRUE))
}
