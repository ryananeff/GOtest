\name{plot.gsea}
\docType{package}
\alias{plot.gsea}
\alias{print.gsea}
\title{Present Gene Set Enrichment Analyses}
\description{
  Print and plot GSEA output.
}
\usage{
\method{plot}{gsea}(x, xlab='Rank', ylab='Running Enrichment Score', col.line='green', 
	col.indicator='red', indicator.bar=c('score','fixed'), ...)
\method{print}{gsea}(x, ...)
}
\arguments{
 \item{x}{a \code{gsea} object as generated from \code{GSEA}.}
 \item{xlab}{a title for the x axis.}
 \item{ylab}{a title for the y axis.}
 \item{col.line}{color for running score line.}
 \item{col.indicator}{color for gene tag.}
 \item{indicator.bar}{denote whether the indicator bar should show the running score or a fixed height.}
 \item{...}{additional parameters for \code{\link{plot}}.}
}
\details{
See function \code{GSEA} for object \code{gsea}.
}
\seealso{\code{\link{GSEA}}}
\examples{
#example usage
#create a random data set
n=10000
x=data.frame(ID=paste0('Gene',1:n),strength=sort(rnorm(n,0,1)),stringsAsFactors=FALSE)
fit1=GSEA(x=x,go=paste0('Gene',sample(n,100)))
plot(fit1)
}
