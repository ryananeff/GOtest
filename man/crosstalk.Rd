\name{crosstalk}
\docType{package}
\alias{crosstalk}
\title{Perform crosstalk analysis of enriched pathways}
\description{
  Pathway Crosstalk
}
\usage{
crosstalk(x, col.signature, col.pathway='MSigDB', col.overlap.size='Overlap.Size',
	col.adj.p='P.adj', col.elements='Overlap.Items', sep=';',
	adj.p.cutoff=0.05, min.size=5, min.common=3, ...)
}
\arguments{
 \item{x}{a data.frame as generated from \code{msigdb.gsea}.}
 \item{col.signature}{id of the column containing the query signatures.}
 \item{col.pathway}{id of the column containing the pathway gene sets.}
 \item{col.overlap.size}{id of the column containing the overlap size.}
 \item{col.adj.p}{id of the column containing the adjusted P value.}
 \item{col.elements}{id of the column containing the overlap elements which are concatenated with character sepcified by \code{sep}.}
 \item{sep}{a character string to concatenate the elements.}
 \item{adj.p.cutoff}{cutoff of adjusted P value to select significant pathways.}
 \item{min.size}{minimum number of elements required for a pathway to be selected.}
 \item{min.common}{minimum number of common elements required between a crosstalk pair.}
 \item{...}{additional parameters. Not used currently.}
}
\details{
See Jia et al (BMC Syst Biol. 2011, 5:S12) for a description of the pathway crosstalk analysis.
}
\examples{
#example usage
}
