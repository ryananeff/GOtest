\name{plotGseaEnrTable}
\alias{plotGseaEnrTable}
\title{Plot GSEA Enrichment Table}
\description{
  This function plots the GSEA enrichment table.
}
\usage{
plotGseaEnrTable(GseaTable, x, species=c('human','mouse'),
 alpha=1, simplify.func=shorten.MSigDB.terms, simplify.func.par=list(),
 type=c('RES','NES'), col.up='#F8766D', col.down='#00BFC4', ...)
}
\arguments{
 \item{GseaTable}{a data.frame, must be (usually a subset of) the GSEA output from \code{\link{msigdb.gsea}}.}
 \item{x}{a data.frame of query gene signatures used for generating \code{GseaTable}.}
 \item{species}{Specifying the species for the input gene signatures.}
 \item{alpha}{power to scale the weights: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 or larger (over-weighted)}
 \item{simplify.func}{function to convert MSigDB terms to a simplified form. Set simplify.func=I if the terms should be kept unchanged. Default \code{\link{shorten.MSigDB.terms}}.}
 \item{simplify.func.par}{a named list to pass parameters to \code{simplify.func}.}
 \item{type}{content type to be plotted, either \code{RES} (running enrichment score) or \code{NES} (normalized enrichment score).}
 \item{col.up}{color for positive score.}
 \item{col.down}{color for negative score.}
 \item{...}{additional arguments used for \code{\link{plot}}.}
}
\details{
See example usage below.
}
\seealso{\code{\link{msigdb.gsea}}, \code{\link{shorten.MSigDB.terms}}}
\examples{
\dontrun{

#load approved symbols for protein-coding genes by HGNC;
#these will be the gene universe for the present toy example analysis
universe=curated.genesets(c('HGNC_universe'))$genesets[[1]]
n=length(universe)
#
#Create a toy example for weighted enrichment using GSEA and logistic regression
toy=data.frame(Gene=universe,Phenotype='XYZ',Z=rnorm(n,0,1),stringsAsFactors=FALSE)

#GSEA
fit2=msigdb.gsea(x=toy,genesets=c('C5.BP','C5.CC','C5.MF'),query.population=universe, background='query', method='GSEA')
head(fit2)
#Plot top 10 enriched GSEA terms
plotGseaEnrTable(GseaTable=fit2[1:10,], x=toy)
}
}
