\name{plotGseaEnrTable}
\alias{plotGseaEnrTable}
\title{Plot GSEA Enrichment Table}
\description{
  This function plots the GSEA enrichment table.
}
\usage{
plotGseaEnrTable(GseaTable, x, go, genesets, species=c('human','mouse'),
 alpha=1, simplify.func=shorten.MSigDB.terms, simplify.func.par=list(),
 type=c('RES','NES'), col.up='#F8766D', col.down='#00BFC4', ...)
}
\arguments{
 \item{GseaTable}{a data.frame, must be (usually a subset of) the GSEA output from \code{\link{msigdb.gsea}} or \code{\link{GOtest}}.}
 \item{x}{a data.frame of query gene signatures used for generating \code{GseaTable}.}
 \item{go}{a data.frame of gene annotations. See \code{Details}.}
 \item{genesets}{a vector of gene annotation set names. See \code{Details}.}
 \item{species}{Specifying the species for the input gene signatures. Only used when \code{GseaTable} is computed from \code{\link{msigdb.gsea}}.}
 \item{alpha}{power to scale the weights: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 or larger (over-weighted)}
 \item{simplify.func}{function to convert MSigDB terms to a simplified form. Default \code{\link{shorten.MSigDB.terms}}.}
 \item{simplify.func.par}{a named list to pass parameters to \code{simplify.func}.}
 \item{type}{content type to be plotted, either \code{RES} (running enrichment score) or \code{NES} (normalized enrichment score).}
 \item{col.up}{color for positive score.}
 \item{col.down}{color for negative score.}
 \item{...}{additional arguments used for \code{\link{plot}}.}
}
\details{
Depending on the original function for calculating \code{GseaTable}, either \code{go} or \code{genesets} must be specified as in function \code{\link{GOtest}} or \code{\link{msigdb.gsea}}. When both are provided, only \code{go} will be used.

See example usage below.
}
\seealso{\code{\link{msigdb.gsea}}, \code{\link{GOtest}}, \code{\link{shorten.MSigDB.terms}}}
\examples{
\dontrun{

#load approved symbols for protein-coding genes by HGNC;
#these will be the gene universe for the present toy example analysis
universe=curated.genesets(c('HGNC_universe'))$Gene
n=length(universe)
#
#Create a toy example for weighted enrichment using GSEA and logistic regression
toy=data.frame(Gene=universe,Phenotype='XYZ',Z=rnorm(n,0,1),stringsAsFactors=FALSE)

#GSEA against MSigDB gene ontology biological process (C5.BP). This may take a long time to run.
fit2=msigdb.gsea(x=toy, genesets=c('C5.BP'), query.population=universe,
	background='query', method='GSEA', permutations=100)
head(fit2)
#Plot the GSEA running enrichment score for the top 10 GSEA terms
plotGseaEnrTable(GseaTable=fit2[1:10,], x=toy, genesets=c('C5.BP'))
}
}
