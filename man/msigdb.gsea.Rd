\name{msigdb.gsea}
\alias{msigdb.gsea}
\alias{merck.gsea}
\title{Gene Set Enrichment Analysis Using MSigDB or Merck Collections}
\description{
  These functions perform gene set enrichment analyses based on MSigDB or Merck collections.
}
\usage{
msigdb.gsea(x, query.population=NULL, genesets=c('C2.CP','C5.BP','C5.CC','C5.MF'),
 background='common', name.x='Input', name.go='MSigDB',
 method=c('hypergeometric','GSEA','logitreg'), adj="BH", species=c('human','mouse'),
 ReportOverlap=TRUE, ncores=1, gsea.alpha=1,
 permutations=ifelse(method=='GSEA',1000,0), iseed=12345)
merck.gsea(x, query.population=NULL, genesets=c('GO CC', 'GO MF', 'GO BP','KEGG PATHWAYS', 
 'INGENUITY PATHWAYS', 'GENEGO PATHWAYS'),
 background='common', name.x='Input', name.go='Merck',
 method=c('hypergeometric','GSEA','logitreg'), adj="BH", species=c('human','mouse','rat'),
 ReportOverlap=TRUE, ncores=1, gsea.alpha=1,
 permutations=ifelse(method=='GSEA',1000,0), iseed=12345)
}
\arguments{
 \item{x}{a data.frame of query gene signatures. See \code{Details}.}
 \item{query.population}{a population of the genes from which the query gene signatures are collected.}
 \item{genesets}{MSigDB/Merck gene set collections to be tested.}
 \item{background}{an integer or a keyword. See explanation in function \code{GOtest}.}
 \item{name.x}{a character string specifing the name of the query.}
 \item{name.go}{a character string specifing the name of the function gene sets.}
 \item{method}{algorithm to performing the enrichment test. See \code{Details}.}
 \item{adj}{approach to correct for multiple tests. See function \code{p.adjust}}
 \item{species}{Specifying the species for the input gene signatures.}
 \item{ReportOverlap}{whether to output overlapping elements.}
 \item{ncores}{number of CPU cores to be used for parallel computing.}
 \item{gsea.alpha}{power to scale the weights: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 or larger (over-weighted)}
 \item{permutations}{number of permutations for computing P value significance and controlling FDR for \code{GSEA} and \code{logitreg} methods. See \code{\link{GOtest}} for details.}
 \item{iseed}{seed for random number generation in permutations.}
}
\details{
For \code{msigdb.gsea}, known annotated gene sets in MSigDB will be loaded by function \code{\link[msigdb]{msigdb.genesets}} and the enrichment analysis is computed using function \code{\link{GOtest}}.

Other details can be found in function \code{\link{GOtest}}.
}
\value{
See output of function \code{GOtest}.
}
\seealso{\code{\link{GOtest}}, \code{\link{plotOverlap}}, \code{\link{p.adjust}}, \code{\link{GSEA}}}
\examples{
\dontrun{

#load approved symbols for protein-coding genes by HGNC;
#these will be the gene universe for the present toy example analysis
universe=curated.genesets(c('HGNC_universe'))$genesets[[1]]
n=length(universe)

###
# Example case 1
#
#load 23 functional gene sets curated by the MacArthur's Lab
MacArthur=curated.genesets(c('MacArthur'))$genesets
#Since variable MacArthur is a list, it needs to be converted to a data.frame
MacArthur.df=data.frame(Gene = unlist(MacArthur), Category = rep(names(MacArthur), sapply(MacArthur, length)))

#Test whether the MacArthur curated gene sets are enriched in MSigDB gene ontology collections using hypergeometric test
fit1=msigdb.gsea(x=MacArthur.df,genesets=c('C5.BP','C5.CC','C5.MF'), query.population=universe, background='query', method='hypergeometric')
#As expected, there are many significant overlaps
head(fit1)

###
# Example case 2
#
#Create a toy example for weighted enrichment using GSEA and logistic regression
toy=data.frame(Gene=universe,Phenotype='XYZ',Z=rnorm(n,0,1),stringsAsFactors=FALSE)

#GSEA
fit2=msigdb.gsea(x=toy,genesets=c('C5.BP','C5.CC','C5.MF'),query.population=universe, background='query', method='GSEA', permutations=10)
head(fit2)
#Plot top 10 enriched GSEA terms
plotGseaEnrTable(GseaTable=fit2[1:10,], x=toy, simplify.func.par=list(maxlen=40))

#logistic regression
#This is very slow
fit3=msigdb.gsea(x=toy,genesets=c('C5.BP','C5.CC','C5.MF'),query.population=universe, background='query',method='logitreg')
head(fit3)
}
}
