\name{GSEA}
\alias{GSEA}
\alias{fGSEA}
\title{Gene Set Enrichment Analyses}
\description{
  This function performs gene set enrichment analysis.
}
\usage{
GSEA(x, go, alpha=1, permutations=1000, ncores=1, iseed=12345)
fGSEA(x, go, alpha=1, permutations=1000, ncores=1, iseed=12345)
}
\arguments{
 \item{x}{data.frame, with at least two columns. See \code{Details}.}
 \item{go}{vector, a set of genes.}
 \item{alpha}{power to scale the weights. See \code{Details}.}
 \item{permutations}{number of permutations to estimate P value of the enrichment score.}
 \item{ncores}{number of CPU cores to be used.}
 \item{iseed}{integer, used as a seed for random number generation.}
}
\details{
\code{x} is a data.frame with at least two columns: the first column is a vector of gene ids and the second column is a vector of phenotype association statistics, eg, log fold change, or t-test statistics, whose sign and absolute values denote the direction and strength of the association.
\code{alpha} is power to scale the weights: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 or larger (over-weighted).
}
\value{
\code{GSEA} produces a \code{gsea} object, which is a list with the elements:
 \item{ES}{enrichment score}
 \item{indicator}{a vector of integers indicating the location of the genes of interest}
 \item{RES}{a vector of running enrichment scores}
 \item{NES}{normalized enrichment score accounted for the gene set size}
 \item{Pvalue}{p value of the enrichment score estimated by permutations}
}
\seealso{\code{\link{plot.gsea}}, \code{\link{print.gsea}}, \code{\link{msigdb.gsea}}}
\examples{
#example usage
#create a random data set
set.seed(1234)
n=10000
x=data.frame(ID=paste0('Gene',1:n),strength=sort(rnorm(n,0,1)),stringsAsFactors=FALSE)
fit1=GSEA(x=x,go=paste0('Gene',sample(n,100)))
plot(fit1)
}
