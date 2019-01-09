\name{curated.genesets}
\alias{curated.genesets}
\title{Curated Gene Sets}
\description{
  This function returns curated gene set collections
}
\usage{curated.genesets(sets=c('MacArthur'), return.data.frame=TRUE)}
\arguments{
 \item{sets}{A vector of possible set collections. See \code{details}.}
 \item{return.data.frame}{Logic, return a data.frame (when TRUE) or a list (when FALSE).}
}
\details{
The collection of genes includes:

\code{HGNC_universe}: approved symbols for 18,991 protein-coding genes according to HGNC as of Feb 9, 2015, (https://github.com/macarthur-lab/gene_lists);

\code{MacArthur}: a list of 23 functional gene sets curated by the MacArthur's Lab (https://github.com/macarthur-lab/gene_lists).

}
\value{
Either a data.frame or a list with three components as retunred by \code{\link[msigdb]{read.gmt}}
}
\examples{
##not run###
#x=curated.genesets(c('MacArthur'))
}
