\name{plotOverlap}
\docType{package}
\alias{plotOverlap}
\title{Plot Overlap With Heatmap}
\description{
  This function plots the overlap with a heatmap.
}
\usage{
plotOverlap(x, col.set1, col.set2, col.set1.size, col.set2.size, 
	col.overlap.size, col.FE, col.Pvalue, col.text, cex.text=1.0,
	cex.lab=1.0, cex.legend=1.0, main='', order.by.size=TRUE,
	max.log10.p=320, cluster.by=NULL, fill=c('count','FE','text'),
	colorLabels=FALSE, show.set.size=TRUE, set.order=NULL, ...)
}
\arguments{
 \item{x}{a data.frame of overlap analysis result.}
 \item{col.set1}{column id of set 1 in \code{x}; used in rows of the heatmap.}
 \item{col.set2}{column id of set 2 in \code{x}; used in columns of the heatmap.}
 \item{col.set1.size}{column id of set 1 size.}
 \item{col.set2.size}{column id of set 2 size.}
 \item{col.overlap.size}{column id of overlap size.}
 \item{col.FE}{column id of fold enrichment.}
 \item{col.Pvalue}{column id of P value.}
 \item{col.text}{column id of text string used for \code{fill}.}
 \item{cex.text}{magnitude of text size.}
 \item{cex.lab}{magnitude of lable size.}
 \item{cex.legend}{magnitude of legend text size.}
 \item{main}{title of the plot.}
 \item{order.by.size}{whether to sort the sets by size.}
 \item{max.log10.p}{maximum absolute log10 P value. Larger values will be fixed to this value.}
 \item{cluster.by}{one of NULL, 'P' and 'FE'. See \code{Details}.}
 \item{fill}{fill the cells with count, fold enrichment (FE) or text. NULL to turn off filling.}
 \item{colorLabels}{whether to color the set labels if the sets are named by colors.}
 \item{show.set.size}{whether to show set size in the row and column labels.}
 \item{set.order}{NULL or a list with two named vectors \code{Set1} and \code{Set2} which specify the odered element IDs of two sets.}
 \item{...}{additional arguments used in function \code{labeledHeatmap2}.}
}
\details{
A modified version of function \code{labeledHeatmap} from package \code{WGCNA} will be called to plot overlap in a heatmap. \code{cluster.by} specifies how to cluster the sets (either by 'P' value or 'FE'), and, if not NULL, it overrides \code{order.by.size}.
}
\seealso{\code{\link{GOtest}}}
\examples{
##not run###
}
