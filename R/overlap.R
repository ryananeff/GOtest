## to compute enrichment of overlaps between set1 and set2 using hypergeometic test
## Input:
## 		Set*:		    2 column data frame, 1st col elements, 2nd col set category groups.
##                      NA or "" should have been removed in a preprocessing step.
##		totalBalls:		background set size.

SetOverlapHypGeoTest <- function(Set1, Set2, totalBalls, ReportOverlapElements=TRUE)
{
	totalElements <- union(Set1[,1],Set2[,1])
	totalElements=sort(totalElements)
	CatList1=unique(Set1[,2])
	CatList2=unique(Set2[,2])
	#elements in the rows and categorie in the cols
	esMat1 <- table2sparse(Set1,totalElements, CatList1)
	esMat2 <- table2sparse(Set2,totalElements, CatList2)	
	## category size	
	CatSize1 <- colSums(esMat1)
	CatSize2 <- colSums(esMat2)

	## overlap matrix
	CountMat <- as.matrix(crossprod(esMat1,esMat2))
	
	## p-value matrix
	totalWhite <- matrix(CatSize1, nrow=length(CatSize1), ncol=length(CatSize2), byrow=FALSE)
	totalBlack <- totalBalls - totalWhite
	totalDrawn <- matrix(CatSize2, nrow=length(CatSize1), ncol=length(CatSize2), byrow=TRUE)
	PMat       <- phyper(CountMat-1, totalWhite, totalBlack, totalDrawn, lower.tail=FALSE)
	FoldMat    <- (CountMat/totalDrawn)/(totalWhite/totalBalls)
	
	set1Mat    <- matrix(CatList1, nrow=length(CatSize1), ncol=length(CatSize2), byrow=FALSE)
	set2Mat    <- matrix(CatList2, nrow=length(CatSize1), ncol=length(CatSize2), byrow=TRUE)

	tabnames <- c("Set1", "Set2", 
			  "Overlap size", "Sampling size", "Positive size", "Background size", 
			  "Fold enrichment", "P value"
			 )
	EnrichTable <- data.frame(Set1=rep(CatList1,length(CatSize2)), Set2=rep(CatList2,each=length(CatSize1)),
		as.vector(CountMat), as.vector(totalDrawn), as.vector(totalWhite), as.vector(totalBalls),
		as.vector(FoldMat), as.vector(PMat),
		stringsAsFactors=FALSE
	)
	colnames(EnrichTable)[3:8] <- c("Overlap size", "Sampling size", "Positive size", "Background size", "Fold enrichment", "P value")
	if(ReportOverlapElements)
	{
		overMat=apply(esMat2,2,function(x,esMat1){
			Y = x * esMat1
			iy=colSums(Y)>0
			s=rep('',ncol(Y))
			if(!any(iy)) return(s)
			s[iy]=apply(Y[,iy,drop=FALSE],2,function(y,g) paste0(g[y>0], collapse=";"), g=rownames(esMat1))
			s
		},esMat1=esMat1)
		#overMat=apply(esMat2,2,function(x,esMat1) apply(x * esMat1,2,function(y,g) paste0(g[y>0], collapse=";"), g=rownames(esMat1)),esMat1=esMat1) #much slower than the above function in most cases
		EnrichTable[,"Overlap elements"]=as.vector(overMat)
	}
	
	EnrichTable[is.na(EnrichTable)] <- 0
	EnrichTable <- EnrichTable[order(EnrichTable[,"P value"]),]
	return(EnrichTable)
}
#convert a two-column table to a design matrix
#the first column of x denotes row entry and the second column of x denotes col entry of the design matrix
table2design=function(x,rows,cols){
	nr=length(rows)
	nc=length(cols)
	Mat1 <- matrix(FALSE, nrow=nr, ncol=nc, dimnames = list(rows, cols))
	i=match(x[,1],rows)+(match(x[,2],cols)-1) * nr
	Mat1[i]=TRUE
	Mat1
}
table2sparse=function(x,rows,cols){
	sparseMatrix(i=match(x[,1],rows),j=match(x[,2],cols), x=1,dims=c(length(rows),length(cols)), dimnames=list(rows, cols))
}
