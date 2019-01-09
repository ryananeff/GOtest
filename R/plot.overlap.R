plotOverlap=function(x, col.set1, col.set2, col.set1.size, col.set2.size, 
	col.overlap.size, col.FE, col.Pvalue, col.text, cex.text=1.0, cex.lab=1.0, cex.legend=1.0, main='',
	order.by.size=TRUE,max.log10.p=320,cluster.by=NULL,fill=c('count','FE','text'),
	colorLabels=FALSE,show.set.size=TRUE,set.order=NULL,...){
#cluster.by, one of NULL, 'P' or 'FE' by which the rows and columns are clustered. When not NULL, it will override order.by.set.size
	if(!is.null(fill)) fill=match.arg(fill)
	x=as.data.frame(x)
	for(i in c(col.set1.size, col.set2.size, col.overlap.size, col.FE, col.Pvalue)) x[,i]=as.numeric(x[,i])
	#
	pmat <- -log10(c3toMatrix(x,col.set1, col.set2, col.Pvalue,1))
	pmat[is.infinite(pmat)] <- 1.3*max(pmat[is.finite(pmat)])
	pmat[pmat > max.log10.p] <- max.log10.p
	FE <- c3toMatrix(x, col.set1, col.set2, col.FE, 1)
	if(!is.null(fill)){
		if(fill=='count'){
			textMatrix <- c3toMatrix(x,col.set1, col.set2, col.overlap.size,0)
		}else if(fill=='text'){
			textMatrix <- c3toMatrix(x,col.set1, col.set2, col.text,'')
		}else{
			textMatrix=apply(FE,2,function(x) as.numeric(sprintf("%.0f",as.numeric(x))))
			rownames(textMatrix)=rownames(FE)
		}
	}else{
		textMatrix=NULL
	}
	#
	modname1=unique(x[,c(col.set1)])
	modname2=unique(x[,c(col.set2)])
	mod1Size <- unique(x[,c(col.set1,col.set1.size)])
	if(nrow(mod1Size) != length(modname1)) stop('Some entries of set1 have multiple set1.size entries\n')
	mod1Size=mod1Size[match(modname1,mod1Size[,1]),2]
	mod2Size <- unique(x[,c(col.set2,col.set2.size)])
	if(nrow(mod2Size) != length(modname2)) stop('Some entries of set2 have multiple set2.size entries\n')
	mod2Size=mod2Size[match(modname2,mod2Size[,1]),2]
	#use the same orders
	if(!is.null(fill)) textMatrix=textMatrix[modname1,modname2]
	FE=FE[modname1,modname2]
	pmat=pmat[modname1,modname2]
	if(!is.null(set.order)){
		if(is.null(set.order$Set1) || is.null(set.order$Set2)) stop('set.order must have two entries named Set1 and Set2\n')
		if(!all(set.order$Set1 %in% modname1)) stop('set.order$Set1 contains invalid entry\n')
		if(!all(set.order$Set2 %in% modname2)) stop('set.order$Set2 contains invalid entry\n')
		mod1Size=mod1Size[match(set.order$Set1,modname1)]
		mod2Size=mod2Size[match(set.order$Set2,modname2)]
		modname1=set.order$Set1
		modname2=set.order$Set2
		textMatrix=textMatrix[modname1,modname2]
		FE=FE[modname1,modname2]
		pmat=pmat[modname1,modname2]
	}else{
		if(!is.null(cluster.by)){
			if(cluster.by == 'P'){
				mod1Order=hclust(dist(pmat))$order
				mod2Order=hclust(dist(t(pmat)))$order
			}else if(cluster.by == 'FE'){
				mod1Order=hclust(dist(FE))$order
				mod2Order=hclust(dist(t(FE)))$order
			}else{
				stop('cluster.by must be one of NULL, "P" or "FE"\n')
			}
			mod1Size <- mod1Size[mod1Order]
			mod2Size <- mod2Size[mod2Order]
			modname1 <- modname1[mod1Order]
			modname2 <- modname2[mod2Order]
			textMatrix <- textMatrix[mod1Order, mod2Order]
			pmat <- matrix(pmat[mod1Order, mod2Order],nrow=dim(pmat)[1],ncol=dim(pmat)[2])
		}else{
			if(order.by.size==TRUE){
				mod1Order <- order(mod1Size, decreasing=T)
				mod2Order <- order(mod2Size, decreasing=T)
				mod1Size <- mod1Size[mod1Order]
				mod2Size <- mod2Size[mod2Order]
				modname1 <- modname1[mod1Order]
				modname2 <- modname2[mod2Order]
				textMatrix <- textMatrix[mod1Order, mod2Order]
				pmat <- matrix(pmat[mod1Order, mod2Order],nrow=dim(pmat)[1],ncol=dim(pmat)[2])
			}
		}
	}
	if(all(pmat == pmat[1,1])) pmat[1,1] = pmat[1,1]-1
	if(show.set.size){
		xSymbols=paste(modname2, ":", mod2Size)
		ySymbols=paste(modname1, ":", mod1Size)
	}else{
		xSymbols=modname2
		ySymbols=modname1
	}
	if(colorLabels==TRUE) {
		xLabels=paste0(" ", modname2)
		yLabels=paste0(" ", modname1)
	} else {
		xLabels=xSymbols
		yLabels=ySymbols
	}
	labeledHeatmap2(Matrix=pmat, colorLabels=FALSE, 
		xLabels=xLabels,
		yLabels=yLabels,
		xSymbols=xSymbols,
		ySymbols=ySymbols,
		main=main, 
		textMatrix=textMatrix,
		colors=blueWhiteRed2(200)[100:200],
		cex.text=cex.text, cex.lab=cex.lab, cex.legend=cex.legend, setStdMargins=FALSE,...)
	return(invisible())
}
c3toMatrix=function(dat,row.col=1,col.col=2,val.col=3,default.value=0){
	#dat is a data.frame or matrix with each row containing a record in the 2D matrix;
	#row.col, col.col and val.col specify the row, column and value for each record of the 2D matrix
	M=matrix(default.value,nrow=length(unique(dat[,row.col])),ncol=length(unique(dat[,col.col])))
	rownames(M)=unique(dat[,row.col])
	colnames(M)=unique(dat[,col.col])
	for(i in 1L:nrow(dat)) M[dat[i,row.col],dat[i,col.col]]=dat[i,val.col]
	M
}
