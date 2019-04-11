#plot the running score of GSEA
plot.gsea=function(x,xlab='Rank',ylab='Running Enrichment Score',col.line='green',col.indicator='red',indicator.bar=c('score','fixed'),...){
	#x is a list returned from GSEA
	stopifnot(inherits(x,"gsea"))
	indicator.bar=match.arg(indicator.bar)
	n=length(x$RES)
	plot(1:n,x$RES,type='l',xlab=xlab,ylab=ylab,col=col.line,...)
	lines(c(0,n),c(0,0))
	if(indicator.bar=='score'){
		for(i in x$indicator) lines(c(i,i),c(0,x$RES[i]),col=col.indicator)
	}else{
		j=max(abs(x$RES))
		for(i in x$indicator) lines(c(i,i),c(-j/10,j/10),col=col.indicator)
	}
}
print.gsea=function(x,...){
	#obj is a list returned from GSEA
	stopifnot(inherits(x,"gsea"))
	cat('An gsea object\n')
	cat('ES = ',x$ES,'\n')
	cat('NES = ',x$NES,'\n')
	cat('P value = ',x$Pvalue,'\n')
}
no.na=function(x){
	x[!is.na(x)]
}
plotGseaEnrTable=function(GseaTable, x, go, genesets, species=c('human','mouse'), alpha=1, simplify.func=shorten.MSigDB.terms, simplify.func.par=list(), type=c('RES','NES'), col.up='#F8766D', col.down='#00BFC4', ...){
	species=match.arg(species)
	type=match.arg(type)
	if(ncol(GseaTable) < 11) stop('Please provide unmodified output from GOtest/msigdb.gsea\n')
	if(type=='NES'){
		v=data.frame(nes=GseaTable$NES,name=simplify.func(GseaTable[,2],parms=simplify.func.par),stringsAsFactors=FALSE)
		v$u = v$nes > 0
		v$anes=abs(v$nes)
		v=v[order(v$u,-v$anes),]
		d <- setNames(v$anes,v$name)
		parms=list(...)
		horiz=FALSE
		NLetters=sapply(v$name,nchar)+countCapitalLetters(v$name)*0.4
		if((!is.null(parms$horiz)) && parms$horiz==TRUE) horiz=TRUE
		if(horiz){
			par(mar=c(4,max(NLetters)*0.45,0.5,0.5))
			xlab='Absolute NES'
			ylab=''
		}else{
			par(mar=c(max(NLetters)*0.45,4,0.5,0.5))
			xlab=''
			ylab='Absolute NES'
		}
		barplot(setNames(v$anes,v$name),col=c(col.down,col.up)[v$u+1],ylab=ylab,xlab=xlab,las=2,...)
		return(invisible())
	}
	if(!missing(go)){
		if((is.data.frame(go) || is.matrix(go)) && ncol(go)==2){
			go=apply(go,2,as.character)
		}
		addcomma=FALSE
	}else{
		if(is.character(genesets)){
			go=msigdb.genesets(sets=genesets,species=species,return.data.frame=TRUE)
		}else{
			stop('Error. Either go or genesets must be provided.\n')
		}
		addcomma=TRUE
	}
	colnames(go)=c('genes','set')
	GseaTable$simplified_terms=simplify.func(GseaTable[,2],parms=simplify.func.par)
	GseaTable$NES=sprintf("%.1f",GseaTable$NES)
	GseaTable$Pvalue=sprintf("%.1E",GseaTable$Pvalue)
	GseaTable$P.adj=sprintf("%.1E",GseaTable$P.adj)
	x=x[x[,2] %in% unique(GseaTable[,3]),]
	if(length(unique(x[,2])) > 1){
		n11=as.vector(table(x[,2]))
		if(!all(n11==n11[1])) stop('Error: input enrichment table contains multiple phenotypes/groups which have different gene background sizes. Please provide results from a single phenotype/group analysis or ensure they have common gene background.\n')
	}else{
		n11=nrow(x)
	}
	x=x[order(x[,2],-x[,3]),]
	NLetters=sapply(GseaTable$simplified_terms,nchar)+countCapitalLetters(GseaTable$simplified_terms)*0.4
	par(mar=c(3,max(NLetters)*0.45,0,10))
	plot(c(1,n11[1]),c(1,nrow(GseaTable)+0.5),type='n',bty='n',xlab='',ylab='',yaxt='n',...)
	par(yaxt='s')
	Res1=vector(mode = "list", length = nrow(GseaTable))
	maxres=0
	for(i in 1:nrow(GseaTable)){
		s1=ifelse(addcomma,paste0(GseaTable[i,1],':',GseaTable[i,2]),GseaTable[i,2])
		go1=go[go[,'set'] == s1,'genes']
		x1=x[x[,2]==GseaTable[i,3],-2]
		x1=data.frame(strength=x1[,2],tag=x1[,1] %in% go1,stringsAsFactors=FALSE)
		Res1[[i]]=GSEAfunc(x1, alpha=alpha, isOrdered=TRUE)
		a=max(abs(Res1[[i]]$RES))
		maxres=ifelse(maxres>=a,maxres,a)
	}
	for(i in 1:nrow(GseaTable)){
		axis(side=2, at=nrow(GseaTable)-i+1, labels=GseaTable$simplified_terms[i],las=2,tick=FALSE)
		axis(side=4, at=nrow(GseaTable)-i+1, labels=paste(GseaTable[i,c('NES','Pvalue','P.adj')],collapse='   '),las=2,tick=FALSE)
		Res11=Res1[[i]]
		Res11$RES=Res11$RES/maxres*0.5
		for(j in Res11$indicator) lines(c(j,j),c(nrow(GseaTable)-i+1,nrow(GseaTable)-i+1+Res11$RES[j]),col=ifelse(Res11$RES[j]>0,col.up,col.down))
	}
}
countCapitalLetters=function(x){
	sapply(regmatches(x, gregexpr("[A-Z]", x, perl=TRUE)), length)
}
shorten.MSigDB.terms=function(x,maxlen=32,Tolower=TRUE,unchange=FALSE,parms=list()){
	if(length(parms)>0){
		if(!is.null(parms$Tolower)) Tolower=parms$Tolower
		if(!is.null(parms$maxlen)) maxlen=parms$maxlen
		if(!is.null(parms$unchange)) unchange=parms$unchange
	}
	if(unchange) return(x)
	if(Tolower) x=tolower(x)
	x=gsub('_',' ',x)
	x=sub('^go ','',x)
	x=sub('^kegg ','',x)
	x=sub('^reactome ','',x)
	x=sub('^biocarta ','',x)
	x=sub('^pid ','',x)
	x=sub('^naba ','',x)
	x=sub('gtpase','GTPase',x)
	x=sub('gpcr','GPCR',x)
	x=sub(' utr',' UTR',x)
	x=sub(' ecm',' ECM',x)
	x=sub('^ecm ','ECM ',x)
	x=ifelse(nchar(x)>maxlen,paste0(substr(x,1,maxlen-2),'..'),x)
	x
}
