#Compute pathway crosstalk as described Jia et al BMC Syst Biol. 2011; 5(Suppl 3): S12. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3287567/)
JaccardCoefficient=function(x,y){
	length(intersect(x,y))/length(union(x,y))
}
OverlapCoefficient=function(x,y){
	length(intersect(x,y))/min(length(unique(x)),length(unique(y)))
}
crosstalk=function(x, col.signature, col.pathway='MSigDB', col.overlap.size='Overlap.Size', col.adj.p='P.adj', col.elements='Overlap.Items', sep=';', adj.p.cutoff=0.05, min.size=5, min.common=3,...){
	x = x[x[,col.overlap.size] >= min.size & x[,col.adj.p] <= adj.p.cutoff, ]
	res=lapply(split(x,x[,col.signature]),function(y,col.signature,col.pathway, col.overlap.size, col.adj.p, col.elements, sep=';',adj.p.cutoff=0.05,min.size=5,min.common=3){
		.crosstalk(y, col.signature=col.signature, col.pathway=col.pathway, col.overlap.size=col.overlap.size, col.adj.p=col.adj.p, col.elements=col.elements, sep=sep,adj.p.cutoff=adj.p.cutoff,min.size=min.size,min.common=min.common)
	},col.signature=col.signature, col.pathway=col.pathway, col.overlap.size=col.overlap.size, col.adj.p=col.adj.p, col.elements=col.elements, sep=sep,adj.p.cutoff=adj.p.cutoff,min.size=min.size,min.common=min.common)
	res=do.call(rbind,res)
	if(is.null(res)) return(NULL)
	res$MeanScore=(res$JC+res$OC)/2
	res=res[order(-res$MeanScore),]
	rownames(res)=NULL
	res
}
.crosstalk=function(x, col.signature, col.pathway, col.overlap.size, col.adj.p, col.elements, sep=';', adj.p.cutoff=0.05, min.size=5, min.common=3, ...){
	#
	#(1) only pathways with at least min.size (=5) DE elements
	#(2) only pathways with adjusted P values < adj.p.cutoff (=0.05)
	#(3) two pathways in crosstalk were required to share at least min.common (3) DE elements
	#
	x = x[x[,col.overlap.size] >= min.size & x[,col.adj.p] <= adj.p.cutoff, ]
	if(nrow(x)<2) return(NULL)
	if(length(unique(x[,col.signature]))>1) stop('Multiple signatures found in input. Please call crosstalk(..) instead for signature specific crosstalk.\n')
	res=setNames(lapply(x[,col.elements],function(z) strsplit(z,sep)[[1]]),x[,col.pathway])
	tab=NULL
	for(i in 1L:(length(res)-1)){
		s1=unique(res[[i]])
		for(j in (i+1):length(res)){
			s2=unique(res[[j]])
			JC=length(intersect(s1,s2))/length(union(s1,s2))
			OC=length(intersect(s1,s2))/min(length(s1),length(s2))
			N=length(intersect(res[[i]],res[[j]]))
			tab=rbind(tab,data.frame(Pathway1=names(res)[i],Pathway2=names(res)[j],N1=length(s1),N2=length(s2),NCommon=N,JC=JC,OC=OC,stringsAsFactors=FALSE))
		}
	}
	tab=data.frame(Signature=x[1,col.signature],tab,stringsAsFactors=FALSE)
	tab=tab[tab$NCommon >= min.common,]
	tab
}
