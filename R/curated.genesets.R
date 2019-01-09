curated.genesets=function(sets=c('MacArthur')){
	gmtfiles=c(MacArthur='MacArthur.gmt',HGNC_universe='HGNC_universe.gmt')
	if(! all(sets %in% names(gmtfiles))) stop('Invalid sets\n')
	gmtfiles=gmtfiles[sets]
	obj=NULL
	for(s in sets){
		f=file.path(system.file("extdata", package="GOtest"),gmtfiles[s])
		res=read.gmt(f)
		res$geneset.names=paste('Curated',res$geneset.names,sep=':')
		if(is.null(obj)){
			obj=res
		}else{
			obj$genesets=c(obj$genesets,res$genesets)
			obj$geneset.names=c(obj$geneset.names,res$geneset.names)
			obj$geneset.description=c(obj$geneset.description,res$geneset.description)
		}
	}
	names(obj$genesets)=obj$geneset.names
	obj
}
