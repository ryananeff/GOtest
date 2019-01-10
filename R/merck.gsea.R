#Author: Minghui Wang
merck.genesets=function(sets=c('GO CC','GO MF','GO BP','KEGG PATHWAYS','INGENUITY PATHWAYS','GENEGO PATHWAYS'), species=c('human','mouse','rat')){
	sets=match.arg(sets,several.ok=TRUE)
	gof=file.path(system.file("extdata", package="GOtest"),paste0('merck.db.',species,'.RDS'))
	go=readRDS(gof)
	go=go[go[,2] != 'GO CC:chromosome',]
	go=go[sub('(.+):(.+)','\\1',go[,2]) %in% sets,]
	go
}
merck.collections=c('GO CC','GO MF','GO BP','KEGG PATHWAYS','INGENUITY PATHWAYS','GENEGO PATHWAYS')
merck.gsea=function(x, query.population=NULL, genesets=c('GO CC','GO MF','GO BP','KEGG PATHWAYS','INGENUITY PATHWAYS','GENEGO PATHWAYS'),
 background='common', name.x='Input', name.go='Merck', method=c('hypergeometric','GSEA','logitreg'),
 adj="BH", species=c('human','mouse','rat'), ReportOverlap=TRUE, ncores=1, gsea.alpha=1, permutations=ifelse(method=='GSEA',1000,0), iseed=12345){
	species=match.arg(species)
	genesets=match.arg(genesets,several.ok=TRUE)
	go=merck.genesets(sets=genesets,species=species)
	GOtest(x=x,go=go,query.population=query.population,background=background,name.x=name.x,name.go=name.go,method=method,adj=adj,ReportOverlap=ReportOverlap,ncores=ncores,gsea.alpha=gsea.alpha,permutations=permutations, iseed=iseed)
}
