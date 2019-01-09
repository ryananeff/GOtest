GOtest=function(x, go, query.population=NULL, background='common', name.x='Input', name.go='Category', method=c('hypergeometric','GSEA','logitreg'), adj="BH", ReportOverlap=TRUE, ncores=1, gsea.alpha=1, permutations=ifelse(method=='GSEA',1000,0), iseed=12345){
#x, a data.frame with first column being genes and second column being groups, and for GSEA method, the third column being measure of association strenght
#go, a data.frame with first column being genes and second column being GO categories. The GO category is written in a format of System:Term.
#query.population, a vector of background genes in where x is derived
	method=match.arg(method)
	if(is.null(name.go) || is.na(name.go) || name.go=='') stop('Invalid argument "name.go"\n')
	if(is.null(name.x) || is.na(name.x) || name.x=='') stop('Invalid argument "name.x"\n')
	## clean up data
	go=cleanData(go,ncols=2)
	if(method=='hypergeometric'){
		x=cleanData(x,ncols=2)
		if(is.character(background)){
			background=match.arg(background,c('common','query','annotation','intersection','union'))
		}
		if(is.null(query.population)){
			if(is.character(background)){
				if(background=='common' || background=='query') stop(paste0('query.population can not be NULL when background="',background,'"\n'))
			}
		}else{
			query.population=query.population[! is.na(query.population)]
			query.population=query.population[query.population != '']
			query.population=unique(query.population)
		}
		if(is.numeric(background)){
			n1=length(unique(union(x[,1],go[,1])))
			if(n1>background) stop('background is smaller than the size of the union of query and annotated gene sets\n')
			BackgroundSize=background
		}else if (background=='annotation'){
			x=x[x[,1] %in% go[,1],,drop=FALSE]
			BackgroundSize=length(unique(go[,1]))
		}else{
			if(background=='common' || background=='intersection'){
				query.population=intersect(query.population,go[,1])
				go=go[go[,1] %in% query.population,,drop=FALSE]
				x=x[x[,1] %in% query.population,,drop=FALSE]
				BackgroundSize=length(query.population)
			}else if(background=='union'){
				query.population=union(query.population,go[,1])
				go=go[go[,1] %in% query.population,,drop=FALSE]
				x=x[x[,1] %in% query.population,,drop=FALSE]
				BackgroundSize=length(query.population)
			}else if(background=='query'){
				go=go[go[,1] %in% query.population,,drop=FALSE]
				x=x[x[,1] %in% query.population,,drop=FALSE]
				BackgroundSize=length(query.population)
			}else{
				stop('Invalid background\n')
			}
		}
		if(ncores>1){
			cl=makeCluster(ncores)
			cat(ncores,'CPU cores recruited\n')
			gps1=unique(x[,2])
			gps2=unique(go[,2])
			if(length(gps1)<length(gps2)){
				EnrichTable <- parLapply(cl,gps1,function(x,x0,go1,BackgroundSize,ReportOverlap){
					x1=x0[x0[,2]==x,,drop=FALSE]
					SetOverlapHypGeoTest(Set1=go1, Set2=x1, totalBalls=BackgroundSize, ReportOverlapElements=ReportOverlap)
				},x0=x,go1=go,BackgroundSize=BackgroundSize,ReportOverlap=ReportOverlap)
			}else{
				EnrichTable <- parLapply(cl,gps2,function(x,x1,go0,BackgroundSize,ReportOverlap){
					go1=go0[go0[,2]==x,,drop=FALSE]
					SetOverlapHypGeoTest(Set1=go1, Set2=x1, totalBalls=BackgroundSize, ReportOverlapElements=ReportOverlap)
				},x1=x,go0=go,BackgroundSize=BackgroundSize,ReportOverlap=ReportOverlap)
			}
			stopCluster(cl)
			EnrichTable=do.call(rbind,EnrichTable)
		}else{
			EnrichTable <- SetOverlapHypGeoTest(Set1=go, Set2=x, totalBalls=BackgroundSize, ReportOverlapElements=ReportOverlap)
		}
	}else if(method=='GSEA'){
		x=cleanData2(x,ncols=3)
		EnrichTable <- callGSEAfunc(Set1=go, Set2=x, alpha=gsea.alpha, permutations=permutations, ncores=ncores, iseed=iseed)
	}else{
		x=cleanData2(x,ncols=ncol(x))
		EnrichTable <- logitregtest(Set1=go, Set2=x, permutations=permutations, ncores=ncores, iseed=iseed)
	}
	if(method=='hypergeometric'){
		cnames0=c(name.go,name.x,'Overlap.Size',paste0(name.x,'.Size'),paste0(name.go,'.Size'),'Background.Size','FE','Pvalue')
		colnames(EnrichTable)=if(ReportOverlap){c(cnames0,'Overlap.Items')}else{cnames0}
		EnrichTable=data.frame(EnrichTable,stringsAsFactors=FALSE)
		if(ReportOverlap) EnrichTable=data.frame(EnrichTable[,- ncol(EnrichTable)],P.adj=NA,Overlap.Items=EnrichTable[,ncol(EnrichTable)],stringsAsFactors=FALSE)
		else EnrichTable$P.adj=NA
		for(i in c('Overlap.Size',paste0(name.x,'.Size'),paste0(name.go,'.Size'),'Background.Size','FE','Pvalue')) EnrichTable[,i]=as.numeric(EnrichTable[,i])
	}else if(method=='GSEA'){
		colnames(EnrichTable)=c(name.go,name.x,'Overlap.Size',paste0(name.x,'.Size'),paste0(name.go,'.Size'),'ES','NES','Pvalue','P.adj','Overlap.Items')
		if(! ReportOverlap) EnrichTable=EnrichTable[,- ncol(EnrichTable)]
	}else{
		colnames(EnrichTable)=c(name.go,name.x,'Overlap.Size',paste0(name.x,'.Size'),paste0(name.go,'.Size'),'log.odds','SE','Z','Pvalue','P.adj','Overlap.Items')		
		if(! ReportOverlap) EnrichTable=EnrichTable[,- ncol(EnrichTable)]
	}
	cnames0=colnames(EnrichTable)
	EnrichTable$System=sub('^([^:]+):(.+)$','\\1',EnrichTable[,name.go])
	EnrichTable[,name.go]=sub('^([^:]+):(.+)$','\\2',EnrichTable[,name.go])
	EnrichTable=EnrichTable[order(EnrichTable[,'Pvalue']),]
	if(method=='hypergeometric' || (method=='logitreg' & permutations==0)) EnrichTable$P.adj=p.adjust(EnrichTable[,'Pvalue'],method=adj)
	return(EnrichTable[,c('System',cnames0)])
}
## clean up data
cleanData <- function(x,ncols=2) {
	## convert datasets to characters
	x=matrix(apply(x[,1:ncols],2,as.character),ncol=ncols);
	## remove NA, and ""
	x[x==""]=NA
	for(i in 1:ncols) x <- x[! is.na(x[,i]), ,drop=FALSE]
	if(nrow(x)==0) stop("No valid data left after removing NAs\n")
	## remove redundant rows
	unique(x)
}
cleanData2 <- function(x,ncols=3) {
	stopifnot(is.data.frame(x))
	stopifnot(is.numeric(x[,3]))
	## convert datasets to characters
	for(i in 1:2) x[,i]=as.character(x[,i]);
	## remove NA, and ""
	x[x[,1]=="",1]=NA
	x[x[,2]=="",2]=NA
	for(i in 1:ncols) x <- x[! is.na(x[,i]), ,drop=FALSE]
	if(nrow(x)==0) stop("No valid data left after removing NAs\n")
	## remove redundant rows
	unique(x)
}
