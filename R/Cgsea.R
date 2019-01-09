fGSEAfunc=function(x,alpha=1,isOrdered=FALSE){
	if(!isOrdered) x=x[order(- x$strength),]
	.C("fGSEAfunc2",as.integer(x$tag),as.numeric(x$strength),nrow(x),as.numeric(alpha),as.numeric(0.0),PACKAGE='GOtest')[[5]]
}
fGSEA=function(x, go, alpha=1, permutations=1000,ncores=1,iseed=12345){
#x, a data.frame with first column gene ids and second column phenotype association measures
#go, a vector of gene ids
	x1=data.frame(strength=x[,2],tag=x[,1] %in% go,stringsAsFactors=FALSE)
	ntag=sum(x1$tag)
	if(ntag<2){
		Res=list(ES=0,indicator=which(x1$tag),RES=rep(0,nrow(x1)),NES=0,Pvalue=1)
		class(Res)='gsea'
		return(Res)
	}
	x1=x1[order(- x1$strength),] #we will keep the order from now on
	isOrdered=TRUE
	Res=GSEAfunc(x1, alpha=alpha, isOrdered=isOrdered)
	x1$tag=FALSE #reset for repeated sampling in permutations
	func1=function(i,x1,ntag,alpha,isOrdered=FALSE){
		x1[sample.int(nrow(x1),ntag,replace=FALSE),2]=TRUE #we will keep the order of column "strength" such that there is no need to sort it again
		fGSEAfunc(x1,alpha=alpha,isOrdered=isOrdered)
	}
	if(ncores>1){
		cl=makeCluster(ncores)
		cat('Recruited',ncores,'cores\n')
		clusterSetRNGStream(cl, iseed)
		clusterExport(cl,varlist=c('fGSEAfunc'),envir=environment())
		ESnull=parSapply(cl,1:permutations,func1,x1=x1,ntag=ntag,alpha=alpha,isOrdered=isOrdered)
		stopCluster(cl)
	}else{
		set.seed(iseed)
		ESnull=sapply(1:permutations,func1,x1=x1,ntag=ntag,alpha=alpha,isOrdered=isOrdered)
	}
	if(Res$ES>=0){
		ESnull=ESnull[ESnull>=0]
		Res$NES=Res$ES/mean(ESnull)
		Res$Pvalue=sum(ESnull >= Res$ES)/length(ESnull)
	}else{
		ESnull=ESnull[ESnull < 0]
		Res$NES=Res$ES/abs(mean(ESnull))
		Res$Pvalue=sum(ESnull <= Res$ES)/length(ESnull)
	}
	class(Res)='gsea'
	Res
}
fcallGSEAfunc=function(Set1, Set2, alpha=1, permutations=1000,ncores=1,iseed = 12345){
#Set1, usually Gene Ontology
#Set2, data.frame of phenotype association signatures
	stopifnot(is.data.frame(Set2))
	colnames(Set2)=c('tag','phenotype','strength')
	gos=split(Set1[,1],Set1[,2])
	Res=NULL
	if(ncores>1){
		cl=makeCluster(ncores)
		cat('Recruit',ncores,'cores\n')
		clusterSetRNGStream(cl, iseed)
	}else{
		set.seed(iseed)
	}
	fun1=function(k,dat1,gos,alpha,isOrdered=FALSE){
		#dat1[,3]=sample(dat1[,3],replace=FALSE) #sample the 1st column rather than the 3rd column to avoid re-sorting in the later step when isOrdered=TRUE
		dat1$tag=sample(dat1$tag,replace=FALSE)
		sapply(gos,function(iGO,dat1,alpha,isOrdered=FALSE){
			dat1$tag=dat1$tag %in% iGO
			fGSEAfunc(dat1,alpha=alpha,isOrdered=isOrdered)
		},dat1=dat1,alpha=alpha,isOrdered=isOrdered)
	}
	for(dat in split(Set2,Set2$phenotype)){ #for each phenotype
		phenotype=dat$phenotype[1]
		dat=dat[order(- dat$strength),colnames(dat)!='phenotype'] #we will sort the data by the strength and keep the order from now on
		isOrdered=TRUE
		tab2=lapply(gos,function(iGO,dat,alpha,isOrdered=FALSE){
			Overlap=dat$tag
			dat$tag = dat$tag %in% iGO
			ESs=0
			if(any(dat$tag)) ESs=fGSEAfunc(dat,alpha=alpha,isOrdered=isOrdered)
			Overlap=sort(Overlap[dat$tag])
			data.frame(Overlap.Size=length(Overlap),Input.Size=nrow(dat),Category.Size=length(iGO),ES=ESs,NES=NA,Pvalue=NA,P.adj=NA,Elements=paste0(Overlap,collapse=';'),stringsAsFactors=FALSE)
		},dat=dat,alpha=alpha,isOrdered=isOrdered)
		tab2=do.call(rbind,tab2)
		#Compute significance and control for multiple tests by permutations
		if(ncores>1){
			ESnull=parLapply(cl,1:permutations,fun1,dat1=dat,gos=gos,alpha=alpha,isOrdered=isOrdered)
		}else{
			ESnull=lapply(1:permutations,fun1,dat1=dat,gos=gos,alpha=alpha,isOrdered=isOrdered)
		}
		ESnull=do.call(cbind,ESnull) #number of GO terms x number of permutations
		NESnull=vector('list',nrow(tab2))
		for(i in 1:nrow(tab2)){
			if(tab2$ES[i]==0){
				NESnull[[i]]=0
				tab2$NES[i]=0
				#tab2$Pvalue[i]=1
				next
			}else if(tab2$ES[i]<0){
				s1=ESnull[i,ESnull[i,]<0]
				#tab2$Pvalue[i]=(sum(s1 <= tab2$ES[i])+1)/ (length(s1)+1)
			}else{
				s1=ESnull[i,ESnull[i,]>0]
				#tab2$Pvalue[i]=(sum(s1 >= tab2$ES[i])+1) / (length(s1)+1)
			}
			if(length(s1)==0) next #this happens when there is insufficient number of permutations
			m1=abs(mean(s1))
			NESnull[[i]]=s1/m1
			tab2$NES[i]=tab2$ES[i]/m1
		}
		NESnull=unlist(NESnull)
		#Method 1: separate positive and negative scores (no longer used)
		#NESnull=sort(NESnull)
		#wpid = tab2$NES >= 0
		#if(nrow(tab2)>1){
			#if(any(wpid)) tab2$P.adj[wpid]=getFDR.0(tab2$NES[wpid],NESnull[NESnull >= 0],truncated.size=0)
			#if(any(! wpid)) tab2$P.adj[ ! wpid ]=getFDR.0(abs(tab2$NES[ ! wpid]),abs(NESnull[NESnull < 0]),truncated.size=0)
		#}
		#Method 2: combine positive and negative
		NESnull=abs(NESnull)
		pq=getFDR.0(abs(no.na(tab2$NES)),NESnull,truncated.size=0,return.p.value=TRUE)
		tab2$Pvalue[!is.na(tab2$NES)]=pq$P.value
		tab2$P.adj[!is.na(tab2$NES)]=pq$P.adj
		Res=rbind(Res,data.frame(Term=names(gos),Input=phenotype,tab2,stringsAsFactors=FALSE))
	}
	if(ncores>1) stopCluster(cl)
	colnames(Res)=c('Category','Input','Overlap.Size','Input.Size','Category.Size','ES','NES','Pvalue','P.adj','Elements')
	rownames(Res)=NULL
	Res
}
