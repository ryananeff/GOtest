#Logistic regression approach for identifying enrichment
logitregtest <- function(Set1, Set2, permutations=1000, ncores=1, iseed=12345){
#Set1, usually Gene Ontology
#Set2, data.frame of phenotype association signatures
	stopifnot(is.data.frame(Set2))
	colnames(Set2)[1:3]=c('tag','phenotype','strength')
	gos=split(Set1[,1],Set1[,2])
	Res=NULL
	if(ncores>1){
		cl=makeCluster(ncores)
		cat('Recruit',ncores,'cores\n')
	}
	fun2=function(iGO,dat){
		Overlap=dat$tag
		dat$tag = dat$tag %in% iGO
		Pvalue=log.odds=Z=SE=NA
		if(any(dat$tag)){
			fit1 <- glm(tag ~.,family=binomial(link='logit'),data=dat)
			log.odds <- coef(fit1)['strength']
			s1=summary(fit1)
			SE=s1$coefficients['strength','Std. Error']
			Z=s1$coefficients['strength','z value']
			Pvalue=s1$coefficients['strength','Pr(>|z|)']
		}
		Overlap=sort(Overlap[dat$tag])
		data.frame(Overlap.Size=length(Overlap),Input.Size=nrow(dat),Category.Size=length(iGO),log.odds=log.odds,SE=SE,Z=Z,Pvalue=Pvalue,P.adj=NA,Elements=paste0(Overlap,collapse=';'),stringsAsFactors=FALSE)
	}
	fun3=function(k,dat,gos){
		dat$strength=sample(dat$strength,replace=FALSE)
		sapply(gos,function(iGO,dat){
			dat$tag = dat$tag %in% iGO
			if(!any(dat$tag)) return(0)
			fit1 <- glm(tag ~.,family=binomial(link='logit'),data=dat)
			#log.odds <- coef(fit1)['strength']
			s1=summary(fit1)
			#SE=s1$coefficients['strength','Std. Error']
			Z=s1$coefficients['strength','z value']
			Z
		},dat=dat)
	}
	for(dat in split(Set2,Set2$phenotype)){ #for each phenotype
		phenotype=dat$phenotype[1]
		dat=dat[order(- dat$strength),colnames(dat)!='phenotype'] #we will sort the data by the strength and keep the order from now on
		if(ncores>1){
			tab2=parLapply(cl,gos,fun2,dat=dat)
		}else{
			tab2=lapply(gos,fun2,dat=dat)
		}
		tab2=do.call(rbind,tab2)
		if(permutations>0) {
			if(ncores>1){
				clusterSetRNGStream(cl, iseed)
				Zs=parLapply(cl,1:permutations,fun3,dat=dat,gos=gos)
			}else{
				set.seed(iseed)
				Zs=lapply(1:permutations,fun3,dat=dat,gos=gos)
			}
			#Method 2: combine positive and negative
			Zs=unlist(Zs)
			pq=getFDR.0(abs(no.na(tab2$Z)),abs(Zs),truncated.size=0,return.p.value=TRUE)
			tab2$Pvalue[!is.na(tab2$Z)]=pq$P.value
			tab2$P.adj[!is.na(tab2$Z)]=pq$P.adj			
		}
		Res=rbind(Res,data.frame(Term=names(gos),Input=phenotype,tab2,stringsAsFactors=FALSE))
	}
	if(ncores>1) stopCluster(cl)
	colnames(Res)=c('Category','Input','Overlap.Size','Input.Size','Category.Size','log.odds','SE','Z','Pvalue','P.adj','Elements')
	rownames(Res)=NULL
	Res
}
