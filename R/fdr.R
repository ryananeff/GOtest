getFDR.0=function(TS,TSN,truncated.size=0,return.p.value=FALSE){
#TS, a vector of non-negative test statistics from observed data. TS could be a truncated set, with small values discarded to save memory. See 'truncated.size'.
#TSN, a vector of non-negative test statistics from permutation
#truncated.size, number of the smallest values discarded from the original TS
#return.p.value, if TRUE, return a list with both P value and FDR; if FALSE, return a vector of FDR
	if(any(TS<0,na.rm=TRUE) || any(TSN<0,na.rm=TRUE)) stop('Input must be non-negative\n')
	n1=length(TS)
	n2=length(TSN)
	TSN=sort(TSN)
	order.TS=order(TS)
	rank.TS=rep(0L,n1)
	rank.TS[order.TS]=1:n1
	p=rep(0,n1)
	p[order.TS]=n2-findInterval(TS[order.TS],TSN)
	p0=(p+1)/(n2+1) #add 1 to both numerator and denominator to get a more conservative estiamte (this modification was made on May 10, 2018)
	p=p0/((n1+1-rank.TS)/(n1+truncated.size))
	p[p>1]=1
	if(n1==1) return(p)
	for(i in 1:(n1-1)) if(p[order.TS[i+1]] > p[order.TS[i]]){p[order.TS[i+1]]=p[order.TS[i]]} #make sure the larger statistics will have a FDR no bigger than that of the smaller statistics
	if(!return.p.value) return(p)
	list(P.value=p0,P.adj=p)
}
