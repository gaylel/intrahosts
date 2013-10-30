psumexpdist<-function(r,t)
{
	p=0 ;
	N<- length(r) ;
	rdiff<-r ;
	for (j in seq(1,N))
	{
		rdiff[j]<-prod(r[-j]-r[j]);
	}
	
	for (j in seq(1,N))
	{
		p=p+exp(-r[j]*t)/rdiff[j] ;
	}
	p=p*prod(r) ;
	return(p) ;
}


cSIRcalcc1<-function(m,w,a,b,N)
{
	if (a==(N-1) && b==1)
	{
			
	}
}

metroph<-function(params,ll.fun,prior.fun,maxit=1000)
{
	nparams<-length(params) ;
	chain<-matrix(NA,nrow=maxit,ncol=nparams) ;
	ll<-vector(length=maxit) ;
	last <- params ;
	chain[1,]<-last ;
	last.ll<- ll.fun()
}


cSIRtransprob<-function(r,s,t,N,br,dr)
{
	# relative removal rate
	rho<-dr/br ;
	
	A<-factorial(N-1) ;
	A<-A*(rho^(N-r-s)) ;
	A<-A/(factorial(r)*factorial(s)) ;
	
	
	
	
}