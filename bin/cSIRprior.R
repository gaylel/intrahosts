cSIR_loadoptsparams <- function(opfile)
{
  source(opfile)
}

prior_init<-function(nHosts,dat)
{
	SN<-dat$SN 
	ST<-dat$ST ;
	
	
	# initialise the parameters
	dr=0.1 ;
	mr=1e-3 ;
	B<-matrix(0.001*runif(nHosts*nHosts),nrow=nHosts,ncol=nHosts) ;
	for (i in seq(1,nHosts))
	{
		B[i,i]=1 ;
	}
	Bcon<-cSIR_Bstruct(ST=ST) ;
	B<-B*Bcon ;
	
	nB<-length(Bcon>0) ;
	l1=1 ;
	l2=1/0.1;
	l3=1/0.1;
	
	iend<-NULL ;
	Bchain<-list() ;
	drchain<-list() ;
	sirchain<-list() ;
	trchain<-list() ;
	llchain<-list() ;
	mrchain<-list() ;
	dvec<-NULL
	drvec<-NULL
	pold<-cSIR_Bdrprior(B,Bcon,dr,l1,l2,l3)  ;
	M=4
	
	
	# initialise 
	ini<-cSIR_upgma(x,SN,ST)
	lo<-ini$lo ;
	s<-ini$s ;
	x<-phyDat(x) ;
	llcur<- log(0) ;
	params<-list(B=B, dr=dr, mr= mr, ll=llcur, lo=lo, s=s, Bcon=Bcon)
	
	return(params)
}