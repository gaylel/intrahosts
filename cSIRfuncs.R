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