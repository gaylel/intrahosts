cSIRsim<-function(NHosts=3, N=1000, SN=c(10,10,10), ST=c(7,11,16), rootname="test")
{
	# simulate data from the model
	params<-list()
	dat.params<-list()
	dat.params$I0<-1
	dat.params$NS<-N
	dat.params$NHosts<-NHosts
	params$B<-matrix(0,nrow=NHosts,ncol=NHosts)
	for (i in seq(1,NHosts))
	{
		params$B[i,i] <- 1
		
	}
	params$B[1,2]<-0.4
	params$B[2,3]<-0.2
	params$dr<-0.1
	params$s<-NULL
	params$lo<-NULL
	a2<-0
	for (i in seq(1,NHosts))
	{
		a1<-a2+1
		a2<-a2+SN[i]
		params$lo<-c(params$lo,rep(i,SN[i]))
	}
	params$lo<-rbind(seq(1,sum(SN)),params$lo)
	
	s<-seq(1:(sum(SN)-1))
	a2<-0
	for (i in seq(1,NHosts-1))
	{
		a1<-a2+1
		a2<-a2+SN[i]
		params$s<-c(params$s,sample(s[a1:a2]))
	}
	params$s<-c(params$s,sample(s[(a2+1):length(s)]))
	
	#params$s<-rev(params$s)
	params$T<-NULL
	print(params)
	dat<-list()
	dat$SN<-SN
	dat$ST<-ST
	for (i in seq(1,NHosts))
	{
		dat$info[[i]]<-list(key=paste("H",i,"S",seq(1,SN[i]),sep="")) ;
	}
	tre<-cSIR_drawTre(params, dat.params, dat, Ntries=1000)
	return(list(tre=tre,params=params, dat.params=dat.params, dat=dat))
}


cSIRsimtree <- function(paramfile)
{
	source(paramfile)
	x <- NULL
	
	# initialise sampling times and number of hosts.
	set.seed(opt$seed)
	NHosts <- length(dat$SN)
	SN<-dat$SN 
  	ST<-dat$ST 
  	dat.params<-list(I0=init$I0, NS=init$NS, NHosts=length(SN))
  	dat$STin = dat$ST - min(dat$ST)
	params <- cSIR_modelinit(x, init, dat, dat.params, mcp)
	params$tr_list$tr$tip.label <- NULL
	params$tr_list$tr$node.label <- NULL
	for (i in seq(1,NHosts))
	{
		params$tr_list$tr$tip.label <- c(params$tr_list$tr$tip.label, paste("H",i,"S",seq(1,SN[i]),sep="")) ;
		
	}
	return(params)
	
}

