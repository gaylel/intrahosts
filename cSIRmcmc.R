mcmc_loadsmps<-function(resdir)
{
	library("mcmc")
	load(paste(resdir,"/out.mcmc",sep=""))
	mr<-smp$mr
	mr<-as.vector(unlist(mr))
	
	NHosts<-ncol(smp$B[[1]])
	
	
}

mcmc_convert<-function(ch)
{
	mcmc_ch<-mcmc(ch)
	return(mcmc_ch)
}

mcmc_stats<-function()
{

}

mcmc_getBij<-function(i,j,B,NHosts)
{
	k<-(j-1)*NHosts + i
	NHosts2=NHosts*NHosts
	b<-(unlist(B))[seq(k,NHosts2*length(B),by=NHosts2)]	
	b<-mcmc(b)
	return(b)
}

mcmc_writetrees<-function(trees, resdir)
{
	lab<-paste("STATE",seq(1:length(trees))-1,sep="_")
	names(trees)<-lab
	write.nexus(trees,file=paste(resdir,"/out.trees",sep=""))
}

