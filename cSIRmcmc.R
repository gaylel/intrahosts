mcmc_loadsmps<-function(resdir)
{
	library("mcmc")
	load(paste(resdir,"/out.mcmc",sep=""))
<<<<<<< HEAD
	mr<-smp$mr
	mr<-as.vector(unlist(mr))
	
	NHosts<-ncol(smp$B[[1]])
	
	
=======
	return(smp)
>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714
}

mcmc_convert<-function(ch)
{
	mcmc_ch<-mcmc(ch)
	return(mcmc_ch)
}

<<<<<<< HEAD
mcmc_stats<-function()
{

=======
mcmc_trace<-function(smp,param,i,j)
{
    ch<-smp[[param]]
    NHosts<-ncol(smp$B[[1]])
    if (param=="B")
    {
      ch<-mcmc_getBij(i,j,ch,NHosts)
    }
    else{
      ch<-as.vector(unlist(ch))
    }
    #ch<-mcmc_convert(ch)
    return(ch)  
>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714
}

mcmc_getBij<-function(i,j,B,NHosts)
{
<<<<<<< HEAD
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

=======
  k<-(j-1)*NHosts + i
  NHosts2=NHosts*NHosts
  b<-as.vector(unlist(B))[seq(k,length(B),by=NHosts2)]
  return(b)
}

mcmc_writetrees<-function(trees,resdir)
{
  lab<-paste("STATE",seq(1:length(trees))-1,sep="_")
  names(trees)<-lab
  write.nexus(trees,file=paste(resdir,"/out.trees",sep=""))
}








>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714
