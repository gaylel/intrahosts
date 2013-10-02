mcmc_loadsmps<-function(resdir)
{
	library("mcmc")
	load(paste(resdir,"/out.mcmc",sep=""))
	return(smp)
}

mcmc_convert<-function(ch)
{
	mcmc_ch<-mcmc(ch)
	return(mcmc_ch)
}

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
}

mcmc_getBij<-function(i,j,B,NHosts)
{
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








