mcmc_loadsmps<-function(resdir)
{
	
	# load variable smp
	load(paste(resdir,"/out.mcmc",sep=""))

	mr<-smp$mr
	mr<-as.vector(unlist(mr))
	
	NHosts<-ncol(smp$B[[1]])
	return(smp)
}

mcmc_loadvar <- function(resdir, var, opt, mcp)
{
  if (opt$saveevery > 0)
  {
    i <- 1
    smp <- NULL
    while (i*opt$saveevery < mcp$Niters)
    {
      load(paste(resdir, "/", var, "_", format(i*opt$saveevery, scientific=FALSE), ".mcmc", sep=""))
      smp<-c(smp,eval(parse(text=var)))
      i <- i + 1
    }
    
  }else{
    load(paste(resdir, "/", var, ".mcmc", sep=""))  
    smp<-eval(parse(text=var))
  }
  smp[[1]] <- NULL
  return(smp)
}

mcmc_loadsmp<-function(resdir,var, opt, mcp)
{
	# when smps are stored individually
	
	vn<-list()
	switch(var,
	"B"={
	  smp <- mcmc_loadvar(resdir, var, opt, mcp)
    NHosts<-ncol(smp[[1]])
    smp<-t(matrix(unlist(smp),nrow=NHosts*NHosts))
    
		
		k=1
		for (i in seq(1,NHosts))
		{
			for (j in seq(1,NHosts))
			{
				vn[[k]]<-paste("B_",j,",",i,sep="")
				k<-k+1
			}
		}
		smp<-mcmc_convert(smp)
	},
	"tr"={
		library("ape")
		smp <- mcmc_loadvar(resdir, var, opt, mcp)
		# mcmc_writetrees(smp[5000:length(smp)], resdir)	
		#tr<-read.nexus(paste(resdir,"/out.trees",sep=""))
		tr<-smp
		#smp<-matrix(0,ncol=length(tr[[1]]$tip.label)-1,nrow=length(tr)/100)
		smp<-NULL
		for (i in seq(1,length(tr),by=100))
		{
			#smp[i,]<-tr[[i]]$edge.length[order(tr[[i]]$edge[,1])]
			#smp[i,]<-branching.times(tr[[i]])
			smp<-rbind(smp,branching.times(tr[[i]]))
			#print(smp)
		}
		
		for (i in seq(1,ncol(smp)))
		{
			vn[[i]]<-paste("bt_",i,sep="")

		}	
		
		smp<-mcmc_convert(smp,100)
		#smp<-mcmc_convert(smp)
	},
  "tr_out"={
    library("ape")
    smp <- mcmc_loadvar(resdir, "tr", opt, mcp)
    smp <- smp[50000:length(smp)]
    mcmc_writetrees(smp[5000:length(smp)], resdir)
    smp <- NULL
  },       
	{
	  smp <- mcmc_loadvar(resdir, var, opt, mcp)
	  smp<-matrix(unlist(smp),ncol=1)
    print(smp)
		vn[[1]]<-var
		smp<-mcmc_convert(smp)
	}
	)
	
	
	varnames(smp)<-vn 
	#smp<-mcmc_process(smp,1000,100)
	#smp<-mcmc_process(smp,3000,1)
	smp<-mcmc_process(smp,0,1)
	return(smp)
}

mcmc_convert<-function(ch,thin=1)
{
	mcmc_ch<-mcmc(ch,thin=thin)
	
	return(mcmc_ch)
}

mcmc_process<-function(ch,burnin,thin)
{
	smp<-window(x=ch,thin=thin,start=burnin+1)
	return(smp)
}

mcmc_tredgelengths<-function(tr)
{
	

}

mcmc_stats<-function()
{

}

mcmc_trace<-function(param)
{

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
	b<-(unlist(B))[seq(k,NHosts2*length(B),by=NHosts2)]	
	return(b)
}

mcmc_writetrees<-function(trees, resdir)
{
	lab<-paste("STATE",seq(1:length(trees))-1,sep="_")
	names(trees)<-lab
	write.nexus(trees,file=paste(resdir,"/out.trees",sep=""))
}

