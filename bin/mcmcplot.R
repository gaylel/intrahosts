args <- commandArgs(TRUE) ;
smpdir <- args[1]
outdir <- args[2]
params <- args[3]

library("ape") 
library("coda") 
source("cSIRmcmc.R") 
source(params)
vars<-c("mr","B","ll","tr_out")
vars<-c("bn")


	
for (i in seq(1,length(vars)))
{
  pdf(paste(outdir,"/", vars[i], ".pdf",sep=""))
	smp<-mcmc_loadsmp(smpdir,vars[i], opt, mcp)
	plot(smp)
	rm(smp)
  dev.off()
}
#smp<-mcmc_loadsmp(smpdir,"tr")
#rm(smp)
warnings()
q()