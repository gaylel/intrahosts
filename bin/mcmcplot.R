args <- commandArgs(TRUE) ;
smpdir<-args[1]
outdir<-args[2]


library("ape") 
library("coda") 
source("cSIRmcmc.R") 

#vars<-c("mr","dr","B","ll","tr")
vars<-c("B")


pdf(paste(outdir,"/traces.pdf",sep=""))
	
for (i in seq(1,length(vars)))
{
	smp<-mcmc_loadsmp(smpdir,vars[i])
	plot(smp)
	rm(smp)
}
dev.off()
#smp<-mcmc_loadsmp(smpdir,"tr")
#rm(smp)
warnings()
q()