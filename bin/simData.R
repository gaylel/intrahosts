args <- commandArgs(TRUE) 
rundir <- args[1] 
library("ape") 
library("phangorn") 
source("cSIRmain.R") 
source("cSIRsim.R")
dyn.load("cSIR.so") 
source("cSIRio.R")
paramfile <- "sim.params" 
source(paste(rundir,"/",paramfile,sep=""))

NHosts <- length(dat$SN)
for (i in seq(1,NHosts))
{
	dat$info[[i]]<-list(key=paste("H",i,"S",seq(1,dat$SN[i]),sep="")) ;
}
# save params
save(dat,file=paste(rundir,"/data/dat.RData", sep="")) 

# load dna data

seqs <- dnaBINmxtolist(paste(rundir,"/data", sep=""))
save(seqs,file=paste(rundir,"/data/seqs.RData", sep="")) 
q()

