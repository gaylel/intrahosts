args <- commandArgs(TRUE) 
rundir <- args[1] 
library("ape") 
library("phangorn") 
source("cSIRmain.R") 
source("cSIRsim.R")
dyn.load("cSIR.so") 
source("cSIRio.R")
paramfile <- "sim.params" 
params <- cSIRsimtree(paste(rundir,"/",paramfile,sep=""))
# save params
save(params,file=paste(rundir,"/data/params.RData", sep="")) 

# save tree in nwk format
write.tree(params$tr_list$tr, paste(rundir,"/data/test.nwk", sep=""))
q()

