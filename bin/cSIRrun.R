args <- commandArgs(TRUE) ;
outdir <- args[1]
datadir <- args[2] 
paramfile <- args[3]

library("ape") 
library("phangorn") 
source("cSIRmain.R") 
source("phylo.R")
dyn.load("cSIR.so") 



cSIR_runmodel(outdir, datadir, paramfile)

#q()
