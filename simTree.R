args <- commandArgs(TRUE) ;
datadir <- args[1] ;
library("ape") ;
library("phangorn") ;
source("cSIRmain.R") ;
dyn.load("cSIR.so") ;
source("cSIRio.R") ;
cSIR_testtree(datadir)
q()

