args <- commandArgs(TRUE) ;
fname <- args[1] ;
datadir<- args[2] ;
library("ape") ;
source("cSIRio.R") ;
cSIR_loaddata(fname,datadir) ;
q()

