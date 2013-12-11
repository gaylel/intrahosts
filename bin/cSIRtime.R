args <- commandArgs(TRUE) ;
outdir <- args[1]
datadir <- args[2] 
paramfile <- args[3]

library("ape") 
library("phangorn") 
source("cSIRmain.R") 
source("phylo.R")
dyn.load("cSIR.so") 

N <- 10
s.ts <- matrix(0, nrow=N, ncol=1)
for (i in seq(1, N))
{
	s.t <- system.time(cSIR_runmodel(outdir, datadir, paramfile))
	s.ts[i, 1] <- s.t[3]
}

write.table(s.ts, file = paste(outdir, "/run.time", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
