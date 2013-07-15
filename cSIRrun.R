args <- commandArgs(TRUE) ;
outdir<-args[1]
datadir<- args[2] ;
library("ape") ;
library("phangorn") ;
source("cSIR.R") ;
dyn.load("cSIR.so") ;

# load in data
load(file=paste(datadir,"/dat.RData",sep="")) ;
load(file=paste(datadir,"/seqs.RData",sep="")) ;

NHosts<-length(dat$SN) ;
smp<-cSIR_SB_metrop(nHosts=NHosts,I0=1,nS=1000,dr=1,ST=dat$ST,SN=dat$SN,N=10)
save(smp,file=paste(outdir,"/out.mcmc",sep="")) ;
q()

