args <- commandArgs(TRUE) ;
outdir<-args[1]
datadir<- args[2] ;
library("ape") ;
library("phangorn") ;
source("cSIRmain.R") ;
dyn.load("cSIR.so") ;

# load in data
load(file=paste(datadir,"/dat.RData",sep="")) ;
load(file=paste(datadir,"/seqs.RData",sep="")) ;

#dat$SN<-10 ;
#dat$ST<-dat$ST[1] ;
#dat$info[[1]]$key<-dat$info[[1]]$key[1:10] ;

<<<<<<< HEAD
dat$SN<-dat$SN[1:20] ;
dat$ST<-dat$ST[1:20] ;
dat$info<-dat$info[1:20] ;
=======
dat$SN<-dat$SN[1:10] ;
dat$ST<-dat$ST[1:10] ;
dat$info<-dat$info[1:10] ;
>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714
#dat$info<-dat$info[[1]] ;

NHosts<-length(dat$SN) ;

# Initialisation
# convert sequences to class PhyBin
#x<-phyDat(seqs[1:sum(dat$SN)]) ;


# calculate UPGMA tree for each within-host population


<<<<<<< HEAD
system.time(smp<-cSIR_SB_metrop2(nHosts=NHosts,I0=1,nS=1000,dr=1,dat=dat,x=seqs,N=20))
=======
system.time(smp<-cSIR_SB_metrop2(nHosts=NHosts,I0=1,nS=1000,dr=1,dat=dat,x=seqs,N=100))
>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714
save(smp,file=paste(outdir,"/out.mcmc",sep="")) ;
warnings()
q()

