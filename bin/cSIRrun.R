args <- commandArgs(TRUE) ;
outdir <- args[1]
datadir <- args[2] 
paramfile <- args[3]

library("ape") 
library("phangorn") 
source("cSIRmain.R") 
source("phylo.R")
source(paramfile)
dyn.load("cSIR.so") 
set.seed(opt$seed)

# load in data
load(file=paste(datadir,"/dat.RData",sep="")) ;
load(file=paste(datadir,"/seqs.RData",sep="")) ;

if (!is.null(opt$firstN))
{
  dat$SN<-dat$SN[opt$firstN] ;
  dat$ST<-dat$ST[opt$firstN] ;
  dat$info<-dat$info[opt$firstN] ;
}

opt$outdir <- outdir
smp <- cSIR_runmcmc(seqs, dat, opt, init, mcp, hp)
vars<-c("mr","dr","B","ll","tr")
for (i in seq(1,length(vars)))
{
	assign(vars[i],smp[[vars[i]]])
	save(list=vars[i],file=paste(outdir,"/",vars[i],".mcmc",sep="")) ;
}
if (mcmc.params$acc.rate==1)
{
  acc_rate <- smp$acc.rate
  save(acc_rate, file=paste(outdir,"/acc_rate.mcmc",sep=""))
}  
warnings()
#q()

