args <- commandArgs(TRUE) 
opt<-args[1]
testdir<-args[2]

#library("ape") 
#library("phangorn") 
source("cSIRmain.R") 
source("phylo.R")
dyn.load("cSIR.so") 
set.seed(2000)

SN<-c(10,10,10)
ST<-c(7,9,11)
dat<-list()
dat$SN<-SN
dat$ST<-ST
dat$info<-list()
NHosts<-3
dat$info[[1]]<-list(key=as.character(seq(1,10)))
dat$info[[2]]<-list(key=as.character(seq(11,20)))
dat$info[[3]]<-list(key=as.character(seq(21,30)))



s<-seq(29,1,by=-1)
lo<-rbind(seq(1,30),c(rep(1,10),rep(2,10),rep(3,10)))

switch(opt,
	"tree.1"={
		sir<-cSIR_demo_1()
		save(sir,file=paste(testdir,"/tree.1.RData",sep=""))
		print(sir)
	},
	"tree.2"={
		load(file=paste(testdir,"/tree.1.RData",sep=""))
		library("ape")
		pdf(paste(testdir,"/tree.2.pdf",sep=""))
		tr_list<-.Call("tree_reconstruct",sir, NHosts, dat) 
		save(tr_list,file=paste(testdir,"/tree.2.RData",sep=""))
		tr_list[[1]]$node.label=as.character(seq(31,59))
		plot.phylo(tr_list[[1]],show.node.label=TRUE)
		dev.off()
		print(tr_list)
	},
	"tree.3"={
		load(file=paste(testdir,"/tree.2.RData",sep=""))
		library("ape")
		pdf(paste(testdir,"/tree.3.pdf",sep=""))
		T=tr_list[[2]]
	
		# get new branching order based on current and leaf order
		bt<-cSIR_stimes(T=T,lo=lo,s=s) ;
	
	
	
		tr2<-cSIR_phylofroms(s=bt$s,lo=bt$lo, T=T, dat=dat) ;
		tr_list2<-list(tr=tr2,s=bt$s, lo=bt$lo, T=T) ;
		
		
		save(tr_list2,file=paste(testdir,"/tree.3.RData",sep=""))
		tr2$node.label=as.character(seq(31,59))
		plot.phylo(tr2,show.node.label=TRUE)
		dev.off()
		print(tr_list2)
	},
	"tree.4"={
		#
		load(file=paste(testdir,"/tree.3.RData",sep=""))
		library("ape")
		pdf(paste(testdir,"/tree.4.pdf",sep=""))
		T<-tr_list2$T
		lo<-tr_list2$lo
		s<-tr_list2$s
		# get new branching order based on current and leaf order
		bt<-cSIR_stimes(T=T,lo=lo,s=s) ;
	
	
	
		tr2<-cSIR_phylofroms(s=bt$s,lo=bt$lo, T=T, dat=dat) ;
		tr_list3<-list(tr=tr2,s=bt$s, lo=bt$lo, T=T) ;
		
		
		save(tr_list3,file=paste(testdir,"/tree.4.RData",sep=""))
		tr2$node.label=as.character(seq(31,59))
		plot.phylo(tr2,show.node.label=TRUE)
		dev.off()
		print(tr_list3)
	},
	"tree.5"={
		load(file=paste(testdir,"/tree.3.RData",sep=""))
		library("ape")
		library("phangorn") 
		pdf(paste(testdir,"/tree.5.pdf",sep=""))
		T<-tr_list2$T
		lo<-tr_list2$lo
		tr<-tr_list2$tr
		s<-tr_list2$s
		for (i in seq(1,2000))
		{
			#res<-cSIR_narrowexchange(tr,SN, lo, s)
			res<-cSIR_ne(tr,SN,lo,s)
			lo<-res$lo
			tr<-res$tr
			s<-res$s
			print(lo)
			#res$tr$node.label=as.character(seq(31,59))
			#res$tr$node.label=as.character(lo[2,match(seq(31,59),lo[1,])])
			#print(res$tr$node.label)
			#print(tr)
			plot.phylo(res$tr,show.node.label=TRUE)
		}
		dev.off()
	},
	"tree.B"={
		source("cSIRmain.R")
		library("ape")
		library("phangorn") 
		load(file=paste(testdir,"/tree.3.RData",sep=""))
		Bchain<-list() ;
		drchain<-list() ;
		sirchain<-list() ;
		trchain<-list() ;
		llchain<-list() ;
		mrchain<-list() ;
		T<-tr_list2$T
		print(T)
		lo<-tr_list2$lo
		tr<-tr_list2$tr
		s<-tr_list2$s
		dr=1 ;
		mr=1.83e-4 ;
		B<-matrix(0.1*runif(NHosts*NHosts),nrow=NHosts,ncol=NHosts) ;
		for (i in seq(1,NHosts))
		{
			B[i,i]=2 ;
		}
		Bcon<-cSIR_Bstruct(ST=ST) ;
		print(Bcon)
		B<-B*Bcon ;
		l1=1/2 ;
		l2=1/0.5;
		l3=1/0.5;
		hyperparams<-list(l1=l1,l2=l2,l3=l3) 
		dat.params<-list(I0=1, NS=1000, NHosts=3,x=NULL)
		params<-list(B=B, dr=dr, mr= mr, ll=1, lo=lo, s=s, Bcon=Bcon, tr=tr)
	 	params$T<-diff(T[,1])
	 	TS<-cSIR_getTstats(T,3)
		#print(TS)
		params$TS<-TS
		Nacc<-1
		for (i in seq(1,1))
		{
			params<-cSIR_Bupdate(params, dat.params, dat, hyperparams, th=0.3)
			trchain[[Nacc]]<-params$tr ;
			llchain[[Nacc]]<-params$ll ;
			drchain[[Nacc]]<-params$dr ;
			Bchain[[Nacc]]<-params$B ;
			mrchain[[Nacc]]<-params$mr ;
			#sirchain[[Nacc]]<-sir ;
			
			print(sprintf("%i loglikelihood   %8.4f   Mutation rate   %8.4f  Death rate %8.4f Move %i",Nacc,params$ll,params$mr,params$dr, 1))
			#print(params$T)
			if ((Nacc %% 10)==0)
			{
					plot.phylo(params$tr,show.node.label=TRUE) ;
			}
			
			#print(params)
			Nacc<-Nacc+1 ;
		}
		smp<-list(B=Bchain,dr=drchain,mr=mrchain,tr=trchain,ll=llchain) ;
		vars<-c("mr","dr","B","ll","tr")
		for (i in seq(1,length(vars)))
		{
			assign(vars[i],smp[[vars[i]]])
			save(list=vars[i],file=paste(testdir,"/",vars[i],".mcmc",sep="")) ;
		}
				
	}	
	)
#q()