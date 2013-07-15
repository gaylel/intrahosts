# input / output functions for the model

cSIR_readmetadata<-function(fname)
{
	# read in metadata
	dat<-read.table(file=fname,fill=TRUE,skip=1) ;

	# Just use Newmarket data?
	dat<-dat[which(dat[,5]=="Newmarket"),] ;
	
	
	NHosts<-nrow(dat) ;
	info<-list() ;
	
	seqs<-gsub("HE","",dat[,9]) ;
	seqs<-strsplit(seqs,"-") ;
	seqs<-matrix(as.numeric(unlist(seqs)),ncol=2,byrow=TRUE) ;

	for (i in seq(1,NHosts))
	{
		l=seqs[i,2]-seqs[i,1] ;
		info[[i]]<-list(isolate=dat[i,8],acc=paste("HE",seq(seqs[i,1],seqs[i,2]-1),sep=""),key=paste("H",i,"S",seq(1,l),sep="")) ;
		
	}
	
	
	SN<-seqs[,2]-seqs[,1];
	SN[SN<1]<-1 ;

	dts<-dat[,3] ;
	ms<-which(is.na(dts==TRUE)) ;
	dts[ms]<-dat[ms,3] ;
	dts<-strptime(dts,format="%Y-%m-%d",tz="GMT") ;
	ST<-c(0,difftime(dts[2:length(dts)],dts[1],units="days")) ;
	ST<-ST+7 ;	
	ret<-list(SN=SN,ST=ST,info=info) ;
	return(ret) ;
}

cSIR_getsequences<-function(dat)
{
	# take in metadata structure and then get sequences from GenBank
	NHosts<-length(dat$info) ;
	acc.names<-NULL ;
	seq.names<-NULL ;
	for (i in seq(1,NHosts))
	{
		acc.names<-c(acc.names,dat$info[[i]]$acc) ;
		seq.names<-c(seq.names,dat$info[[i]]$key) ;
		
	}
	
	x<-read.GenBank(access.nb=acc.names,seq.names=seq.names) ;
	return(x) ;
}

cSIR_loaddata<-function(fname,datadir)
{
	dat<-cSIR_readmetadata(fname) ;
	seqs<-cSIR_getsequences(dat) ;
	
	save(dat,file=paste(datadir,"/dat.RData",sep="")) ;
	save(seqs,file=paste(datadir,"/seqs.RData",sep="")) ;
}