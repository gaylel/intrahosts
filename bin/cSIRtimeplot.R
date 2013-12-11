args <- commandArgs(TRUE) ;
simdir <- args[1]

hosts <- seq(1,6)
r <- matrix(0,10,length(hosts))
for (i in hosts) 
{
	r[,i] <- as.matrix(read.table(file=paste(simdir,"/simtest.",i,"/run.time",sep="")))
}

print(r)
fout <- paste(simdir,"/run.time.eps",sep="")
postscript(fout, horizontal = FALSE)
boxplot(as.data.frame(r),names=hosts,xlab="# hosts (50 sequences each)", ylab="Time to run 1000 iterations (s)")
dev.off()

