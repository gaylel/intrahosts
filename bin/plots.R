skyride_plots<-function(ph,sr,psout=out.ps)
{
  postscript(psout,horizontal=FALSE)
  par(mfrow=c(2,1))
  phpr<-plot.phylo(ph,show.tip.label=FALSE,edge.width=1)
  par(bty="n")
  plot(sr[,1],sr[,3],type="l",xlim=phpr$x.lim,ylim=c(0,0.12),xlab="",ylab="N_E")
  #polygon(c(sr[,1],rev(sr[,1])),c(sr[,4],rev(sr[,5])),col="grey",border= NA)
  lines(sr[,1],sr[,4],lty=2)
  lines(sr[,1],sr[,5],lty=2)
  lines(sr[,1],sr[,3],lwd=2)
  dev.off()
}