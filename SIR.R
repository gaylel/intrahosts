reed_frost<-function(I0=1,N=100, q=0.5,niters=100)
{
  I<-I0 ;
  S<-N-I ;
  R<-0 ;
  
  for (n in seq(1,niters))
  {
  p<-1-q^I ;
  I<-rbinom(1,S,p) ; 
  S<-S-I ;
  print(I) ;
  }
}



reed_frost_multi<-function(I0=1, N=100, NHosts=2, q0=0.99, q1=0.99)
{
  # evolution of viral population sizes in multiple hosts
  niters<-50;
  I<-matrix(0,nrow=50,ncol=NHosts) ;
  I[1,1]<-I0 ;
  It<-matrix(0,nrow=NHosts, ncol=NHosts+1) ;
  
  S <- matrix(N,nrow=1,ncol=NHosts) ;
  S[1,1] <- N - I[1,1] ;
  p <- matrix(0,nrow=1,ncol=NHosts+1) ;
 
  for (n in seq(2,niters))
  {
       for (j in seq(1,NHosts))
       {
	 p[1,NHosts+1]=1 ;
	 for ( i in seq(1,NHosts))
	  {
	      # probability of infection between i and j
	      q<-ifelse(i==j,q0,q1) ; 
	      p[1,i]<-1-q^I[n-1,i] ;
	      p[1,NHosts+1] = p[1,NHosts+1]*(q^I[n-1,i]) ;
	  }
	  p<-p/sum(p) ;
	  print(p) ;
	  It[j,]<-rmultinom(1,S[1,j],p) ; 
	  
	  I[n,j]<- sum(It[j,1:NHosts]);
       }
       

       
       for ( i in seq(1,NHosts))
       {
	  S[1,i]<-It[i,NHosts+1];
       }
       
       #print(S[1,]);
       #print(I[n,]) ;
       #print(It);
  }
  return(list(I=I));
}

reed_frost_multi.2<-function(I0=1, N, NHosts, q0=0.99, q1=0.99)
{
  # evolution of viral population sizes in multiple hosts
  niters<-50;
  I<-matrix(0,nrow=50,ncol=NHosts) ;
  I[1,1]<-I0 ;
  Q <- matrix(q1, nrow=NHosts,ncol = NHosts) ;
  for (i in seq(1,NHosts)){
  		Q[i,i] <-q0 ;
  	
  } 
  S <- matrix(N,nrow=1,ncol=NHosts) ;
  S[1,1] <- N - I[1,1] ;
  p <- matrix(0,nrow=1,ncol=NHosts+1) ;
  .C("sample_SIR2",as.integer(NHosts),as.integer(niters),as.double(as.vector(Q)),as.integer(as.vector(S)),as.integer(as.vector(I))) ;
  print(S) ;
  return(list(I=I));
}

reed_frost_multi_demo<-function(N=100, NHosts=10)
{
	dyn.load("SIR.so") ;
	rf<-reed_frost_multi.2(N,NHosts,q0=0.999,q1=0.9999) ;
	return(list(I=rf$I)) ;
}

plot_reed_frost<-function(rf,psfile)
{
  postscript(psfile,horizontal=FALSE)
  iters<- nrow(rf$I) ;
  nHosts<-6#ncol(rf$I) ;
  par(mfrow=c(nHosts,1)) ;

  for (i in seq(1,nHosts))
  {
    plot(1:iters,rf$I[,i],type="b",col="black",xlab="t",ylab=i) ;
  }
  dev.off() 
}
