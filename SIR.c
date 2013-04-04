#include <R.h>
#include <Rmath.h>

void sample_SIR(int *NHosts, double **Q, double **S, double **I)
{
  
  int t, niters, ;
  for (n in seq(2,niters))
  {
       for (j in seq(1,NHosts))
       {
	 p[1,NHosts+1]=1 ;
	 for ( i in seq(1,NHosts))
	  {
	      // probability of infection between i and j
	      // q<-ifelse(i==j,q0,q1) ; 
	      p[1,i]<-1-Q[i,j]^I[n-1,i] ;
	      p[1,NHosts+1] = p[1,NHosts+1]*(Q[i,j]^I[n-1,i]) ;
	  }
	  p<-p/sum(p) ;
	  //print(p) ;
	  It[j,]<- rmultinom(1,S[1,j],p) ; 
	  I[n,j]<- sum(It[j,1:NHosts]);
       }
       
       for ( i in seq(1,NHosts))
       {
	  S[1,i]<-It[i,NHosts+1];
       }
       
       //print(S[1,]);
       //print(I[n,]) ;
       //print(It);
  }
}