#include <R.h>
#include <Rmath.h>


void sample_SIR(int *NH, int *niters, double *Q, int *S, int *I)
{
  
  int t, i, j, *It,  NHosts=*NH, jNH;
  double sump , *p ;
  p = calloc(sizeof(double),NHosts+1) ;
  It = calloc(sizeof(int),NHosts+1) ;
  for (t = 1; t < 50 ; t++)
  {
       for (j = 0; j < NHosts ; j++)
       {
	   		sump = 0 ;
	   		p[NHosts]=1 ;
	   		jNH = j*NHosts ; 
	 		for ( i = 0; i < NHosts ; i++)
	  		{
	      		// probability of infection between i and j
	      	 
	      		p[i] = 1-pow(Q[i + jNH],I[t-1 + i**niters]) ;
	      		p[NHosts] = p[NHosts]*pow(Q[i + jNH],I[t-1 + i**niters]) ;
	      		sump += p[i] ;
	  		}
	  		sump += p[NHosts] ;
	  		
	  		for ( i = 0 ; i<= NHosts ; i++)
	  		{
	        	p[i] = p[i]/sump ;	
	    	}
	    	
	  		rmultinom(S[j],p, NHosts+1, It) ; 	  		
	  		I[t + j**niters]=0 ; 
	  		
	  		for ( i = 0; i<NHosts; i++ )
	  		{
	  			I[t + j**niters]+= It[i] ;
	  		} 
	  		
	  		S[j] = It[NHosts] ;
       }
       
       
  }
  free(It) ;
  free(p) ;
}