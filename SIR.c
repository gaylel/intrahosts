#include <R.h>
#include <Rinternals.h>

#include <Rmath.h>

void sample_SIR(int *NH, int *niters, double *Q, int *S, int *I) ;
int * sample_SIR2(int *NH, int *niters, double *Q, int *S, int *I, int **I2, int *Irows) ;
SEXP sample_SIR2_R(SEXP R_I0, SEXP R_N, SEXP R_NHosts, SEXP R_Q) ;



void sample_SIR(int *NH, int *niters, double *Q, int *S, int *I)
{
  
  int t, i, j, *It,  NHosts=*NH, jNH;
  double sump , *p ;
  p = Calloc(NHosts+1, double) ;
  It = Calloc(NHosts+1, int) ;
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

int* sample_SIR2(int *NH, int *niters, double *Q, int *S, int *I, int **I2, int *Irows)
{
  // also return I2, a matrix of interhost infections, where each row is:
  // time	from i	to j	k
  
  Rprintf("%i,%i,%8.4f,%i\n",NH[0],niters[0],Q[0],I[0]) ;
  int t, i, j, ind, NHosts=*NH, jNH, Itcount=0, **rval, *It2;
  double sump ;
  double p[NHosts+1] ;
  int It[NHosts+1] ;
  It2 = Calloc(4, int) ;
  
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
	  	    Rprintf("%i\n", It[0]) ;
	  		
	  		for ( i = 0 ; i< NHosts ; i++)
	  		{
	  			if (i!=j && It[i]>0)
	  			{
	  				ind = Itcount*4 ; 
	  				Itcount++;
	  				It2 = Realloc(It2, 4*Itcount, int) ;
	  				It2[ind] = t ;
	  				It2[ind + 1] = i ; 
	  				It2[ind + 2] = j ;
	  				It2[ind + 3] = It[i] ;
	  			}
	  		}
	  		 	  		
	  		I[t + j**niters]=0 ; 
	  		
	  		for ( i = 0; i<NHosts; i++ )
	  		{
	  			I[t + j**niters]+= It[i] ;
	  			
	  		} 
	  		
	  		
	  		S[j] = It[NHosts] ;
       }
       
  }
 // free(It) ;
 // free(p) ;
  I2[0]=It2 ;
  Irows[0] = Itcount ;
  Rprintf("%i\n",Itcount) ;
  return(I) ;
}

SEXP sample_SIR2_R(SEXP R_I0, SEXP R_N, SEXP R_NHosts, SEXP R_Q)
{
	SEXP R_Qdim, R_I, R_I2, R_list ;
	int NHosts, niters=50, I0, N, *Irows=Calloc(1,int), i, j ;
	int *I , *I2, **p_I2 = Calloc(1,int*), *pR_I2, *S, *Iout;
	
	double *Q ;
	
	R_I0 = coerceVector(R_I0, INTSXP) ;
	R_N = coerceVector(R_N, INTSXP) ;
	R_NHosts = coerceVector(R_NHosts, INTSXP) ;
	R_Q = coerceVector(R_Q, REALSXP) ;
	
	I0 = INTEGER(R_I0)[0] ;
	N = INTEGER(R_N)[0] ;
	NHosts = INTEGER(R_NHosts)[0] ;
	Q = REAL(R_Q) ;

	R_Qdim = getAttrib(R_Q, R_DimSymbol) ;
	
	S = Calloc(NHosts, int) ;
	PROTECT(R_I = allocMatrix(INTSXP, niters, NHosts));
	I = INTEGER(R_I) ; 
	
	for (i=0 ; i<NHosts ; i++)
	{
		S[i] = N ;
		for (j=0 ; j<niters ; j++)
		{
			I[i*niters + j]=0;
		}
	}
	I[0] = I0 ;
	S[0] = S[0] - I[0] ; 
	
	Iout = sample_SIR2(&NHosts,&niters, Q, S, I, p_I2, Irows) ;
	
	// Copy return arguments to a list
	PROTECT(R_I2=allocMatrix(INTSXP,Irows[0],4)) ;
	pR_I2=INTEGER(R_I2) ;
	I2 = p_I2[0] ;
	for (i=0; i< Irows[0] ; i++)
	{
		for (j=0; j<4; j++)
		{
			pR_I2[i + j*Irows[0]]=I2[i*4 + j] ;
		}
	}
	
	PROTECT(R_list = allocVector(VECSXP ,2)) ;
	SET_VECTOR_ELT(R_list, 0, R_I) ;
	SET_VECTOR_ELT(R_list, 1, R_I2) ;
	UNPROTECT(3) ;
	Rprintf("%i,%i\n", niters, NHosts) ;
	return(R_list);
}