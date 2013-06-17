#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

SEXP sample_cSIR_R(SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr) ;
SEXP sample_cSIR_S_R(SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN) ;
void cSIR_iter(int n, int nS, int NHosts, double *B, double dr, int **p1, double **p2) ;
void cSIR_iter_ST(int *np, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN) ;
void cSIR_iters(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2) ;
void cSIR_iters_ST(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN) ;
int check_infectives(int *I, int NHosts, int t) ;
int check_sampled(double *ST, int NHosts, int t_n, int t_n1 ) ;

SEXP sample_cSIR_R(SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr)
{
	
	int I0, NS, NHosts, *S, *I, *T, *R, **p_RVAL1 = Calloc(2,int*),*n=Calloc(1,int),i,j;
	int *PR_I, *PR_S;
	double *B, **p_RVAL2 = Calloc(1,double*), dr, *PR_T ;
	
	
	SEXP R_Bdim, R_I, R_T, R_S, R_list ;
	
	R_I0 = coerceVector(R_I0, INTSXP) ;
	R_NS = coerceVector(R_NS, INTSXP) ;
	R_NHosts = coerceVector(R_NHosts, INTSXP) ;
	R_B = coerceVector(R_B, REALSXP) ;
	R_dr = coerceVector(R_dr, REALSXP) ;
	
	I0 = INTEGER(R_I0)[0] ;
	NS = INTEGER(R_NS)[0] ;
	NHosts = INTEGER(R_NHosts)[0] ;
	B = REAL(R_B) ;
	dr = REAL(R_dr)[0] ;

	R_Bdim = getAttrib(R_B, R_DimSymbol) ;
	n[0]=0 ;
	
	cSIR_iters(n, I0, NS, NHosts, B, dr, p_RVAL1, p_RVAL2) ;
	
	
	
	
	PROTECT(R_I=allocMatrix(INTSXP,n[0],NHosts)) ;
	PROTECT(R_S=allocMatrix(INTSXP,n[0],NHosts)) ;
	PROTECT(R_T=allocMatrix(REALSXP,n[0],4)) ;
	PR_I = INTEGER(R_I) ;
	PR_S = INTEGER(R_S) ;
	PR_T = REAL(R_T) ;
	
	for (i=0; i<n[0]; i++)
	{
		for (j=0; j< NHosts ; j++)
		{
			PR_I[i+j*n[0]]=p_RVAL1[1][i*NHosts + j] ;
			PR_S[i+j*n[0]]=p_RVAL1[0][i*NHosts + j] ;
		}
	}
	
	for (i=0; i<(n[0]) ; i++)
	{
		for (j=0; j<4; j++)
		{
			PR_T[i+j*n[0]] = p_RVAL2[0][i*4+j] ;
		}
	}
	
	
	PROTECT(R_list = allocVector(VECSXP ,3)) ;
	SET_VECTOR_ELT(R_list, 0, R_I) ;
	SET_VECTOR_ELT(R_list, 1, R_S) ;
	SET_VECTOR_ELT(R_list, 2, R_T) ;
	UNPROTECT(4) ;
	return(R_list);
}

SEXP sample_cSIR_S_R(SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN)
{
	
	int I0, NS, NHosts, *S, *I, *T, *R, **p_RVAL1 = Calloc(2,int*),*n=Calloc(1,int),i,j;
	int *PR_I, *PR_S, *SN;
	double *B, **p_RVAL2 = Calloc(1,double*), dr, *PR_T, *ST ;
	
	
	SEXP R_Bdim, R_I, R_T, R_S, R_list ;
	
	// initial number of infective particles
	R_I0 = coerceVector(R_I0, INTSXP) ;
	
	// number of particles in each host	
	R_NS = coerceVector(R_NS, INTSXP) ;
	
	// number of hosts
	R_NHosts = coerceVector(R_NHosts, INTSXP) ;
	
	// matrix of transmission coefficients
	R_B = coerceVector(R_B, REALSXP) ;
	
	// death rate
	R_dr = coerceVector(R_dr, REALSXP) ;
	
	// vector of sampling times
	R_ST = coerceVector(R_ST, REALSXP) ;
	
	// vector of viral population sizes
	R_SN = coerceVector(R_SN, INTSXP) ;
	
	
	I0 = INTEGER(R_I0)[0] ;
	NS = INTEGER(R_NS)[0] ;
	NHosts = INTEGER(R_NHosts)[0] ;
	B = REAL(R_B) ;
	dr = REAL(R_dr)[0] ;
	ST = REAL(R_ST) ;
	SN = INTEGER(R_SN) ;


	R_Bdim = getAttrib(R_B, R_DimSymbol) ;
	n[0]=0 ;
	
	//cSIR_iters(n, I0, NS, NHosts, B, dr, p_RVAL1, p_RVAL2) ;
	cSIR_iters_ST(n, I0, NS, NHosts, B, dr, p_RVAL1, p_RVAL2, ST, SN) ;
	PROTECT(R_I=allocMatrix(INTSXP,n[0],NHosts)) ;
	PROTECT(R_S=allocMatrix(INTSXP,n[0],NHosts)) ;
	PROTECT(R_T=allocMatrix(REALSXP,n[0],4)) ;
	PR_I = INTEGER(R_I) ;
	PR_S = INTEGER(R_S) ;
	PR_T = REAL(R_T) ;
	
	for (i=0; i<n[0]; i++)
	{
		for (j=0; j< NHosts ; j++)
		{
			PR_I[i+j*n[0]]=p_RVAL1[1][i*NHosts + j] ;
			PR_S[i+j*n[0]]=p_RVAL1[0][i*NHosts + j] ;
		}
	}
	
	for (i=0; i<(n[0]) ; i++)
	{
		for (j=0; j<4; j++)
		{
			PR_T[i+j*n[0]] = p_RVAL2[0][i*4+j] ;
		}
	}
	
	
	PROTECT(R_list = allocVector(VECSXP ,3)) ;
	SET_VECTOR_ELT(R_list, 0, R_I) ;
	SET_VECTOR_ELT(R_list, 1, R_S) ;
	SET_VECTOR_ELT(R_list, 2, R_T) ;
	UNPROTECT(4) ;
	return(R_list);
}

int check_infectives(int *I, int NHosts, int t)
{
	int i=0, ep_end=1;
	while (ep_end==1 && i<NHosts)
	{
		if (I[t*NHosts + i]>0) ep_end=0;
		i++ ;
	}
	return(ep_end);
}

int check_sampled(double *ST, int NHosts, int t_n, int t_n1 )
{
	// return index of host which was sampled first in the interval (t_n, t_n1)
	int i, rv=-1 ;
	double tmin=t_n1+1;
	
	
	for (i=0 ; i<NHosts ; i++)
	{
		if (ST[i] > t_n && ST[i]<=t_n1 && ST[i] <= tmin)
		{
			tmin=ST[i] ;
			rv=i ;
		}
	}
	return(rv) ;
}

void cSIR_iters(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2)
{
	int *pp1[2];
	double *pp2[1] ;
	
	pp1[0]= Calloc(NHosts,int) ;
	pp1[1]= Calloc(NHosts,int) ;
	pp2[0]= Calloc(4,double) ;
	
	int i ;
	for (i=0; i<NHosts ; i++)
	{
		pp1[0][i] = nS ;
		pp1[1][i] = 0 ;
	}
	
	pp1[0][0]-=I0 ;
	pp1[1][0]+=I0 ;
	
	i=1 ;
	while (check_infectives(pp1[1],NHosts,i-1)==0 && i<100)
	{
		cSIR_iter(i, nS, NHosts, B, dr, pp1, pp2) ;
		i++ ;
	}
	p1[0]=pp1[0] ;
	p1[1]=pp1[1] ;
	p2[0]=pp2[0] ;
	pp1[0]=0;
	pp1[1]=0;
	pp2[0]=0;
	n[0]=i	;
}

void cSIR_iters_ST(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN)
{
	int *pp1[2];
	double *pp2[1] ;
	
	pp1[0]= Calloc(NHosts,int) ;
	pp1[1]= Calloc(NHosts,int) ;
	pp2[0]= Calloc(4,double) ;
	
	int i ;
	for (i=0; i<NHosts ; i++)
	{
		pp1[0][i] = nS ;
		pp1[1][i] = 0 ;
	}
	
	pp1[0][0]-=I0 ;
	pp1[1][0]+=I0 ;
	
	i=1 ;
	n[0]=1 ;
	
	while (check_infectives(pp1[1],NHosts,n[0]-1)==0 && n[0]<10000)
	{
		cSIR_iter_ST(n, nS, NHosts, B, dr, pp1, pp2, ST, SN) ;
		//cSIR_iter(i, nS, NHosts, B, dr, pp1, pp2) ;
		//i++ ;
		n[0]++ ;
	}
	p1[0]=pp1[0] ;
	p1[1]=pp1[1] ;
	p2[0]=pp2[0] ;
	pp1[0]=0;
	pp1[1]=0;
	pp2[0]=0;
	//n[0]=i	;
}


void cSIR_time_samples(double *tpts)
{
	// sample at 
}

void cSIR_iter(int n, int nS, int NHosts, double *B, double dr, int **p1, double **p2)
{
	double *R, sumR=0, t_n, Rtot, *T, rand_e;	// rate matrix
	R = Calloc(NHosts * (NHosts+1),double) ;
	int i,j,e,ec,er, *S, *I;
	
	S=p1[0];
	I=p1[1];
	T=p2[0];
	
	for (i=0; i<NHosts; i++)
	{
		R[i + NHosts*NHosts]=I[(n-1)*NHosts + i]*dr ;
		sumR+= R[i + NHosts*NHosts] ;
		
		for (j=0; j<NHosts; j++)
		{
			R[i*NHosts + j]=I[(n-1)*NHosts + j]*S[(n-1)*NHosts + i]*B[i*NHosts + j]/(double)(nS) ;
			sumR+= R[i*NHosts + j] ;
		}
	}
	
	for (i=0; i<(NHosts*(NHosts+1));i++)
	{
		R[i]=R[i]/sumR ;
	}
	
	// Draw time of next event
	t_n = rexp(1.0/sumR) ;
	
	// Draw event that happens next from R
	rand_e = unif_rand();
	
	e=0;
	Rtot=R[0] ;
	while (Rtot < rand_e)
	{
		e++;
		Rtot+=R[e] ;
	}	
	
	ec =(int) (floor((double)e/NHosts)) ;
	er =(int) (e-ec*NHosts);
	
	// increase size of S,I,T
	S=Realloc(S,(n+1)*NHosts ,int) ;
	I=Realloc(I,(n+1)*NHosts ,int) ;
	T=Realloc(T,(n+1)*4 ,double) ;
	
	for (i=0 ; i< NHosts; i++)
	{
		S[i+(n*NHosts)] = S[i+((n-1)*NHosts)] ;
		I[i+(n*NHosts)] = I[i+((n-1)*NHosts)] ;
		
	}
	
	
	if (ec == NHosts)
	{
		I[n*NHosts + er] = I[n*NHosts + er] - 1;
		T[n*4+3] = -1 ;
		T[n*4+1] = er ;
		T[n*4+2] = er ;
		
	}
	else
	{
		I[n*NHosts + ec] = I[n*NHosts + ec] + 1;
		S[n*NHosts + ec] = S[n*NHosts + ec] - 1;
		T[n*4+3] = 1.0 ;
		T[n*4+1] = er ;
		T[n*4+2] = ec ;
	}
	T[n*4] = t_n+T[(n-1)*4];
	p1[0]=S ;
	p1[1]=I ;
	p2[0]=T ;
	S=0;
	I=0;
	T=0;
}

void cSIR_iter_ST(int *np, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN)
{
	double *R, sumR=0, t_n, Rtot, *T, rand_e, *ST2, T_n;	// rate matrix
	R = Calloc(NHosts * (NHosts+1),double) ;
	ST2 = Calloc(NHosts,double) ;
	int i,j,e,ec,er, *S, *I, t_n1, n, h_i;
	
	S=p1[0];
	I=p1[1];
	T=p2[0];
	n=np[0] ;
	
	for (i=0; i<NHosts; i++)
	{
		ST2[i] = ST[i] ;
		R[i + NHosts*NHosts]=I[(n-1)*NHosts + i]*dr ;
		sumR+= R[i + NHosts*NHosts] ;
		
		for (j=0; j<NHosts; j++)
		{
			R[i*NHosts + j]=I[(n-1)*NHosts + j]*S[(n-1)*NHosts + i]*B[i*NHosts + j]/(double)(nS) ;
			sumR+= R[i*NHosts + j] ;
		}
	}
	
	for (i=0; i<(NHosts*(NHosts+1));i++)
	{
		R[i]=R[i]/sumR ;
	}
	
	// Draw time to next event
	t_n = rexp(1.0/sumR) ;
	
	// Check to see if any host has been sampled
	T_n = T[(n-1)*4] ;
	t_n1 = t_n + T_n ;
	
	h_i = 1 ;
	while (h_i >= 0)
	{
		h_i = check_sampled(ST2, NHosts, T_n, t_n1 ) ;
		if (h_i>=0)
		{
			
			// update t_n
			//t_n	= ST[h_i] ;
			ST2[h_i] = DBL_MAX ;
			// update S, I T
			S=Realloc(S,(n+1)*NHosts ,int) ;
			I=Realloc(I,(n+1)*NHosts ,int) ;
			T=Realloc(T,(n+1)*4 ,double) ;
	
			for (i=0 ; i< NHosts; i++)
			{
				S[i+(n*NHosts)] = S[i+((n-1)*NHosts)] ;
				I[i+(n*NHosts)] = I[i+((n-1)*NHosts)] ;		
			}
			I[h_i+(n*NHosts)] = 0 ;
			S[h_i+(n*NHosts)] = 0 ;
			T[n*4+3] = - I[h_i+((n-1)*NHosts)];
			T[n*4+1] = h_i ;
			T[n*4+2] = h_i ;	
			T[n*4] = ST[h_i] ;				
			n++ ;
		}
	}
	
	sumR=0 ;
	for (i=0; i<NHosts; i++)
	{
		R[i + NHosts*NHosts]=I[(n-1)*NHosts + i]*dr ;
		sumR+= R[i + NHosts*NHosts] ;
		
		for (j=0; j<NHosts; j++)
		{
			R[i*NHosts + j]=I[(n-1)*NHosts + j]*S[(n-1)*NHosts + i]*B[i*NHosts + j]/(double)(nS) ;
			sumR+= R[i*NHosts + j] ;
		}
	}
	
	for (i=0; i<(NHosts*(NHosts+1));i++)
	{
		R[i]=R[i]/sumR ;
	}
	// Draw event that happens next from R
	rand_e = unif_rand();
	
	e=0;
	Rtot=R[0] ;
	while (Rtot < rand_e)
	{
		e++;
		Rtot+=R[e] ;
	}	
	
	ec =(int) (floor((double)e/NHosts)) ;
	er =(int) (e-ec*NHosts);
	
	// increase size of S,I,T
	S=Realloc(S,(n+1)*NHosts ,int) ;
	I=Realloc(I,(n+1)*NHosts ,int) ;
	T=Realloc(T,(n+1)*4 ,double) ;
	
	for (i=0 ; i< NHosts; i++)
	{
		S[i+(n*NHosts)] = S[i+((n-1)*NHosts)] ;
		I[i+(n*NHosts)] = I[i+((n-1)*NHosts)] ;
		
	}
	
	
	if (ec == NHosts)
	{
		I[n*NHosts + er] = I[n*NHosts + er] - 1;
		T[n*4+3] = -1 ;
		T[n*4+1] = er ;
		T[n*4+2] = er ;
		
	}
	else
	{
		I[n*NHosts + ec] = I[n*NHosts + ec] + 1;
		S[n*NHosts + ec] = S[n*NHosts + ec] - 1;
		T[n*4+3] = 1.0 ;
		T[n*4+1] = er ;
		T[n*4+2] = ec ;
	}
	T[n*4] = t_n+T[(n-1)*4];
	p1[0]=S ;
	p1[1]=I ;
	p2[0]=T ;
	S=0;
	I=0;
	T=0;
	np[0] = n ;
}