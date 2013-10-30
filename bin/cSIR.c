#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "llists.h"
#include "phylo.h"
#include "leaf.h"

SEXP sample_cSIR_R(SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr) ;
SEXP sample_cSIR_S_R(SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN) ;
void cSIR_iter(int n, int nS, int NHosts, double *B, double dr, int **p1, double **p2) ;
void cSIR_iter_ST(int *np, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN) ;
void cSIR_iters(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2) ;
void cSIR_iters_ST(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN) ;
int check_infectives(int *I, int NHosts, int t) ;
int check_sampled(double *ST, int NHosts, double t_n, double t_n1 ) ;
SEXP tree_reconstruct(SEXP R_sir, SEXP R_NHosts, SEXP R_dat) ;


struct hnode {
	item_i *n ;
	item_d *t ;
} ;

typedef struct hnode hnode ;


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
	
	int I0, NS, NHosts, *S, *I, *T, *R, **p_RVAL1 = Calloc(3,int*),*n=Calloc(1,int),i,j;
	int *PR_I, *PR_S, *PR_Iend, *SN;
	double *B, **p_RVAL2 = Calloc(1,double*), dr, *PR_T, *ST ;
	
	
	SEXP R_Bdim, R_I, R_T, R_S, R_Iend, R_list ;
	
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
	PROTECT(R_Iend=allocMatrix(INTSXP,1,NHosts)) ;
	PR_I = INTEGER(R_I) ;
	PR_S = INTEGER(R_S) ;
	PR_T = REAL(R_T) ;
	PR_Iend = INTEGER(R_Iend) ;
	
	for (j=0; j< NHosts ; j++)
	{
		PR_Iend[j] = p_RVAL1[2][j] ;
		for (i=0; i<n[0]; i++)
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
	
	
	PROTECT(R_list = allocVector(VECSXP ,4)) ;
	SET_VECTOR_ELT(R_list, 0, R_I) ;
	SET_VECTOR_ELT(R_list, 1, R_S) ;
	SET_VECTOR_ELT(R_list, 2, R_T) ;
	SET_VECTOR_ELT(R_list, 3, R_Iend) ;
	UNPROTECT(5) ;
	Free(p_RVAL1[0]) ;
	Free(p_RVAL1[1]) ;
	Free(p_RVAL1[2]) ;
	Free(p_RVAL1) ;
	Free(p_RVAL2[0]) ;
	Free(p_RVAL2) ;
	Free(n) ;
	
	return(R_list);
}

//SEXP tree_branchingtimes(SEXP tr, SEXP R_dat)
//{
//	double** phylo_bt(phylo* tr, double *tipinfo, int *hostinfo, int NHosts)	
//}

SEXP tree_reconstruct(SEXP R_sir, SEXP R_NHosts, SEXP R_dat)
{
	int *PR_I, *PR_S, *PR_Iend, *SN, NHosts, i, j, k,st=0, st2=-1, TN, ha, hb, *Tend, *NNodes, *minNodes, mn, r1, r2, *ch, ntips=0, ei, intNode;
	double *PR_T, ***PR_T2, *ST, *PR_bt, **bt ;
	phylo *tr ;
	item_i * tst ;
	
	SEXP R_ST, R_SN, R_Tdim, R_tr, R_T, R_List, R_bt ;
	hnode *Nodes ;
	R_NHosts = coerceVector(R_NHosts, INTSXP) ;
	NHosts = INTEGER(R_NHosts)[0] ;
	R_Tdim = getAttrib(VECTOR_ELT(R_sir,2), R_DimSymbol) ;
	
	TN = INTEGER(R_Tdim)[0] ;
	//Rprintf("dims = %i, %i\n",INTEGER(R_Tdim)[0], INTEGER(R_Tdim)[1]) ;
	
	PR_I = INTEGER(VECTOR_ELT(R_sir,0)) ;
	PR_S = INTEGER(VECTOR_ELT(R_sir,1)) ;	
	PR_T = REAL(coerceVector(VECTOR_ELT(R_sir,2),REALSXP)) ;	
	PR_Iend = INTEGER(VECTOR_ELT(R_sir,3)) ;

	R_SN = VECTOR_ELT(R_dat,0) ;
	R_ST = VECTOR_ELT(R_dat,1) ;
	R_ST = coerceVector(R_ST, REALSXP) ;
	R_SN = coerceVector(R_SN, INTSXP) ;
	ST = REAL(R_ST) ;
	SN = INTEGER(R_SN) ;
	
	for (i=0 ; i<NHosts ; i++)
	{
		ntips+= SN[i] ;
	}
	
	Nodes = Calloc(NHosts,hnode) ; 
	Tend = Calloc(NHosts,int) ;
	NNodes = Calloc(NHosts,int) ;
	minNodes = Calloc(NHosts,int) ;
	ch = Calloc(2,int) ;
	tr=phylo_create(ntips) ;
	ei = tr->Nedge - 1 ;
	intNode = 2*tr->NNode  ;
	
	PR_T2 = Calloc(4,double**) ;
	k=0 ;
	for (j=0 ; j<4 ; j++)
	{
		PR_T2[j] = Calloc(TN,double*) ;
	
		for (i=0 ; i<TN ; i++)
		{	
		 	PR_T2[j][i] = &PR_T[k++];
			//Rprintf("%8.4f, ",*PR_T2[j][i]) ; 	
			
		}
		
	}
	
	
	
	st=0 ;
	for (i=0 ; i<NHosts ; i++)
	{
		NNodes[i] = PR_Iend[i] ;
		//Nodes[i].t = Calloc(PR_Iend[i],double) ;
		//Nodes[i].n = Calloc(PR_Iend[i],int) ;
		
		Nodes[i].t = NULL ;
		Nodes[i].n = NULL ;
		
		
		st2=-1 ;
		minNodes[i] = st ;
		
		for (j=0 ; j<SN[i] ; j++)
		{
			//Nodes[i].n[j] = st ;
			//Nodes[i].t[j] = ST[i] ;
			Nodes[i].n = llist_add_el_i(st, Nodes[i].n) ;
			Nodes[i].t = llist_add_el_d(ST[i], Nodes[i].t) ;
			//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;
			st++ ;	
		}
		
		for (j=SN[i] ; j<PR_Iend[i] ; j++)
		{
			//Nodes[i].n[j] = st2 ;
			//Nodes[i].t[j] = ST[i] ;
			Nodes[i].n = llist_add_el_i(st2, Nodes[i].n) ;
			Nodes[i].t = llist_add_el_d(ST[i], Nodes[i].t) ;
			//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
			st2=st2-1 ;
		}
		
		if (SN[i]<PR_Iend[i])
		{
			minNodes[i] = st2+1 ;  
		}
		//st=st+PR_Iend[i]-SN[i] ;
		
		
	}
	// Rprintf("Nodes 0 %i, ", (Nodes[0].n)->val) ;
	// Find last time point for each host
	for (i=0 ; i< TN ; i++)
	{
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		if (ha == hb)
		{
			Tend[ha] = i ; 
		}
	}
	
	
	
	
	
	
	for (i=TN-1 ; i>-1 ; i--)
	{
		
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		
		if (!(ha == hb && i>=Tend[ha]))
		{
			if ((int) *PR_T2[3][i] == -1 )
			{
				// death but birth in reverse.
				//Rprintf("%i\n",1) ;
				mn=(NNodes[ha]==0) ? 0 : minNodes[ha] ;
				if (mn<= -1) 
				{
					mn=mn-1 ;
				}
				else
				{
					mn= -1 ;
				}
				
				// add new node
				
				//Nodes[ha].n = Realloc(Nodes[ha].n,NNodes[ha]+1, int) ;
				//Nodes[ha].t = Realloc(Nodes[ha].t,NNodes[ha]+1, double) ;
				//Nodes[ha].n[NNodes[ha]] = mn ;
				//Nodes[ha].t[NNodes[ha]] = *PR_T2[0][i] ;
				
				Nodes[ha].n = llist_add_el_i(mn, Nodes[ha].n) ;
				
				Nodes[ha].t = llist_add_el_d(*PR_T2[0][i], Nodes[ha].t) ;
				NNodes[ha]++ ;
				minNodes[ha] = mn ;
				
			} 
			// need to update minnodes
			
			if ((int) *PR_T2[3][i] == 1 )
			{
				// birth (coalescent event)
				if (ha == hb)
				{
					//Rprintf("%i\n",2) ;
					// sampling without replacement
					r1 = rand() % NNodes[ha] ;
					r2 = rand() % (NNodes[ha] - 1) ;
					r2 = (r2>=r1) ? r2+1 : r2 ;
					//Rprintf("r1,r2 %i,%i",r1,r2) ;
					ch[0] = llist_get_el_i(r1, Nodes[ha].n) ;
					ch[1] = llist_get_el_i(r2, Nodes[hb].n) ;
					//Rprintf("r1,r2 %i,%i",ch[0],ch[1]) ;
					
					
				}
				else{
					//Rprintf("%i\n",3) ;
					r1 = rand() % NNodes[ha] ;
					r2 = rand() % NNodes[hb] ;
					ch[0]=(NNodes[ha]>1) ? llist_get_el_i(r1, Nodes[ha].n) :  llist_get_el_i(0, Nodes[ha].n);
					ch[1]=(NNodes[hb]>1) ? llist_get_el_i(r2, Nodes[hb].n) :  llist_get_el_i(0, Nodes[hb].n);
					
				}
				
				// check for extinct nodes.
				if (ch[1] < 0)
				{
					Nodes[hb].n = llist_delete_el_i(r2, Nodes[hb].n) ;
					Nodes[hb].t = llist_delete_el_d(r2, Nodes[hb].t) ;
					minNodes[hb] = (ch[1]==minNodes[hb]) ? llist_min(Nodes[hb].n) : minNodes[hb] ; 
					NNodes[hb]-- ;
				}
				else{
					if (ch[0] < 0)
					{
						// move node from hb into ha.
						Nodes[ha].n = llist_delete_el_i(r1, Nodes[ha].n) ;
						Nodes[ha].t = llist_delete_el_d(r1, Nodes[ha].t) ;
						r2 = llist_get_ind_i(ch[1], Nodes[hb].n) ;
						Nodes[ha].n = llist_add_el_i(ch[1], Nodes[ha].n) ;
						Nodes[ha].t = llist_add_el_d(llist_get_el_d(r2, Nodes[hb].t), Nodes[ha].t) ;
						
						Nodes[hb].n = llist_delete_el_i(r2, Nodes[hb].n) ;
						Nodes[hb].t = llist_delete_el_d(r2, Nodes[hb].t) ;
						//minNodes[ha] = llist_min(Nodes[ha].n) ;  
						minNodes[hb] = llist_min(Nodes[hb].n) ; 
						NNodes[hb]--;
					}
					else{
						// coalescent event
						// create pair of edges with new internal node.
						tr->edge[ei][0] = intNode ;
						tr->edge[ei][1] = ch[0] ;
						tr->el[ei] = llist_get_el_d(r1, Nodes[ha].t) - *PR_T2[0][i] ;
						//Rprintf("event %i %i %8.4f %8.4f %8.4f\n",intNode, ch[0], llist_get_el_d(r1, Nodes[ha].t), *PR_T2[0][i], tr->el[ei]) ;
						ei--;
						tr->edge[ei][0] = intNode ;
						tr->edge[ei][1] = ch[1] ; 
						tr->el[ei] = llist_get_el_d(r2, Nodes[hb].t) - *PR_T2[0][i] ;
						//Rprintf("event %i %i %8.4f %8.4f %8.4f\n",intNode, ch[1],llist_get_el_d(r2, Nodes[hb].t), *PR_T2[0][i], tr->el[ei]) ;
						ei--;
						tr->nodelabel[intNode-ntips] = ha ;
						
					
						
						Nodes[ha].n = llist_delete_el_i(r1, Nodes[ha].n) ;
						Nodes[ha].t = llist_delete_el_d(r1, Nodes[ha].t) ;
						r2 = llist_get_ind_i(ch[1], Nodes[hb].n) ;
						Nodes[hb].n = llist_delete_el_i(r2, Nodes[hb].n) ;
						Nodes[hb].t = llist_delete_el_d(r2, Nodes[hb].t) ;
						Nodes[ha].n = llist_add_el_i(intNode, Nodes[ha].n) ;
						Nodes[ha].t = llist_add_el_d(*PR_T2[0][i], Nodes[ha].t) ;
						minNodes[hb] = llist_min(Nodes[hb].n) ; 
						minNodes[ha] = llist_min(Nodes[ha].n) ; 
						NNodes[hb]-- ;
						intNode -- ;
					}
				}
			}
		}
	}
	PROTECT(R_List = allocVector(VECSXP ,2)) ;
	PROTECT(R_bt=allocMatrix(REALSXP,tr->NNode,4)) ;
	PR_bt = REAL(R_bt) ;
	bt = phylo_bt(tr, ST, SN, NHosts) ;
	for (j=0 ; j<tr->NNode ; j++)
	{
		for (i=0;i<4;i++)
		{
			PR_bt[j+i*tr->NNode]=(i>0) ? (bt[j][i] + 1) : bt[j][i] ;
			
		}
	}
	R_tr = phylo_to_R(tr) ;
	SET_VECTOR_ELT(R_List, 0, R_tr) ;
	SET_VECTOR_ELT(R_List, 1, R_bt) ;
	UNPROTECT(2) ;
	
	for (j= 0; j< 4 ; j++)
	{
	
			free(PR_T2[j]) ;
		
	}
	
	for (j=0 ; j<tr->NNode ; j++)
	{
		Free(bt[j]) ;
		
	}
	Free(bt) ;
	
	for (j=0 ; j<NHosts ; j++)
	{
		llist_destroy_i(Nodes[j].n)  ;
		llist_destroy_d(Nodes[j].t)  ;
	}
	free(Nodes) ; 
	free(Tend) ;  
	free(NNodes) ;
	free(minNodes) ;
	free(ch) ;
	free(PR_T2) ;
	return R_List;
}

int check_infectives(int *I, int NHosts, int t)
{
	int i=0, ep_end=1;
	while (ep_end==1 && i<NHosts)
	{
		if (I[t*NHosts + i]!=0) ep_end=0;
		i++ ;
	}
	return(ep_end);
}

int check_sampled(double *ST, int NHosts, double t_n, double t_n1 )
{
	// return index of host which was sampled first in the interval (t_n, t_n1)
	int i, rv=-1 ;
	double tmin=t_n1+1;
	
	
	for (i=0 ; i<NHosts ; i++)
	{
		if (ST[i] >= t_n && ST[i]<=t_n1 && ST[i] <= tmin)
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
	while (check_infectives(pp1[1],NHosts,i-1)==0 && i<10000)
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
	int **pp1,i,ot=0;
	double **pp2,tmax=-1 ;
	for (i=0 ; i<NHosts; i++)
	{
	 if (tmax<ST[i]) tmax = ST[i] ;
	}
	
	pp1=Calloc(3,int*) ;
	pp2=Calloc(1,double*) ;
	pp1[0]= Calloc(NHosts,int) ;
	pp1[1]= Calloc(NHosts,int) ;
	pp1[2] = Calloc(NHosts, int) ;
	pp2[0]= Calloc(4,double) ;
	
	
	for (i=0; i<NHosts ; i++)
	{
		pp1[0][i] = nS ;
		pp1[1][i] = 0 ;
		pp1[2][i] = 0 ;
	}
	
	pp1[0][0]-=I0 ;
	pp1[1][0]+=I0 ;
	
	i=1 ;
	n[0]=1 ;
	
	//while (check_infectives(pp1[1],NHosts,n[0]-1)==0 && n[0]<100 && ot==0)
	while (check_infectives(pp1[1],NHosts,n[0]-1)==0 && ot==0)
	{
		cSIR_iter_ST(n, nS, NHosts, B, dr, pp1, pp2, ST, SN) ;
		
		if (pp2[0][(n[0])*4]>tmax) ot=1 ;
		
		//cSIR_iter(i, nS, NHosts, B, dr, pp1, pp2) ;
		//i++ ;
		n[0]++ ;
	}
	p1[0]=pp1[0] ;
	p1[1]=pp1[1] ;
	p1[2]=pp1[2] ;
	p2[0]=pp2[0] ;
	pp1[0]=0;
	pp1[1]=0;
	pp1[2]=0;
	pp2[0]=0;
	n[0]--	;
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
	double *R, sumR=0, t_n, Rtot, *T, rand_e, *ST2, T_n, t_n1;	// rate matrix
	R = Calloc(NHosts * (NHosts+1),double) ;
	ST2 = Calloc(NHosts,double) ;
	int i,j,e,ec,er, *S, *I, *Iend, n, h_i;
	
	S=p1[0];
	I=p1[1];
	Iend=p1[2] ;
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
			//Rprintf("%8.4f\t",R[i*NHosts+j]) ;
		}
		//Rprintf("\n") ;
	}
	for (i=0; i<NHosts; i++)
	{
		//Rprintf("%8.4f\t",R[i+NHosts*NHosts]) ;
	}
	//Rprintf("\n\n\n") ;
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
			Iend[h_i] = I[h_i+((n-1)*NHosts)]; 
			T[n*4+1] = h_i ;
			T[n*4+2] = h_i ;	
			T[n*4] = ST[h_i] ;		
					
			
			
			n++ ;
			T_n=T[(n-1)*4] ;
			t_n=t_n1-T_n ;
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
		R[i]=(sumR==0) ? 1.0/(double)(NHosts*(NHosts+1)) : R[i]/sumR ;
	}
	
	// if there are any more infected left ..
	
	// increase size of S,I,T
		S=Realloc(S,(n+1)*NHosts ,int) ;
		I=Realloc(I,(n+1)*NHosts ,int) ;
		T=Realloc(T,(n+1)*4 ,double) ;
	for (i=0 ; i< NHosts; i++)
	{
			S[i+(n*NHosts)] = S[i+((n-1)*NHosts)] ;
			I[i+(n*NHosts)] = I[i+((n-1)*NHosts)] ;		
	}	
		
	if (check_infectives(&I[0], NHosts, n-1)!=1) 
	{
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
	
		
		
	/*
		for (i=0 ; i< NHosts; i++)
		{
			S[i+(n*NHosts)] = S[i+((n-1)*NHosts)] ;
			I[i+(n*NHosts)] = I[i+((n-1)*NHosts)] ;
		
		}
	*/
	
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
		
	}
	T[n*4] = t_n+T[(n-1)*4];
	
	p1[0]=S ;
	p1[1]=I ;
	p1[2]=Iend ;
	p2[0]=T ;
	S=0;
	I=0;
	T=0;
	Iend=0;
	np[0] = n ;
	Free(R) ;
	Free(ST2) ;
}
