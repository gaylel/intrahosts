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
SEXP sample_cSIR_S_R(SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN, SEXP R_bnprob, SEXP R_K) ;
void cSIR_iter(int n, int nS, int NHosts, double *B, double dr, int **p1, double **p2) ;
void cSIR_iter_ST(int *np, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN, int * bnsizes) ;
void cSIR_iters(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2) ;
void cSIR_iters_ST(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN, double bnprob, int K) ;
int check_infectives(int *I, int NHosts, int t) ;
int check_sampled(double *ST, int NHosts, double t_n, double t_n1 ) ;
SEXP tree_reconstruct(SEXP R_sir, SEXP R_NHosts, SEXP R_dat) ;
SEXP tree_reconstruct_withll(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_ktt) ;
SEXP tree_reconstruct_old(SEXP R_sir, SEXP R_NHosts, SEXP R_dat) ;
void bnsizeupdate(int * bnsizes, int NHosts, int opt, int * opt_val, double opt_val2, int nS, int K) ;
SEXP tree_reconstruct_with_partialll(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS, SEXP R_OBS2, SEXP R_tr_old) ;
SEXP tree_reconstruct_with_partialll2(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS, SEXP R_OBS2, SEXP R_tr_old, SEXP R_mig2) ;
double ** calculateSIR_rates(double ** lambdao, double *B, int *PR_I, int *PR_S, int ha, int hb, int i, int kt, int NHosts, int NS, int TN, int *Anc, int *OBS2) ;
double calculatelambdaosum(double ** lambdao, int NHosts) ;


struct hnode {
	item_i *n ;
	item_d *t ;
} ;

typedef struct hnode hnode ;

struct hnode2 {
	item_i *n ;
	item_d *t ;
	int minNode ;
	int N ;
} ;

typedef struct hnode2 hnode2 ;

hnode2  Nodes_addelement(int val, double t, hnode2 Nodes) ;
hnode2  Nodes_initialise(hnode2 Nodes) ;
int * Nodes_pick_pair(int *ch, hnode2 Nodesa, hnode2 Nodesb, int ha, int hb) ;
hnode2 Nodes_deleteelement(int val, hnode2 Nodes) ;
int * Nodes_pick_sets(int *ov, hnode2 Nodesa1, hnode2 Nodesa2, hnode2 Nodesb1, hnode2 Nodesb2, int ha, int hb) ;
int * Nodes_pick_children(int *ch, int *ov, hnode2 Nodesa1, hnode2 Nodesa2, hnode2 Nodesb1, hnode2 Nodesb2, int ha, int hb) ;


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


SEXP sample_cSIR_S_R(SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN, SEXP R_bnprob, SEXP R_K)
{
	
	int I0, NS, NHosts, *S, *I, *T, *R, **p_RVAL1 = Calloc(3,int*),*n=Calloc(1,int),i,j;
	int *PR_I, *PR_S, *PR_Iend, *SN, K;
	double *B, **p_RVAL2 = Calloc(1,double*), dr, *PR_T, *ST, bnprob ;
	
	
	SEXP R_Bdim, R_I, R_T, R_S, R_Iend, R_list ;
	
	// initial number of infective particles
	R_I0 = coerceVector(R_I0, INTSXP) ;
	
	// number of particles in each host	
	R_NS = coerceVector(R_NS, INTSXP) ;
  
  // parameters for bottleneck size
  	R_bnprob = coerceVector(R_bnprob, REALSXP) ;
	
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
	
	K = INTEGER(coerceVector(R_K, INTSXP))[0] ;
	I0 = INTEGER(R_I0)[0] ;
	NS = INTEGER(R_NS)[0] ;
  bnprob = REAL(R_bnprob)[0] ;
	NHosts = INTEGER(R_NHosts)[0] ;
	B = REAL(R_B) ;
	dr = REAL(R_dr)[0] ;
	ST = REAL(R_ST) ;
	SN = INTEGER(R_SN) ;


	R_Bdim = getAttrib(R_B, R_DimSymbol) ;
	n[0]=0 ;
	
	//cSIR_iters(n, I0, NS, NHosts, B, dr, p_RVAL1, p_RVAL2) ;
	cSIR_iters_ST(n, I0, NS, NHosts, B, dr, p_RVAL1, p_RVAL2, ST, SN, bnprob, K) ;
	
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

SEXP tree_reconstruct_old(SEXP R_sir, SEXP R_NHosts, SEXP R_dat)
{
	int *PR_I, *PR_S, *PR_Iend, *SN, NHosts, i, j, k,st=0, st2=-1, TN, ha, hb, *Tend, *NNodes, *minNodes, mn, r1, r2, *ch, ntips=0, ei, intNode;
  int bn ;
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
		Rprintf("%8.4f\t%i\t%i\t", ha+1, hb+1, *PR_T2[0][i] ) ;
		if (!(ha == hb && i>=Tend[ha]))
		{
			if ((int) *PR_T2[3][i] == -1 )
			{
				// death but birth in reverse.
				Rprintf("death\n") ;
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
			
			if ((int) *PR_T2[3][i] > 0 )
			{
				Rprintf("birth\t") ;
        		for (bn=0 ; bn < (int)*PR_T2[3][i] ; bn++)
        		{
          			// birth (coalescent event)
          			// choose  children
  			  		if (ha == hb)
					{
					  Rprintf("within,") ;
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
						Rprintf("between,") ;
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % NNodes[hb] ;
					  	ch[0]=(NNodes[ha]>1) ? llist_get_el_i(r1, Nodes[ha].n) :  llist_get_el_i(0, Nodes[ha].n);
					  	ch[1]=(NNodes[hb]>1) ? llist_get_el_i(r2, Nodes[hb].n) :  llist_get_el_i(0, Nodes[hb].n);
				  	}
				  	// check for extinct nodes.
				  	if (ch[1] < 0) // -
				  	{
				  		
					  	Nodes[hb].n = llist_delete_el_i(r2, Nodes[hb].n) ;
					  	Nodes[hb].t = llist_delete_el_d(r2, Nodes[hb].t) ;
					  	minNodes[hb] = (ch[1]==minNodes[hb]) ? llist_min(Nodes[hb].n) : minNodes[hb] ; 
					  	NNodes[hb]-- ;
				  	}
				  	else{
					  	if (ch[0] < 0) // +-
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
					  	else{ //++
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


SEXP tree_reconstruct_with_partialll(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS, SEXP R_OBS2, SEXP R_tr_old)
{
   int  st=0, st2=-1, *minNodes, mn, r1, r2, *ch, ei, intNode;
  	int bn , m, ll_calc=1, hi, hj, Ii, Sj, *Anc, kt, ktt, Ancsum=0, ce = 0, ce1= 0;
  	double *PR_bt, **bt, v1, v2, *vvec, ll=0.0, lambdao_sum=0.0, *PR_ll, ll_i ;

	item_i * tst ;
	
	SEXP R_tr, R_T, R_List, R_bt, R_ll ;
	
	
	/******* Interface with R to construct tree from trajectory **************************
	INPUT:	R_sir		list containing I (Infected) S (Susceptible) T (times of events)
			R_NHosts 	Number of hosts
			R_dat		list containing SN (Number of sequences for each host) ST (time of sampling for each host)
			R_B			matrix of rates
			R_NS		number in SIR population for each host
			R_BN		
			R_OBS		old number of observed
			R_OBS2		additional number of observed
			R_tr_old
	
	
	*************************************************************************************/
	
	// from R_sir
	SEXP R_Tdim = getAttrib(VECTOR_ELT(R_sir,2), R_DimSymbol) ;
	int TN = INTEGER(R_Tdim)[0] ;		// number of time points of trajectory
	int *PR_I = INTEGER(VECTOR_ELT(R_sir,0)) ;
	int *PR_S = INTEGER(VECTOR_ELT(R_sir,1)) ;	
	int *PR_Iend = INTEGER(VECTOR_ELT(R_sir,3)) ;
	double *PR_T = REAL(coerceVector(VECTOR_ELT(R_sir,2),REALSXP)) ;	
	
	//from R_NHosts
	int NHosts = INTEGER(coerceVector(R_NHosts, INTSXP))[0] ;
	
	//from R_dat
	int *SN = INTEGER(coerceVector(VECTOR_ELT(R_dat,0), INTSXP)) ;
	double *ST = REAL(coerceVector(VECTOR_ELT(R_dat,1), REALSXP)) ;
	
	double *B = REAL(coerceVector(R_B, REALSXP)) ;
	int NS = INTEGER(coerceVector(R_NS, INTSXP))[0] ;
	double BN = REAL(coerceVector(R_BN, REALSXP))[0] ;
	int *OBS = INTEGER(coerceVector(R_OBS, INTSXP)) ;
	int *OBS2 = INTEGER(coerceVector(R_OBS2, INTSXP)) ;
	int sumOBS2=0, i1 ;
	phylo *tr_old ;
	for (i1=0 ; i1 <NHosts ; i1++)
	{
		sumOBS2 += OBS2[i1] ;
	}
	Rprintf("%i\n", sumOBS2) ;
	if (sumOBS2 > 0)
	{
		tr_old = R_to_phylo(R_tr_old) ;
		// calculate branching times
		double **bt_old = phylo_bt(tr_old, ST, OBS2, NHosts) ;
	}
	
	
	/*************** INTERNAL VARIABLES FOR ALGORITHM ***********************************/
	
	int i, j, k=0, ha, hb ;
	double ***PR_T2 ;
	PR_T2 = Calloc(4,double**) ;	//reshape(PR_T)
	for (j=0 ; j<4 ; j++)
	{
		PR_T2[j] = Calloc(TN,double*) ;
		for (i=0 ; i<TN ; i++)
		{	
		 	PR_T2[j][i] = &PR_T[k++];
		}
	}
	int *Tend = Calloc(NHosts,int) ; // times at which hosts are sampled
	for (i=0 ; i< TN ; i++)
	{
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		if (ha == hb && (*PR_T2[0][i] <= ST[ha]))
		{
			Tend[ha] = i ; 
		}
	}
	
	int *SNsum = calloc(NHosts, sizeof(int)) ;
	int ntips = 0;
	double **lambdao = calloc(NHosts, sizeof(double*)) ;
	for (i=0 ; i<NHosts ; i++)
	{
		ntips+= OBS[i];
		lambdao[i] = calloc(NHosts, sizeof(double)) ;
	}
	hnode * Nodes = Calloc(NHosts,hnode) ; 
	hnode * Nodes2 = Calloc(NHosts, hnode) ;
	phylo *tr = phylo_create(ntips) ;
	
	
	int * NNodes = Calloc(NHosts,int) ;
	int * NNodes2 = Calloc(NHosts, int) ;
	minNodes = Calloc(NHosts,int) ;
	ch = Calloc(2,int) ;
	
	ei = tr->Nedge - 1 ;
	intNode = 2*tr->NNode  ;
	Anc = calloc(NHosts, sizeof(int)) ;
	vvec = calloc(NS, sizeof(double)) ;
	
	
	
	
  	j = 0;
  	// maybe still some infected left at end:
  
	for (i=0 ; i<NHosts ; i++)
	{
		Nodes[i].t = NULL ;
		Nodes[i].n = NULL ;
   		NNodes[i] = 0 ;
    	minNodes[i] = 0 ;
    	//SNsum[i] = SN[i] + j ;
    	//j+= SN[i] ;
    	SNsum[i] = OBS[i] + j ;
    	j+= OBS[i] ;
    	st2 = -1 ;
    	//Rprintf("Iend in host[%i] = %i\n", i, PR_I[TN * (i+1) - 1]) ;
    	for (k=0 ; k< PR_I[TN * (i+1) - 1] ; k++)
    	{
  	    	Nodes[i].n = llist_add_el_i(st2, Nodes[i].n) ;
			Nodes[i].t = llist_add_el_d(ST[NHosts - 1], Nodes[i].t) ;
			//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
        	NNodes[i]++ ;
			st2-- ;
		  	minNodes[i] -- ;         
		}
		Anc[i] = 0 ;
	}


/*******************************ALGORITHM************************************************/

	for (i=TN-1 ; i>-1 ; i--)
	{
		ll_i = 0 ;
		ce = 0 ;
		ce1 = 0 ;
    // for each event get hosts
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		
   		
    // if terminal event
   		if (i==Tend[ha])
    	{
      		st = (ha>0 ? SNsum[ha-1] : 0) ;
      		//NNodes[ha] += -*PR_T2[3][i] ;
      		st2 = minNodes[ha] - 1;
		  	//Rprintf("Host sampled\n") ;
		  	
		  	
      		for (j=0 ; j<OBS[ha] ; j++)
		  	{
				Nodes[ha].n = llist_add_el_i(st, Nodes[ha].n) ;
			  	Nodes[ha].t = llist_add_el_d(ST[ha], Nodes[ha].t) ;
			  	NNodes[ha]++ ;
       		  	st++ ;	
		  	}
		  	for (j=OBS[ha] ; j< -*PR_T2[3][i] ; j++)
		  	{
				//Rprintf("%i, ",st2) ;
        	  	Nodes[ha].n = llist_add_el_i(st2, Nodes[ha].n) ;
			  	Nodes[ha].t = llist_add_el_d(ST[ha], Nodes[ha].t) ;
			  	//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
              	NNodes[ha]++ ;
			  	st2=st2-1 ;
		  	}
		  	if (OBS[ha]< -*PR_T2[3][i])
		  	{
				minNodes[ha] = st2+1 ;  
		  	}
		  	Anc[ha]+=OBS[ha] ;   
      	}
    	else
		{
			// calculate rates for each host
			Ancsum = 0 ;
			for (kt = 0 ; kt < NHosts ; kt++)
			{
				Ancsum += Anc[kt] ; 
			}
			
			//if (i > 0 & Ancsum > 1)
			if (i > 0 & (intNode >= ntips))
			{
			lambdao_sum = 0 ;
			kt = *PR_T2[3][i] ;	
			// hi == hj
			for (hj=0 ; hj < NHosts ; hj++)
			{
				Ii = PR_I[hj*TN + i] ;
				Sj = PR_S[hj*TN + i] ;
				if (hj == hb)
				{
					Ii = Ii - kt ;
					Sj = (kt > 0) ? Sj + kt : Sj;
				}
				lambdao[hj][hj] = B[(NHosts * hj) + hj]  *(double) (Ii * Sj) / NS ; 
				lambdao[hj][hj] = (Anc[hj] >0) ? choose((double)(Anc[hj]), (double)(2)) / choose((double)(PR_I[hj*TN + i]), (double)(2)) *lambdao[hj][hj] : 0;
				lambdao_sum += lambdao[hj][hj] ;
			}
			//Rprintf("here\n") ;
			// hi != hj
			for (hj=0 ; hj < NHosts ; hj++)
			{
				Sj = PR_S[hj*TN + i] ;
				ktt = 100 ;
				if (PR_I[hj * TN + i] >= ktt)
				{
					
					vvec[0] = 1 - choose((double)(PR_I[hj*TN + i] - Anc[hj]), (double)(ktt)) / choose((double)(PR_I[hj*TN + i]), (double)(ktt)) ; 
				//	vvec[0] = (double) ktt / PR_I[hj * TN + i] ;
				}
				else 
				vvec[0] = 0;
					
				for (hi=0 ; hi < NHosts ; hi++)
				{
				
					if (hi != hj)
					{
						Ii = PR_I[hi*TN + i] ;
						if (hj == hb)
						{
							
							Sj = (kt > 0) ? Sj + kt : Sj;
						}
					    //vvec[0] = (Ii >= ktt) ? vvec[0] : 0 ;
						lambdao[hj][hi] = B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / (NS * ktt) ; 
						lambdao[hj][hi] = lambdao[hj][hi] * ((Ii >= ktt) ? vvec[0] : 0) ;
						lambdao_sum += lambdao[hj][hi] ;
					}
					else
					{
						//vvec[0] = 0 ;
					}
		//			Rprintf("%8.4f %8.4f ", vvec[0], B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / NS) ;
					//Rprintf("%8.4f ", vvec[0], B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / NS) ;
				
				}
		//		Rprintf("\n") ;
				
			}
			/*
			for (hj = 0 ; hj<NHosts ; hj++)
			{
				Rprintf(" %i ", PR_I[hj*TN + i]) ;
			}
			Rprintf("\n") ;
			*/
			
			//Rprintf("%8.4f\n", 1/lambdao_sum) ;
			//Rprintf("%8.4f\n", ll) ;
			
			//Rprintf("%i\n", kt) ;
		
	/*	
			for (hj=0 ; hj<NHosts ; hj++)
			{
				for (hi=0 ; hi<NHosts ; hi++)
				{
					Rprintf("%8.4f ", lambdao[hj][hi]) ;
				}
				Rprintf("\n") ;
			}
			
			Rprintf("\n") ;  */
			ll_i = -lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
			//ll += ll_i ; //-lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
			
		//Rprintf("here\n") ;
    		if ((int) *PR_T2[3][i] == -1 )
			{
				// death but birth in reverse.
				//Rprintf("death\n") ;
				// can calculate the probability of no coalescence has occurred
				
				// minnode
        		mn = (NNodes[ha]==0) ? 0 : minNodes[ha] ;
				if (mn<= -1) 
				{
					mn=mn-1 ;
				}
				else
				{
					mn= -1 ;
				}
				// add new node
				Nodes[ha].n = llist_add_el_i(mn, Nodes[ha].n) ;
				Nodes[ha].t = llist_add_el_d(*PR_T2[0][i], Nodes[ha].t) ;
				NNodes[ha]++ ;
				minNodes[ha] = mn ;
			} 
			if ((int) *PR_T2[3][i] > 0 )
			{
        		//Rprintf("birth\n") ;
        		for (bn=0 ; bn < (int)*PR_T2[3][i] ; bn++)
        		{
        		//Rprintf("here %i of %i\n", bn, (int)*PR_T2[3][i]) ;
          			// birth (coalescent event)
          			// choose  children
  			  		if (ha == hb)
				  	{
            			//Rprintf("got here\n" ) ;  
					  	//Rprintf("\twithin, ") ;
					  	// sampling without replacement
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % (NNodes[ha] - 1) ;
					  	r2 = (r2>=r1) ? r2+1 : r2 ;
					 	// Rprintf("r1,r2 %i,%i",r1,r2) ;
					  	ch[0] = llist_get_el_i(r1, Nodes[ha].n) ;
					  	ch[1] = llist_get_el_i(r2, Nodes[hb].n) ;
					  	//Rprintf("r1,r2 %i,%i",ch[0],ch[1]) ;
        			}
				  	else{
					  	//Rprintf("\tbetween, ") ;
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % NNodes[hb] ;
					  	ch[0]=(NNodes[ha]>1) ? llist_get_el_i(r1, Nodes[ha].n) :  llist_get_el_i(0, Nodes[ha].n);
					  	ch[1]=(NNodes[hb]>1) ? llist_get_el_i(r2, Nodes[hb].n) :  llist_get_el_i(0, Nodes[hb].n);
					  	//Rprintf("%i %i %i %i\n", ch[0], ch[1], NNodes[ha], NNodes[hb]) ;
				  	}
          
				  	// check for extinct nodes.
				  	if (ch[1] < 0) // -
				  	{
				  		//Rprintf("one extinct 1\n") ;
					  	Nodes[hb].n = llist_delete_el_i(r2, Nodes[hb].n) ;
					  	Nodes[hb].t = llist_delete_el_d(r2, Nodes[hb].t) ;
					  	minNodes[hb] = (ch[1]==minNodes[hb]) ? llist_min(Nodes[hb].n) : minNodes[hb] ; 
					  	NNodes[hb]-- ;
				  	}
				  	else{
					  	if (ch[0] < 0) // +-
					  	{
					  		//Rprintf("one extinct 2\n") ; 
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
						  	if (ha != hb) ce1++ ;
						  	//if (ha != hb) ll += log(lambdao[hb][ha]) ;
						  	Anc[hb]-- ;
						  	Anc[ha]++ ;
					  	}
					  	else{ //++
							// coalescent event
							// create pair of edges with new internal node.
							tr->edge[ei][0] = intNode ;
							tr->edge[ei][1] = ch[0] ;
							tr->el[ei] = llist_get_el_d(r1, Nodes[ha].t) - *PR_T2[0][i] ;
							//Rprintf("event %i %i %8.4f %8.4f %8.4f\n",intNode, ch[0], llist_get_el_d(r1, Nodes[ha].t), *PR_T2[0][i], tr->el[ei]) ;
							//Rprintf("event %i %i %i %8.4f\n",ha+1, ha+1, hb+1, *PR_T2[0][i]) ;
						
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
							//ll += log(lambdao[hb][ha]) ;
							intNode -- ;
							Anc[hb] -- ;
							ce++ ;
						}
				  }  
        	}		
			}
			if (ce > 0 || ce1 > 0 )
			{
				Ii = PR_I[hb*TN + i] ;
				//Ii = Anc[hb]  + 1 ;
				//if (ha == hb)
				//	Ii = Ii - kt ;
				//Sj = (kt > 0) ? Sj + kt : Sj;
						
				//ll_i+= log(B[(NHosts * hb) + ha]  *(double) (Ii * Sj) / NS ); 
			if (ha == hb)
			{
				ll_i+= log(lambdao[hb][ha]) - log((Anc[hb]+1)*(Anc[hb])/2) ;
				//Rprintf("%8.4f\n", ll_i) ;
			}else{
				ll_i+= log(lambdao[hb][ha]) ;
				//Rprintf("%8.4f\n", ll_i) ;
				ktt = *PR_T2[3][i] ;
				Ii = PR_I[ha*TN + i] ;
				
				ll_i+= log(choose(Ii-(Anc[ha]-ce1), ce1 + ce - ce)*gamma(ce1 + 1))  ;
				ll_i -= log( choose(Anc[hb] + ce1 + ce, ce1 + ce) ) ;
				ll_i -= log( choose(Ii, Ii - (ce1 + ce)) * gamma(ce1 + ce + 1) ) ;
	//			Rprintf("Ii=%i, Ij=%i, ni=%i, m=%i, m'=%i, nj=%i\n", Ii, PR_I[hb*TN + i], Anc[ha]-ce1, ce1+ ce, ce, Anc[hb]+ce1+ce) ;
			 	//ll_i-= (log(choose((double)(Ii), (double)(ktt))  -  choose((double)(Ii - Anc[hb] - 1), (double)(ktt))) ) ;
			 	//ll_i -= (log(choose((double)(PR_I[ha*TN + i]), (double)(ktt)) * gamma(ktt + 1)))   ;
			
			//ll_i+= dhyper(1, Anc[hb]+1, Ii - Anc[hb] - 1, ktt, 1) ;
			
				//ll_i-=  log(choose(Ii, *PR_T2[3][i])-choose(Ii - (Anc[hb] + 1), *PR_T2[3][i])) ; 
				//ll_i+= log(*PR_T2[3][i] / Ii) ;
				//ll_i-= log(PR_I[ha*TN + i]) ;
			//	Rprintf("%8.4f\t%i\t%i\t%8.4f\n", ll_i, Ii, Anc[hb], *PR_T2[3][i]) ;
				//ll_i-=  log(choose(PR_I[ha*TN + i], PR_I[ha*TN + i] - *PR_T2[3][i] )) ; 
				//Rprintf("%8.4f\n", ll_i) ;
			}
			//ll_i = log(1 - exp(ll_i)) + log(lambdao[hb][ha]) - log(lambdao_sum);
			}
			//printf("ll_i = %8.4f\n", ll_i) ; 
			//if (Ancsum > 1) ll+= ll_i ;
			ll += ll_i ;
			}
			else
			{
				//Rprintf("end \n") ;
			}
			
			
		//	Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[ha], Anc[hb]) ;
			
		
		}
		//Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[0], Anc[1], Anc[2], Anc[3], Anc[4], NNodes[0], NNodes[1], Ancsum, intNode) ;
			
		
}


/*******************************RESULT TO OUTPUT*****************************************/

	PROTECT(R_List = allocVector(VECSXP ,3)) ;
	PROTECT(R_bt=allocMatrix(REALSXP,tr->NNode,4)) ;
	PROTECT(R_ll = allocVector(REALSXP, 1)) ;
	PR_bt = REAL(R_bt) ;
	bt = phylo_bt(tr, ST, OBS, NHosts) ;
	for (j=0 ; j<tr->NNode ; j++)
	{
		for (i=0;i<4;i++)
		{
			PR_bt[j+i*tr->NNode]=(i>0) ? (bt[j][i] + 1) : bt[j][i] ;
			
		}
	}
	R_tr = phylo_to_R(tr) ;
	PR_ll = REAL(R_ll) ;
	PR_ll[0] = ll ;
	
	SET_VECTOR_ELT(R_List, 0, R_tr) ;
	SET_VECTOR_ELT(R_List, 1, R_bt) ;
	SET_VECTOR_ELT(R_List, 2, R_ll) ;
	UNPROTECT(3) ;
	
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
		free(lambdao[j]) ;
	}
	free(Nodes) ; 
	free(Tend) ;  
	free(NNodes) ;
	free(minNodes) ;
	free(ch) ;
	free(PR_T2) ;
 	free(SNsum) ;
 	free(Anc) ;
 	free(lambdao) ;
 	free(vvec) ;
	return R_List;
}

double ** calculateSIR_rates(double ** lambdao, double *B, int *PR_I, int *PR_S, int ha, int hb, int i, int kt, int NHosts, int NS, int TN, int *Anc, int *OBS2)
{
	/******************************* INPUT ********************************************
	
	
	
	
	**********************************************************************************/

	int hi, hj ;
	int Ii ;
	int Sj ;
	int ktt=100 ;
	double vvec ;
	
	for (hj=0 ; hj < NHosts ; hj++)
	{
		Ii = PR_I[hj*TN + i] ;
		Sj = PR_S[hj*TN + i] ;
		if (hj == hb)
		{
			Ii = Ii - kt ;
			Sj = (kt > 0) ? Sj + kt : Sj;
		}
		lambdao[hj][hj] = B[(NHosts * hj) + hj]  *(double) (Ii * Sj) / NS ; 
		lambdao[hj][hj] = (Anc[hj] >0) ? (choose((double)(Anc[hj]), (double)(2))  - choose((double)(OBS2[hj]), (double)(2)) ) / ((PR_I[hj*TN + i] * (PR_I[hj*TN + i] - 1)) /2) *lambdao[hj][hj] : 0;
	}
	for (hj=0 ; hj < NHosts ; hj++)
	{
		Sj = PR_S[hj*TN + i] ;
		if (PR_I[hj * TN + i] >= ktt)
		{				
			//vvec[0] = 1 - choose((double)(PR_I[hj*TN + i] - Anc[hj]), (double)(ktt)) / choose((double)(PR_I[hj*TN + i]), (double)(ktt)) ; 
			vvec = (choose((double)(PR_I[hj*TN + i] - OBS2[hj]), (double)(ktt)) - choose((double)(PR_I[hj*TN + i] - Anc[hj]), (double)(ktt)))/ choose((double)(PR_I[hj*TN + i]), (double)(ktt)) ; 
		}
		else 
		{
			vvec = 0;
		}		
		for (hi=0 ; hi < NHosts ; hi++)
		{
			if (hi != hj)
			{
				Ii = PR_I[hi*TN + i] ;
				if (hj == hb)
				{				
					Sj = (kt > 0) ? Sj + kt : Sj;
				}
			    //vvec[0] = (Ii >= ktt) ? vvec[0] : 0 ;
				lambdao[hj][hi] = B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / (NS * ktt) ; 
				lambdao[hj][hi] = lambdao[hj][hi] * ((Ii >= ktt) ? vvec : 0) ;
				
			}
		}
	}
	return lambdao ;
}

double  calculatelambdaosum(double ** lambdao, int NHosts)
{
	double lambdaosum=0 ;
	int i, j ;
	for (i = 0 ; i<NHosts ; i++)
	{
		for (j = 0 ; j<NHosts ; j++)
		{
			lambdaosum += lambdao[i][j] ;
		}
	}
	return lambdaosum ;
}

hnode2 Nodes_addelement(int val, double t, hnode2 Nodes)
{
	Nodes.n = llist_add_el_i(val, Nodes.n) ;
   	Nodes.t = llist_add_el_d(t, Nodes.t) ;
	Nodes.N++ ;
	if (val < Nodes.minNode)
	{		
		Nodes.minNode = val ;
	}
	return Nodes ; 
} 

hnode2 Nodes_deleteelement(int val, hnode2 Nodes)
{
	int val_i ;
	val_i = llist_get_ind_i(val, Nodes.n) ;
	Nodes.n = llist_delete_el_i(val_i, Nodes.n) ;
	Nodes.t = llist_delete_el_d(val_i, Nodes.t) ;
	Nodes.N-- ;
	Nodes.minNode = (val == Nodes.minNode) ? llist_min(Nodes.n) : Nodes.minNode ; 
	Nodes.minNode = (Nodes.N == 0) ? 0 : Nodes.minNode ;
	return Nodes ;				  	
}
 
 
hnode2 Nodes_initialise(hnode2 Nodes)
{
	Nodes.n = NULL ;
	Nodes.t = NULL ;
	Nodes.N = 0 ;
	Nodes.minNode = 0 ;
}

int * Nodes_pick_pair(int *ch, hnode2 Nodesa, hnode2 Nodesb, int ha, int hb)
{
	// sampling without replacement
	int r1 ;
	int r2 ;
	r1 = (Nodesa.N > 1) ? rand() % Nodesa.N : 0;
	ch[0] = llist_get_el_i(r1, Nodesa.n) ;
	if (ha != hb)
	{
		r2 = (Nodesb.N > 1) ? rand() % Nodesb.N : 0 ;
	}
	else
	{
		r2 = (Nodesb.N - 1 > 1) ? rand() % (Nodesb.N - 1) : 0 ;
		r2 = (r2>=r1) ? r2+1 : r2 ;
	}
	ch[1] = llist_get_el_i(r2, Nodesb.n) ;
	ch[2] = r1 ;
	ch[3] = r2 ;
	return ch ;				  	
}

int * Nodes_pick_sets(int *ov, hnode2 Nodesa1, hnode2 Nodesa2, hnode2 Nodesb1, hnode2 Nodesb2, int ha, int hb)
{
	ov[0] = 0;
	ov[1] = 0;
	if (rbinom(1, (double) Nodesa2.N/(Nodesa2.N + Nodesa1.N)) == 1)
	{
		ov[0] = 1 ; // 1,0
	}
	else{
		if (rbinom(1, (double) Nodesb2.N/(Nodesb2.N + Nodesb1.N)) == 1)
		{
			ov[1] = 1 ; // 0,1
		}
	}
	// 0,0
	return ov ;
					  	
}

int * Nodes_pick_children(int *ch, int *ov, hnode2 Nodesa1, hnode2 Nodesa2, hnode2 Nodesb1, hnode2 Nodesb2, int ha, int hb)
{
// sampling without replacement
	if (ov[0] == 0)        		  	
	{		  		
        if (ov[1] == 0)    	
        {
        	ch = Nodes_pick_pair(ch, Nodesa1, Nodesb1, ha, hb) ;	  		
        }
        else
 		{            		  		
        	ch = Nodes_pick_pair(ch, Nodesa1, Nodesb2, ha, hb) ;
        
        }
    }else
        {
            ch = Nodes_pick_pair(ch, Nodesa2, Nodesb1, ha, hb) ;
    }
    
    return ch ;
}

SEXP tree_reconstruct_with_partialll2(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS2, SEXP R_OBS, SEXP R_tr_old, SEXP R_mig2)
{
   int  st=0, st2=-1, *minNodes, mn, r1, r2, *ch, ei, intNode;
  	int bn , m, ll_calc=1, hi, hj, Ii, Sj, *Anc, kt, ktt, Ancsum=0, ce = 0, ce1= 0;
  	double *PR_bt, **bt, v1, v2, *vvec, ll=0.0, lambdao_sum=0.0, *PR_ll, ll_i ;

	item_i * tst ;
	
	SEXP R_tr, R_T, R_List, R_bt, R_ll, R_mig ;
	
	
	/******* Interface with R to construct tree from trajectory **************************
	INPUT:	R_sir		list containing I (Infected) S (Susceptible) T (times of events)
			R_NHosts 	Number of hosts
			R_dat		list containing SN (Number of sequences for each host) ST (time of sampling for each host)
			R_B			matrix of rates
			R_NS		number in SIR population for each host
			R_BN		
			R_OBS2		old number of observed
			R_OBS		additional number of observed
			R_tr_old
			R_mig2		old migration matrix
	
	
	*************************************************************************************/
	
	// from R_sir
	SEXP R_Tdim = getAttrib(VECTOR_ELT(R_sir,2), R_DimSymbol) ;
	int TN = INTEGER(R_Tdim)[0] ;		// number of time points of trajectory
	int *PR_I = INTEGER(VECTOR_ELT(R_sir,0)) ;
	int *PR_S = INTEGER(VECTOR_ELT(R_sir,1)) ;	
	int *PR_Iend = INTEGER(VECTOR_ELT(R_sir,3)) ;
	double *PR_T = REAL(coerceVector(VECTOR_ELT(R_sir,2),REALSXP)) ;	
	
	//from R_NHosts
	int NHosts = INTEGER(coerceVector(R_NHosts, INTSXP))[0] ;
	
	//from R_dat
	int *SN = INTEGER(coerceVector(VECTOR_ELT(R_dat,0), INTSXP)) ;
	double *ST = REAL(coerceVector(VECTOR_ELT(R_dat,1), REALSXP)) ;
	
	double *B = REAL(coerceVector(R_B, REALSXP)) ;
	int NS = INTEGER(coerceVector(R_NS, INTSXP))[0] ;
	double BN = REAL(coerceVector(R_BN, REALSXP))[0] ;
	int *OBS = INTEGER(coerceVector(R_OBS, INTSXP)) ;
	int *OBS2 = INTEGER(coerceVector(R_OBS2, INTSXP)) ;
	int sumOBS=0, i1 ;
	phylo *tr_old ;
	for (i1=0 ; i1 <NHosts ; i1++)
	{
		sumOBS += OBS2[i1] ;
	}
	Rprintf("%i\n", sumOBS) ;
	double **bt_old ;
	double *mig2 ;
	double ttmp ;
	int mig2_n ;
	int bt_n ;
	if (sumOBS > 0)
	{
		tr_old = R_to_phylo(R_tr_old) ;
		// calculate branching times
		bt_old = phylo_bt(tr_old, ST, OBS, NHosts) ;
		bt_n = tr_old->NNode - 1 ;
		mig2 = REAL(coerceVector(R_mig2, REALSXP)) ; 
		mig2_n = INTEGER(getAttrib(R_mig2, R_DimSymbol))[0] ;
		Rprintf("here\n") ;
	}
	
	
	/*************** INTERNAL VARIABLES FOR ALGORITHM ***********************************/
	
	int i, j, k=0, ha, hb, oa, ob, Ntrans ;
	double ***PR_T2 ;
	PR_T2 = Calloc(4,double**) ;	//reshape(PR_T)
	for (j=0 ; j<4 ; j++)
	{
		PR_T2[j] = Calloc(TN,double*) ;
		for (i=0 ; i<TN ; i++)
		{	
		 	PR_T2[j][i] = &PR_T[k++];
		}
	}
	int *Tend = Calloc(NHosts,int) ; // times at which hosts are sampled
	for (i=0 ; i< TN ; i++)
	{
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		if (ha == hb && (*PR_T2[0][i] <= ST[ha]))
		{
			Tend[ha] = i ; 
		}
	}
	
	int *SNsum = calloc(NHosts, sizeof(int)) ;
	int ntips = 0;
	double **lambdao = calloc(NHosts, sizeof(double*)) ;
	for (i=0 ; i<NHosts ; i++)
	{
		ntips+= OBS2[i] + OBS[i];
		lambdao[i] = calloc(NHosts, sizeof(double)) ;
	}
	hnode2 * Nodes = Calloc(NHosts,hnode2) ; 
	hnode2 * Nodes2 = Calloc(NHosts, hnode2) ;
	hnode2 hnode2v ;
	
	for (i=0 ; i<NHosts ; i++)
	{
		Nodes[i] = Nodes_initialise(Nodes[i]) ;
		Nodes2[i] = Nodes_initialise(Nodes2[i]) ;
	}
	phylo *tr = phylo_create(ntips) ;
	
	
	int * NNodes = Calloc(NHosts,int) ;
	int * NNodes2 = Calloc(NHosts, int) ;
	minNodes = Calloc(NHosts,int) ;
	
	ch = Calloc(4, int) ;
	int * ov = Calloc(2, int) ;
	
	ei = tr->Nedge - 1 ;
	intNode = 2*tr->NNode  ;
	Anc = calloc(NHosts, sizeof(int)) ;
	vvec = calloc(NS, sizeof(double)) ;
	
	// matrix to store the number of migrating lineages
	double ** mig ;
	int mig_n = 1 ;
	
	
	
  	j = 0;
  	// maybe still some infected left at end:
  
	for (i=0 ; i<NHosts ; i++)
	{
		//SNsum[i] = SN[i] + j ;
    	//j+= SN[i] ;
    	SNsum[i] = OBS[i] + j ;
    	j+= OBS[i] ;
    	st2 = -1 ;
    	for (k=0 ; k< PR_I[TN * (i+1) - 1] ; k++)
    	{
    		Nodes[i] = Nodes_addelement(st2, ST[NHosts - 1], Nodes[i]) ;
  	    	st2-- ;
		}
		Anc[i] = 0 ;
	}


/*******************************ALGORITHM************************************************/

	for (i=TN-1 ; i>-1 ; i--)
	{
		
		ll_i = 0 ;
		ce = 0 ;
		ce1 = 0 ;
    // for each event get hosts
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		
   		
    // if terminal event
   		if (i==Tend[ha])
    	{
    		Rprintf("Terminal event %i\n", i) ;
      		st = (ha>0 ? SNsum[ha-1] : 0) ;
      		st2 = Nodes[ha].minNode - 1;
		  	
		  	for (j=0 ; j<OBS2[ha] ; j++)
		  	{
		  		Nodes2[ha] = Nodes_addelement(st, ST[ha], Nodes2[ha]) ;
				st++ ;	
		  	}
		  	for (j=0 ; j<OBS[ha] ; j++)
		  	{
		  		Nodes[ha] = Nodes_addelement(st, ST[ha], Nodes[ha]) ;
		  		st++ ;
			}		  	
		  	for (j=OBS[ha] + OBS2[ha]; j< -*PR_T2[3][i] ; j++)
		  	{
				//Rprintf("%i, ",st2) ;
				Nodes[ha] = Nodes_addelement(st2, ST[ha], Nodes[ha]) ;
        	  	st2=st2-1 ;
		  	}
		  	//if ((OBS[ha] + OBS2[ha])< -*PR_T2[3][i])
		  	//{
			//	minNodes[ha] = st2+1 ;  
		  	//}
		  	Anc[ha]+=OBS[ha] + OBS2[ha];   
      	}
    	else
		{
			Ancsum = 0 ;
			for (kt = 0 ; kt < NHosts ; kt++)
			{
				Ancsum += Anc[kt] ; 
			}
			
			//if (i > 0 & Ancsum > 1)
			if ((i > 0) & (intNode >= ntips))
			{
				// calculate rates for each host
				kt = *PR_T2[3][i] ;	
				lambdao =  calculateSIR_rates(lambdao, B, PR_I, PR_S, ha, hb, i, kt, NHosts, NS, TN, Anc, OBS2) ;
				lambdao_sum = calculatelambdaosum(lambdao, NHosts) ;
				ll_i = -lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
			
				if ((int) *PR_T2[3][i] == -1 )
				{
					// death but birth in reverse.
					// can calculate the probability of no coalescence has occurred
				
				// minnode
        			mn = (Nodes[ha].N==0) ? 0 : Nodes[ha].minNode ;
					if (mn<= -1) 
					{
						mn=mn-1 ;
					}
					else
					{
						mn= -1 ;
					}
					// add new node
					//minNodes[ha] = mn ;
					Nodes[ha] = Nodes_addelement(mn, *PR_T2[0][i], Nodes[ha]) ;
				} 
				
				//if ((int) *PR_T2[3][i] == 1 )
				//{
					// 1 coalescent event within host
				//	ov = Nodes_pick_sets(ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], ha, hb) ;

				//}
				
				if ((int) *PR_T2[3][i] > 0 )
				{
					Ntrans = (int)*PR_T2[3][i] ;
					// if coincides with coalescent interval among old observed tree
					if (sumOBS > 0) //oa,ob=11
					{			
						if (bt_old[bt_n][0] == *PR_T2[0][i])
						{
							// create coalescent event in new tree
        					Rprintf("birth\n") ;
        					ch[0] = tr_old->edge[bt_n*2 + 1][1] ;
        					ch[1] = tr_old->edge[bt_n*2][1] ;
        					ch[2] = llist_get_ind_i(ch[0], Nodes2[ha].n) ;
        					ch[3] = llist_get_ind_i(ch[1], Nodes2[hb].n) ;
        					
        					// update new tree
        					tr = phylo_addnewedge(tr, ei, intNode, ch[0], llist_get_el_d(ch[2], Nodes2[ha].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr = phylo_addnewedge(tr, ei, intNode, ch[1], llist_get_el_d(ch[3], Nodes2[hb].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr->nodelabel[intNode-ntips] = ha ;
							
							Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
							Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
							Nodes2[ha] = Nodes_addelement(intNode, *PR_T2[0][i], Nodes2[ha]) ;
							intNode -- ;
							Anc[hb] -- ;
							ce++ ;
        				
        				
        				
        				
						
						if (bt_n > 0)
						{
							// replace intNode in tr_old
							for (j = (bt_n-1)*2 ; j>=0 ; j--)
							{
								if (tr_old->edge[j][0]==tr_old->edge[bt_n*2][0])
								{
									tr_old->edge[j][0] = intNode ;			
								}
								if (tr_old->edge[j][1]==tr_old->edge[bt_n*2][0])
								{
									tr_old->edge[j][1] = intNode ;
								}
							}
						}
						bt_n -- ;	
						Ntrans -- ;
						Rprintf("birth\n") ;
        		
        			
        			}
        			
        			// look at migrations
        			}
        		// separate two different kind of transmissions
        		//Ntrans = (sumOBS > 0 & ha!=hb) ? Ntrans - (int) mig2[mig2_n * 3 + (mig_n - 1)] : Ntrans ;
        		for (bn=0 ; bn < Ntrans ; bn++)
        		{
   		  		  	ov = Nodes_pick_sets(ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], ha, hb) ;
            	  	ch = Nodes_pick_children(ch, ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], ha, hb) ;
					//Rprintf("%i\t%i\n", ov[0], ov[1]) ;
					// check for extinct nodes.
				  	if (ch[1] < 0) // -
				  	{
				  		// oa,ob 10 00
				  		Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
				  		
				  	}
				  	else{
				  		// node from hb is observed
					  	if (ch[0] < 0) // +-
					  	{
					  		//oa,ob 00 01
					  		// node from ha is unobserved
					  		//Rprintf("one extinct 2\n") ; 
						  	// move node from hb into ha.
						  	Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
						  	//ttmp = llist_get_el_d(llist_get_ind_i(ch[1], Nodes[hb].n), Nodes[hb].t) ;
						  	//hnode2v = (ov[1] == 0) ? Nodes[hb] : Nodes2[hb] ;
							//hnode2v = Nodes_deleteelement(ch[1], hnode2v) ;
							if (ov[1]==0)
							{
								ttmp = llist_get_el_d(llist_get_ind_i(ch[1], Nodes[hb].n), Nodes[hb].t) ;
						  		
								Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
							}
							else
							{
								ttmp = llist_get_el_d(llist_get_ind_i(ch[1], Nodes2[hb].n), Nodes2[hb].t) ;
						  		
								Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
							
							}
							Nodes[ha] = Nodes_addelement(ch[1], ttmp, Nodes[ha] ) ;
						 	
						 	
						 	if (ha != hb) ce1++ ;
						  	//if (ha != hb) ll += log(lambdao[hb][ha]) ;
						  	Anc[hb]-- ;
						  	Anc[ha]++ ;
						  	//OBS[hb]-- ;
						  	//OBS[ha]++ ;
					  	}
					  	else{ //++
					  		//oa,ob = 01, 10, 00
							// coalescent event
							// create pair of edges with new internal node.
							tr = phylo_addnewedge(tr, ei, intNode, ch[0], llist_get_el_d(ch[2], Nodes[ha].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr = phylo_addnewedge(tr, ei, intNode, ch[1], llist_get_el_d(ch[3], Nodes[hb].t) - *PR_T2[0][i]) ;
							ei-- ;
							
							tr->nodelabel[intNode-ntips] = ha ;
							if (ov[0] == 0)
							{
								Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
							}
							else
							{
								Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
							}
							if (ov[1] == 0)
							{
								Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
							}
							Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
							Nodes[ha] = Nodes_addelement(intNode, *PR_T2[0][i], Nodes[ha]) ;
							intNode -- ;
							Anc[hb] -- ;
							ce++ ;
							//OBS[hb] -- ;
						}
				  }  
        	}		
			}
			if (ce > 0 || ce1 > 0 )
			{
				Ii = PR_I[hb*TN + i] ;
				if (ha == hb)
				{
					ll_i+= log(lambdao[hb][ha]) - log((Anc[hb]+1)*(Anc[hb])/2) ;
					//Rprintf("%8.4f\n", ll_i) ;
				}
				else{
					ll_i+= log(lambdao[hb][ha]) ;
					ktt = *PR_T2[3][i] ;
					Ii = PR_I[ha*TN + i] ;
					ll_i+= log(choose(Ii-(Anc[ha]-ce1), ce1 + ce - ce)*gamma(ce1 + 1))  ;
					ll_i -= log( choose(Anc[hb] + ce1 + ce, ce1 + ce) ) ;
					ll_i -= log( choose(Ii, Ii - (ce1 + ce)) * gamma(ce1 + ce + 1) ) ;
		//			Rprintf("Ii=%i, Ij=%i, ni=%i, m=%i, m'=%i, nj=%i\n", Ii, PR_I[hb*TN + i], Anc[ha]-ce1, ce1+ ce, ce, Anc[hb]+ce1+ce) ;
				}
			}
			/*
			if (ha != hb)
			{
				if (mig_n == 1)
				{
					mig = Calloc(mig_n, double *) ;
					mig[0] = Calloc(5, double) ;
					
				}
				else
				{
					mig = Realloc(mig, mig_n, double *) ;
					mig[mig_n - 1] = Calloc(5, double) ;
				}
				mig[mig_n - 1][0] = *PR_T2[0][i] ; 
				mig[mig_n - 1][1] = ha + 1; 
				mig[mig_n - 1][2] = hb + 1; 
				mig[mig_n - 1][3] = ce ;
				mig[mig_n - 1][4] = ce1 ;
				mig_n++ ;
			}*/
			ll += ll_i ;
			}
			else
			{
				//Rprintf("end \n") ;
			}
			
			
		//	Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[ha], Anc[hb]) ;
			
		
		}
		Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[0][i], *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[0], Anc[1], Anc[2], Anc[3], Anc[4], Nodes[0].N, Nodes[1].N, Ancsum, intNode) ;
			
		
}


/*******************************RESULT TO OUTPUT*****************************************/

	PROTECT(R_List = allocVector(VECSXP ,4)) ;
	PROTECT(R_bt=allocMatrix(REALSXP,tr->NNode,4)) ;
	PROTECT(R_ll = allocVector(REALSXP, 1)) ;
	PROTECT(R_mig = allocMatrix(REALSXP, mig_n-1, 5)) ;
	PR_bt = REAL(R_bt) ;
	bt = phylo_bt(tr, ST, OBS, NHosts) ;
	for (j=0 ; j<tr->NNode ; j++)
	{
		for (i=0;i<4;i++)
		{
			PR_bt[j+i*tr->NNode]=(i>0) ? (bt[j][i] + 1) : bt[j][i] ;
			
		}
	}
	/*
	for (j=0 ; j<(mig_n - 1) ; j++)
	{
		for (i=0;i<5;i++)
		{
			REAL(R_mig)[j+i*(mig_n - 1)]=  mig[j][i] ;
			
		}
	}*/
	
	R_tr = phylo_to_R(tr) ;
	PR_ll = REAL(R_ll) ;
	PR_ll[0] = ll ;
	
	SET_VECTOR_ELT(R_List, 0, R_tr) ;
	SET_VECTOR_ELT(R_List, 1, R_bt) ;
	SET_VECTOR_ELT(R_List, 2, R_ll) ;
	SET_VECTOR_ELT(R_List, 3, R_mig) ;
	UNPROTECT(4) ;
	
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
		free(lambdao[j]) ;
	}
	free(Nodes) ; 
	free(Tend) ;  
	free(NNodes) ;
	free(minNodes) ;
	free(ch) ;
	free(PR_T2) ;
 	free(SNsum) ;
 	free(Anc) ;
 	free(lambdao) ;
 	free(vvec) ;
	return R_List;
}


SEXP tree_reconstruct_withll(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_ktt)
{
   int  st=0, st2=-1, *NNodes, *minNodes, mn, r1, r2, *ch, ei, intNode;
  	int bn , m, ll_calc=1, hi, hj, Ii, Sj, *Anc, kt, Ancsum=0, ce = 0, ce1= 0;
  	double *PR_bt, **bt, v1, v2, *vvec, ll=0.0, lambdao_sum=0.0, *PR_ll, ll_i ;

	item_i * tst ;
	
	SEXP R_tr, R_T, R_List, R_bt, R_ll ;
	hnode *Nodes ;
	
	/******* Interface with R to construct tree from trajectory **************************
	INPUT:	R_sir		list containing I (Infected) S (Susceptible) T (times of events)
			R_NHosts 	Number of hosts
			R_dat		list containing SN (Number of sequences for each host) ST (time of sampling for each host)
			R_B			matrix of rates
			R_NS		number in SIR population for each host
			R_BN		
	
	
	
	*************************************************************************************/
	
	// from R_sir
	SEXP R_Tdim = getAttrib(VECTOR_ELT(R_sir,2), R_DimSymbol) ;
	int TN = INTEGER(R_Tdim)[0] ;		// number of time points of trajectory
	int *PR_I = INTEGER(VECTOR_ELT(R_sir,0)) ;
	int *PR_S = INTEGER(VECTOR_ELT(R_sir,1)) ;	
	int *PR_Iend = INTEGER(VECTOR_ELT(R_sir,3)) ;
	double *PR_T = REAL(coerceVector(VECTOR_ELT(R_sir,2),REALSXP)) ;	
	
	//from R_NHosts
	int NHosts = INTEGER(coerceVector(R_NHosts, INTSXP))[0] ;
	
	//from R_dat
	int *SN = INTEGER(coerceVector(VECTOR_ELT(R_dat,0), INTSXP)) ;
	double *ST = REAL(coerceVector(VECTOR_ELT(R_dat,1), REALSXP)) ;
	
	double *B = REAL(coerceVector(R_B, REALSXP)) ;
	int NS = INTEGER(coerceVector(R_NS, INTSXP))[0] ;
	double BN = REAL(coerceVector(R_BN, REALSXP))[0] ;
	int ktt = INTEGER(coerceVector(R_ktt, INTSXP))[0] ; 
	
	/*************** INTERNAL VARIABLES FOR ALGORITHM ***********************************/
	
	int i, j, k=0, ha, hb ;
	double ***PR_T2 ;
	PR_T2 = Calloc(4,double**) ;	//reshape(PR_T)
	for (j=0 ; j<4 ; j++)
	{
		PR_T2[j] = Calloc(TN,double*) ;
		for (i=0 ; i<TN ; i++)
		{	
		 	PR_T2[j][i] = &PR_T[k++];
		}
	}
	int *Tend = Calloc(NHosts,int) ; // times at which hosts are sampled
	for (i=0 ; i< TN ; i++)
	{
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		if (ha == hb && (*PR_T2[0][i] <= ST[ha]))
		{
			Tend[ha] = i ; 
		}
	}
	
	int *SNsum = calloc(NHosts, sizeof(int)) ;
	int ntips = 0;
	double **lambdao = calloc(NHosts, sizeof(double*)) ;
	for (i=0 ; i<NHosts ; i++)
	{
		ntips+= SN[i] ;
		lambdao[i] = calloc(NHosts, sizeof(double)) ;
	}
	Nodes = Calloc(NHosts,hnode) ; 
	phylo *tr = phylo_create(ntips) ;
	
	
	NNodes = Calloc(NHosts,int) ;
	minNodes = Calloc(NHosts,int) ;
	ch = Calloc(2,int) ;
	
	ei = tr->Nedge - 1 ;
	intNode = 2*tr->NNode  ;
	Anc = calloc(NHosts, sizeof(int)) ;
	vvec = calloc(NS, sizeof(double)) ;
	
	
	
	
  	j = 0;
  	// maybe still some infected left at end:
  
	for (i=0 ; i<NHosts ; i++)
	{
		Nodes[i].t = NULL ;
		Nodes[i].n = NULL ;
   		NNodes[i] = 0 ;
    	minNodes[i] = 0 ;
    	SNsum[i] = SN[i] + j ;
    	j+= SN[i] ;
    	st2 = -1 ;
    	//Rprintf("Iend in host[%i] = %i\n", i, PR_I[TN * (i+1) - 1]) ;
    	for (k=0 ; k< PR_I[TN * (i+1) - 1] ; k++)
    	{
  	    	Nodes[i].n = llist_add_el_i(st2, Nodes[i].n) ;
			Nodes[i].t = llist_add_el_d(ST[NHosts - 1], Nodes[i].t) ;
			//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
        	NNodes[i]++ ;
			st2-- ;
		  	minNodes[i] -- ;         
		}
		Anc[i] = 0 ;
	}


/*******************************ALGORITHM************************************************/

	for (i=TN-1 ; i>-1 ; i--)
	{
		ll_i = 0 ;
		ce = 0 ;
		ce1 = 0 ;
    // for each event get hosts
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		
   		
    // if terminal event
   		if (i==Tend[ha])
    	{
      		st = (ha>0 ? SNsum[ha-1] : 0) ;
      		//NNodes[ha] += -*PR_T2[3][i] ;
      		st2 = minNodes[ha] - 1;
		  	//Rprintf("Nodes in host %i are ", ha) ;
		  	//Rprintf("Host sampled\n") ;
      		for (j=0 ; j<SN[ha] ; j++)
		  	{
				Nodes[ha].n = llist_add_el_i(st, Nodes[ha].n) ;
			  	Nodes[ha].t = llist_add_el_d(ST[ha], Nodes[ha].t) ;
			  	NNodes[ha]++ ;
       			//Rprintf("%i, ",st) ;
			  	st++ ;	
		  	}
		  	for (j=SN[ha] ; j< -*PR_T2[3][i] ; j++)
		  	{
				//Rprintf("%i, ",st2) ;
        	  	Nodes[ha].n = llist_add_el_i(st2, Nodes[ha].n) ;
			  	Nodes[ha].t = llist_add_el_d(ST[ha], Nodes[ha].t) ;
			  	//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
              	NNodes[ha]++ ;
			  	st2=st2-1 ;
		  	}
		  	//Rprintf("\n") ;
		  	if (SN[ha]< -*PR_T2[3][i])
		  	{
				minNodes[ha] = st2+1 ;  
		  	}
		  	Anc[ha]+= SN[ha] ;   
      	}
    	else
		{
			// calculate rates for each host
			Ancsum = 0 ;
			for (kt = 0 ; kt < NHosts ; kt++)
			{
				Ancsum += Anc[kt] ; 
			}
			
			if ((i > 0) & intNode >= ntips)
			//if (i > 0)
			{
			lambdao_sum = 0 ;
			kt = *PR_T2[3][i] ;	
			// hi == hj
			for (hj=0 ; hj < NHosts ; hj++)
			{
				Ii = PR_I[hj*TN + i] ;
				Sj = PR_S[hj*TN + i] ;
				if (hj == hb)
				{
					Ii = Ii - kt ;
					Sj = (kt > 0) ? Sj + kt : Sj;
				}
				lambdao[hj][hj] = B[(NHosts * hj) + hj]  *(double) (Ii * Sj) / NS ; 
				lambdao[hj][hj] = (Anc[hj] >0) ? choose((double)(Anc[hj]), (double)(2)) / choose((double)(PR_I[hj*TN + i]), (double)(2)) *lambdao[hj][hj] : 0;
				lambdao_sum += lambdao[hj][hj] ;
			}
			
			// hi != hj
			for (hj=0 ; hj < NHosts ; hj++)
			{
				Sj = PR_S[hj*TN + i] ;
				
				if (PR_I[hj * TN + i] >= ktt)
				{
					
					vvec[0] = 1 - choose((double)(PR_I[hj*TN + i] - Anc[hj]), (double)(ktt)) / choose((double)(PR_I[hj*TN + i]), (double)(ktt)) ; 
			//		Rprintf("%8.4f\n", vvec[0]) ;
				//	vvec[0] = (double) ktt / PR_I[hj * TN + i] ;
				}
				else 
				vvec[0] = 0;
					
				for (hi=0 ; hi < NHosts ; hi++)
				{
					if (hi != hj)
					{
						Sj = PR_S[hj*TN + i] ;
						Ii = PR_I[hi*TN + i] ;
						if (hj == hb)
						{
							
							Sj = (kt > 0) ? Sj + kt : Sj;
						}
					    //vvec[0] = (Ii >= ktt) ? vvec[0] : 0 ;
						lambdao[hj][hi] = B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / (NS) ; 
						lambdao[hj][hi] = lambdao[hj][hi] * ((Ii >= ktt) ? vvec[0] : 0) ;
						lambdao_sum += lambdao[hj][hi] ;
					}
					else
					{
						//vvec[0] = 0 ;
					}
		//			Rprintf("%8.4f %8.4f ", vvec[0], B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / NS) ;
		//			Rprintf("%8.4f %8.4f ", vvec[0], B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / NS) ;
				
				}
		//		Rprintf("\n") ;
				
			}
			/*
			for (hj = 0 ; hj<NHosts ; hj++)
			{
				Rprintf(" %i ", PR_I[hj*TN + i]) ;
			}
			Rprintf("\n") ;
			*/
			
			//Rprintf("%8.4f\n", 1/lambdao_sum) ;
			//Rprintf("%8.4f\n", ll) ;
			
			//Rprintf("%i\n", kt) ;
		
	/*	
			for (hj=0 ; hj<NHosts ; hj++)
			{
				for (hi=0 ; hi<NHosts ; hi++)
				{
					Rprintf("%8.4f ", lambdao[hj][hi]) ;
				}
				Rprintf("\n") ;
			}
			
			Rprintf("\n") ;  */
			ll_i = -lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
			//ll += ll_i ; //-lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
			
		
    		if ((int) *PR_T2[3][i] == -1 )
			{
				// death but birth in reverse.
				//Rprintf("death\n") ;
				// can calculate the probability of no coalescence has occurred
				
				// minnode
        		mn = (NNodes[ha]==0) ? 0 : minNodes[ha] ;
				if (mn<= -1) 
				{
					mn=mn-1 ;
				}
				else
				{
					mn= -1 ;
				}
				// add new node
				Nodes[ha].n = llist_add_el_i(mn, Nodes[ha].n) ;
				Nodes[ha].t = llist_add_el_d(*PR_T2[0][i], Nodes[ha].t) ;
				NNodes[ha]++ ;
				minNodes[ha] = mn ;
			} 
			if ((int) *PR_T2[3][i] > 0 )
			{
        		//Rprintf("birth\n") ;
        		for (bn=0 ; bn < (int)*PR_T2[3][i] ; bn++)
        		{
          			// birth (coalescent event)
          			// choose  children
  			  		if (ha == hb)
				  	{
            			//Rprintf("got here\n" ) ;  
					  	//Rprintf("\twithin, ") ;
					  	// sampling without replacement
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % (NNodes[ha] - 1) ;
					  	r2 = (r2>=r1) ? r2+1 : r2 ;
					 	// Rprintf("r1,r2 %i,%i",r1,r2) ;
					  	ch[0] = llist_get_el_i(r1, Nodes[ha].n) ;
					  	ch[1] = llist_get_el_i(r2, Nodes[hb].n) ;
					  	//Rprintf("r1,r2 %i,%i",ch[0],ch[1]) ;
        			}
				  	else{
					  	//Rprintf("\tbetween, ") ;
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % NNodes[hb] ;
					  	ch[0]=(NNodes[ha]>1) ? llist_get_el_i(r1, Nodes[ha].n) :  llist_get_el_i(0, Nodes[ha].n);
					  	ch[1]=(NNodes[hb]>1) ? llist_get_el_i(r2, Nodes[hb].n) :  llist_get_el_i(0, Nodes[hb].n);
				  	}
          
				  	// check for extinct nodes.
				  	if (ch[1] < 0) // -
				  	{
				  		//Rprintf("one extinct 1\n") ;
					  	Nodes[hb].n = llist_delete_el_i(r2, Nodes[hb].n) ;
					  	Nodes[hb].t = llist_delete_el_d(r2, Nodes[hb].t) ;
					  	minNodes[hb] = (ch[1]==minNodes[hb]) ? llist_min(Nodes[hb].n) : minNodes[hb] ; 
					  	NNodes[hb]-- ;
				  	}
				  	else{
					  	if (ch[0] < 0) // +-
					  	{
					  		//Rprintf("one extinct 2\n") ; 
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
						  	if (ha != hb) ce1++ ;
						  	//if (ha != hb) ll += log(lambdao[hb][ha]) ;
						  	Anc[hb]-- ;
						  	Anc[ha]++ ;
					  	}
					  	else{ //++
							// coalescent event
							// create pair of edges with new internal node.
							tr->edge[ei][0] = intNode ;
							tr->edge[ei][1] = ch[0] ;
							tr->el[ei] = llist_get_el_d(r1, Nodes[ha].t) - *PR_T2[0][i] ;
							//Rprintf("event %i %i %8.4f %8.4f %8.4f\n",intNode, ch[0], llist_get_el_d(r1, Nodes[ha].t), *PR_T2[0][i], tr->el[ei]) ;
							//Rprintf("event %i %i %i %8.4f\n",ha+1, ha+1, hb+1, *PR_T2[0][i]) ;
						
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
							//ll += log(lambdao[hb][ha]) ;
							intNode -- ;
							Anc[hb] -- ;
							ce++ ;
						}
				  }  
        	}		
			}
			if (ce > 0 || ce1 > 0 )
			{
				Ii = PR_I[hb*TN + i] ;
				//Ii = Anc[hb]  + 1 ;
				//if (ha == hb)
				//	Ii = Ii - kt ;
				//Sj = (kt > 0) ? Sj + kt : Sj;
						
				//ll_i+= log(B[(NHosts * hb) + ha]  *(double) (Ii * Sj) / NS ); 
			if (ha == hb)
			{
				ll_i+= log(lambdao[hb][ha]) - log((Anc[hb]+1)*(Anc[hb])/2) ;
				//Rprintf("%8.4f\n", ll_i) ;
			}else{
				ll_i+= log(lambdao[hb][ha]) ;
				//Rprintf("%8.4f\n", ll_i) ;
				//ktt = *PR_T2[3][i] ;
				Ii = PR_I[ha*TN + i] ;
				
				ll_i+= log(choose(Ii-(Anc[ha]-ce1), ce1 + ce - ce)*gamma(ce1 + 1))  ;
				ll_i -= log( choose(Anc[hb] + ce1 + ce, ce1 + ce) ) ;
				ll_i -= log( choose(Ii, Ii - (ce1 + ce)) * gamma(ce1 + ce + 1) ) ;
	//			Rprintf("Ii=%i, Ij=%i, ni=%i, m=%i, m'=%i, nj=%i\n", Ii, PR_I[hb*TN + i], Anc[ha]-ce1, ce1+ ce, ce, Anc[hb]+ce1+ce) ;
			 	//ll_i-= (log(choose((double)(Ii), (double)(ktt))  -  choose((double)(Ii - Anc[hb] - 1), (double)(ktt))) ) ;
			 	//ll_i -= (log(choose((double)(PR_I[ha*TN + i]), (double)(ktt)) * gamma(ktt + 1)))   ;
			
			//ll_i+= dhyper(1, Anc[hb]+1, Ii - Anc[hb] - 1, ktt, 1) ;
			
				//ll_i-=  log(choose(Ii, *PR_T2[3][i])-choose(Ii - (Anc[hb] + 1), *PR_T2[3][i])) ; 
				//ll_i+= log(*PR_T2[3][i] / Ii) ;
				//ll_i-= log(PR_I[ha*TN + i]) ;
			//	Rprintf("%8.4f\t%i\t%i\t%8.4f\n", ll_i, Ii, Anc[hb], *PR_T2[3][i]) ;
				//ll_i-=  log(choose(PR_I[ha*TN + i], PR_I[ha*TN + i] - *PR_T2[3][i] )) ; 
				//Rprintf("%8.4f\n", ll_i) ;
			}
			//ll_i = log(1 - exp(ll_i)) + log(lambdao[hb][ha]) - log(lambdao_sum);
			}
			//printf("ll_i = %8.4f\n", ll_i) ; 
			//if (Ancsum > 1) ll+= ll_i ;
			ll += ll_i ;
			}
			else
			{
				//Rprintf("end \n") ;
			}
			
			
			//Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[ha], Anc[hb], PR_I[0*TN + i], PR_I[1*TN + i], PR_I[2*TN + i], PR_I[3*TN + i], PR_I[4*TN + i] , ktt ) ;
			
		
		}
		//Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[0], Anc[1], Anc[2], Anc[3], Anc[4], Ancsum, intNode) ;
			
		
}


/*******************************RESULT TO OUTPUT*****************************************/

	PROTECT(R_List = allocVector(VECSXP ,3)) ;
	PROTECT(R_bt=allocMatrix(REALSXP,tr->NNode,4)) ;
	PROTECT(R_ll = allocVector(REALSXP, 1)) ;
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
	PR_ll = REAL(R_ll) ;
	PR_ll[0] = ll ;
	
	SET_VECTOR_ELT(R_List, 0, R_tr) ;
	SET_VECTOR_ELT(R_List, 1, R_bt) ;
	SET_VECTOR_ELT(R_List, 2, R_ll) ;
	UNPROTECT(3) ;
	
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
		free(lambdao[j]) ;
	}
	free(Nodes) ; 
	free(Tend) ;  
	free(NNodes) ;
	free(minNodes) ;
	free(ch) ;
	free(PR_T2) ;
 	free(SNsum) ;
 	free(Anc) ;
 	free(lambdao) ;
 	free(vvec) ;
	return R_List;
}

SEXP tree_reconstruct_withll2(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN)
{
   int  st=0, st2=-1, *NNodes, *minNodes, mn, r1, r2, *ch, ei, intNode;
  	int bn , m, ll_calc=1, hi, hj, Ii, Sj, *Anc, kt, ktt, Ancsum=0, ce = 0, ce1= 0;
  	double *PR_bt, **bt, v1, v2, *vvec, ll=0.0, lambdao_sum=0.0, *PR_ll, ll_i ;

	item_i * tst ;
	
	SEXP R_tr, R_T, R_List, R_bt, R_ll ;
	hnode *Nodes ;
	
	/******* Interface with R to construct tree from trajectory **************************
	INPUT:	R_sir		list containing I (Infected) S (Susceptible) T (times of events)
			R_NHosts 	Number of hosts
			R_dat		list containing SN (Number of sequences for each host) ST (time of sampling for each host)
			R_B			matrix of rates
			R_NS		number in SIR population for each host
			R_BN		
	
	
	
	*************************************************************************************/
	
	// from R_sir
	SEXP R_Tdim = getAttrib(VECTOR_ELT(R_sir,2), R_DimSymbol) ;
	int TN = INTEGER(R_Tdim)[0] ;		// number of time points of trajectory
	int *PR_I = INTEGER(VECTOR_ELT(R_sir,0)) ;
	int *PR_S = INTEGER(VECTOR_ELT(R_sir,1)) ;	
	int *PR_Iend = INTEGER(VECTOR_ELT(R_sir,3)) ;
	double *PR_T = REAL(coerceVector(VECTOR_ELT(R_sir,2),REALSXP)) ;	
	
	//from R_NHosts
	int NHosts = INTEGER(coerceVector(R_NHosts, INTSXP))[0] ;
	
	//from R_dat
	int *SN = INTEGER(coerceVector(VECTOR_ELT(R_dat,0), INTSXP)) ;
	double *ST = REAL(coerceVector(VECTOR_ELT(R_dat,1), REALSXP)) ;
	
	double *B = REAL(coerceVector(R_B, REALSXP)) ;
	int NS = INTEGER(coerceVector(R_NS, INTSXP))[0] ;
	double BN = REAL(coerceVector(R_BN, REALSXP))[0] ;
	
	/*************** INTERNAL VARIABLES FOR ALGORITHM ***********************************/
	
	int i, j, k=0, ha, hb ;
	double ***PR_T2 ;
	PR_T2 = Calloc(4,double**) ;	//reshape(PR_T)
	for (j=0 ; j<4 ; j++)
	{
		PR_T2[j] = Calloc(TN,double*) ;
		for (i=0 ; i<TN ; i++)
		{	
		 	PR_T2[j][i] = &PR_T[k++];
		}
	}
	int *Tend = Calloc(NHosts,int) ; // times at which hosts are sampled
	for (i=0 ; i< TN ; i++)
	{
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		if (ha == hb && (*PR_T2[0][i] <= ST[ha]))
		{
			Tend[ha] = i ; 
		}
	}
	
	int *SNsum = calloc(NHosts, sizeof(int)) ;
	int ntips = 0;
	double **lambdao = calloc(NHosts, sizeof(double*)) ;
	for (i=0 ; i<NHosts ; i++)
	{
		ntips+= SN[i] ;
		lambdao[i] = calloc(NHosts, sizeof(double)) ;
	}
	Nodes = Calloc(NHosts,hnode) ; 
	phylo *tr = phylo_create(ntips) ;
	
	
	NNodes = Calloc(NHosts,int) ;
	minNodes = Calloc(NHosts,int) ;
	ch = Calloc(2,int) ;
	
	ei = tr->Nedge - 1 ;
	intNode = 2*tr->NNode  ;
	Anc = calloc(NHosts, sizeof(int)) ;
	vvec = calloc(NS, sizeof(double)) ;
	
	
	
	
  	j = 0;
  	// maybe still some infected left at end:
  
	for (i=0 ; i<NHosts ; i++)
	{
		Nodes[i].t = NULL ;
		Nodes[i].n = NULL ;
   		NNodes[i] = 0 ;
    	minNodes[i] = 0 ;
    	SNsum[i] = SN[i] + j ;
    	j+= SN[i] ;
    	st2 = -1 ;
    	//Rprintf("Iend in host[%i] = %i\n", i, PR_I[TN * (i+1) - 1]) ;
    	for (k=0 ; k< PR_I[TN * (i+1) - 1] ; k++)
    	{
  	    	Nodes[i].n = llist_add_el_i(st2, Nodes[i].n) ;
			Nodes[i].t = llist_add_el_d(ST[NHosts - 1], Nodes[i].t) ;
			//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
        	NNodes[i]++ ;
			st2-- ;
		  	minNodes[i] -- ;         
		}
		Anc[i] = 0 ;
	}


/*******************************ALGORITHM************************************************/

	for (i=TN-1 ; i>-1 ; i--)
	{
		ll_i = 0 ;
		ce = 0 ;
		ce1 = 0 ;
    // for each event get hosts
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		
   		
    // if terminal event
   		if (i==Tend[ha])
    	{
      		st = (ha>0 ? SNsum[ha-1] : 0) ;
      		//NNodes[ha] += -*PR_T2[3][i] ;
      		st2 = minNodes[ha] - 1;
		  	//Rprintf("Nodes in host %i are ", ha) ;
		  	//Rprintf("Host sampled\n") ;
      		for (j=0 ; j<SN[ha] ; j++)
		  	{
				Nodes[ha].n = llist_add_el_i(st, Nodes[ha].n) ;
			  	Nodes[ha].t = llist_add_el_d(ST[ha], Nodes[ha].t) ;
			  	NNodes[ha]++ ;
       			//Rprintf("%i, ",st) ;
			  	st++ ;	
		  	}
		  	for (j=SN[ha] ; j< -*PR_T2[3][i] ; j++)
		  	{
				//Rprintf("%i, ",st2) ;
        	  	Nodes[ha].n = llist_add_el_i(st2, Nodes[ha].n) ;
			  	Nodes[ha].t = llist_add_el_d(ST[ha], Nodes[ha].t) ;
			  	//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
              	NNodes[ha]++ ;
			  	st2=st2-1 ;
		  	}
		  	//Rprintf("\n") ;
		  	if (SN[ha]< -*PR_T2[3][i])
		  	{
				minNodes[ha] = st2+1 ;  
		  	}
		  	Anc[ha]+= SN[ha] ;   
      	}
    	else
		{
			// calculate rates for each host
			Ancsum = 0 ;
			for (kt = 0 ; kt < NHosts ; kt++)
			{
				Ancsum += Anc[kt] ; 
			}
			
			if (i > 0 & Ancsum > 1)
			//if (i > 0)
			{
			lambdao_sum = 0 ;
			kt = *PR_T2[3][i] ;	
			// hi == hj
			for (hj=0 ; hj < NHosts ; hj++)
			{
				Ii = PR_I[hj*TN + i] ;
				Sj = PR_S[hj*TN + i] ;
				if (hj == hb)
				{
					Ii = Ii - kt ;
					Sj = (kt > 0) ? Sj + kt : Sj;
				}
				lambdao[hj][hj] = B[(NHosts * hj) + hj]  *(double) (Ii * Sj) / NS ; 
				lambdao[hj][hj] = (Anc[hj] >0) ? choose((double)(Anc[hj]), (double)(2)) / choose((double)(PR_I[hj*TN + i]), (double)(2)) *lambdao[hj][hj] : 0;
				lambdao_sum += lambdao[hj][hj] ;
			}
			
			// hi != hj
			for (hj=0 ; hj < NHosts ; hj++)
			{
				Sj = PR_S[hj*TN + i] ;
				ktt = 100 ;
				if (PR_I[hj * TN + i] >= ktt)
				{
					
					vvec[0] = 1 - choose((double)(PR_I[hj*TN + i] - Anc[hj]), (double)(ktt)) / choose((double)(PR_I[hj*TN + i]), (double)(ktt)) ; 
				//	vvec[0] = (double) ktt / PR_I[hj * TN + i] ;
				}
				else 
				vvec[0] = 0;
					
				for (hi=0 ; hi < NHosts ; hi++)
				{
				
					if (hi != hj)
					{
						Ii = PR_I[hi*TN + i] ;
						if (hj == hb)
						{
							
							Sj = (kt > 0) ? Sj + kt : Sj;
						}
					    //vvec[0] = (Ii >= ktt) ? vvec[0] : 0 ;
						lambdao[hj][hi] = B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / (NS * ktt) ; 
						lambdao[hj][hi] = lambdao[hj][hi] * ((Ii >= ktt) ? vvec[0] : 0) ;
						lambdao_sum += lambdao[hj][hi] ;
					}
					else
					{
						//vvec[0] = 0 ;
					}
		//			Rprintf("%8.4f %8.4f ", vvec[0], B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / NS) ;
					//Rprintf("%8.4f ", vvec[0], B[(NHosts * hj) + hi]  *(double) (Ii * Sj) / NS) ;
				
				}
		//		Rprintf("\n") ;
				
			}
			/*
			for (hj = 0 ; hj<NHosts ; hj++)
			{
				Rprintf(" %i ", PR_I[hj*TN + i]) ;
			}
			Rprintf("\n") ;
			*/
			
			//Rprintf("%8.4f\n", 1/lambdao_sum) ;
			//Rprintf("%8.4f\n", ll) ;
			
			//Rprintf("%i\n", kt) ;
		
	/*	
			for (hj=0 ; hj<NHosts ; hj++)
			{
				for (hi=0 ; hi<NHosts ; hi++)
				{
					Rprintf("%8.4f ", lambdao[hj][hi]) ;
				}
				Rprintf("\n") ;
			}
			
			Rprintf("\n") ;  */
			ll_i = -lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
			//ll += ll_i ; //-lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
			
		
    		if ((int) *PR_T2[3][i] == -1 )
			{
				// death but birth in reverse.
				//Rprintf("death\n") ;
				// can calculate the probability of no coalescence has occurred
				
				// minnode
        		mn = (NNodes[ha]==0) ? 0 : minNodes[ha] ;
				if (mn<= -1) 
				{
					mn=mn-1 ;
				}
				else
				{
					mn= -1 ;
				}
				// add new node
				Nodes[ha].n = llist_add_el_i(mn, Nodes[ha].n) ;
				Nodes[ha].t = llist_add_el_d(*PR_T2[0][i], Nodes[ha].t) ;
				NNodes[ha]++ ;
				minNodes[ha] = mn ;
			} 
			if ((int) *PR_T2[3][i] > 0 )
			{
        		//Rprintf("birth\n") ;
        		for (bn=0 ; bn < (int)*PR_T2[3][i] ; bn++)
        		{
          			// birth (coalescent event)
          			// choose  children
  			  		if (ha == hb)
				  	{
            			//Rprintf("got here\n" ) ;  
					  	//Rprintf("\twithin, ") ;
					  	// sampling without replacement
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % (NNodes[ha] - 1) ;
					  	r2 = (r2>=r1) ? r2+1 : r2 ;
					 	// Rprintf("r1,r2 %i,%i",r1,r2) ;
					  	ch[0] = llist_get_el_i(r1, Nodes[ha].n) ;
					  	ch[1] = llist_get_el_i(r2, Nodes[hb].n) ;
					  	//Rprintf("r1,r2 %i,%i",ch[0],ch[1]) ;
        			}
				  	else{
					  	//Rprintf("\tbetween, ") ;
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % NNodes[hb] ;
					  	ch[0]=(NNodes[ha]>1) ? llist_get_el_i(r1, Nodes[ha].n) :  llist_get_el_i(0, Nodes[ha].n);
					  	ch[1]=(NNodes[hb]>1) ? llist_get_el_i(r2, Nodes[hb].n) :  llist_get_el_i(0, Nodes[hb].n);
				  	}
          
				  	// check for extinct nodes.
				  	if (ch[1] < 0) // -
				  	{
				  		//Rprintf("one extinct 1\n") ;
					  	Nodes[hb].n = llist_delete_el_i(r2, Nodes[hb].n) ;
					  	Nodes[hb].t = llist_delete_el_d(r2, Nodes[hb].t) ;
					  	minNodes[hb] = (ch[1]==minNodes[hb]) ? llist_min(Nodes[hb].n) : minNodes[hb] ; 
					  	NNodes[hb]-- ;
				  	}
				  	else{
					  	if (ch[0] < 0) // +-
					  	{
					  		//Rprintf("one extinct 2\n") ; 
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
						  	if (ha != hb) ce1++ ;
						  	//if (ha != hb) ll += log(lambdao[hb][ha]) ;
						  	Anc[hb]-- ;
						  	Anc[ha]++ ;
					  	}
					  	else{ //++
							// coalescent event
							// create pair of edges with new internal node.
							tr->edge[ei][0] = intNode ;
							tr->edge[ei][1] = ch[0] ;
							tr->el[ei] = llist_get_el_d(r1, Nodes[ha].t) - *PR_T2[0][i] ;
							//Rprintf("event %i %i %8.4f %8.4f %8.4f\n",intNode, ch[0], llist_get_el_d(r1, Nodes[ha].t), *PR_T2[0][i], tr->el[ei]) ;
							//Rprintf("event %i %i %i %8.4f\n",ha+1, ha+1, hb+1, *PR_T2[0][i]) ;
						
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
							//ll += log(lambdao[hb][ha]) ;
							intNode -- ;
							Anc[hb] -- ;
							ce++ ;
						}
				  }  
        	}		
			}
			if (ce > 0 || ce1 > 0 )
			{
				Ii = PR_I[hb*TN + i] ;
				//Ii = Anc[hb]  + 1 ;
				//if (ha == hb)
				//	Ii = Ii - kt ;
				//Sj = (kt > 0) ? Sj + kt : Sj;
						
				//ll_i+= log(B[(NHosts * hb) + ha]  *(double) (Ii * Sj) / NS ); 
			if (ha == hb)
			{
				ll_i+= log(lambdao[hb][ha]) - log((Anc[hb]+1)*(Anc[hb])/2) ;
				//Rprintf("%8.4f\n", ll_i) ;
			}else{
				ll_i+= log(lambdao[hb][ha]) ;
				//Rprintf("%8.4f\n", ll_i) ;
				ktt = *PR_T2[3][i] ;
				Ii = PR_I[ha*TN + i] ;
				
				ll_i+= log(choose(Ii-(Anc[ha]-ce1), ce1 + ce - ce)*gamma(ce1 + 1))  ;
				ll_i -= log( choose(Anc[hb] + ce1 + ce, ce1 + ce) ) ;
				ll_i -= log( choose(Ii, Ii - (ce1 + ce)) * gamma(ce1 + ce + 1) ) ;
	//			Rprintf("Ii=%i, Ij=%i, ni=%i, m=%i, m'=%i, nj=%i\n", Ii, PR_I[hb*TN + i], Anc[ha]-ce1, ce1+ ce, ce, Anc[hb]+ce1+ce) ;
			 	//ll_i-= (log(choose((double)(Ii), (double)(ktt))  -  choose((double)(Ii - Anc[hb] - 1), (double)(ktt))) ) ;
			 	//ll_i -= (log(choose((double)(PR_I[ha*TN + i]), (double)(ktt)) * gamma(ktt + 1)))   ;
			
			//ll_i+= dhyper(1, Anc[hb]+1, Ii - Anc[hb] - 1, ktt, 1) ;
			
				//ll_i-=  log(choose(Ii, *PR_T2[3][i])-choose(Ii - (Anc[hb] + 1), *PR_T2[3][i])) ; 
				//ll_i+= log(*PR_T2[3][i] / Ii) ;
				//ll_i-= log(PR_I[ha*TN + i]) ;
			//	Rprintf("%8.4f\t%i\t%i\t%8.4f\n", ll_i, Ii, Anc[hb], *PR_T2[3][i]) ;
				//ll_i-=  log(choose(PR_I[ha*TN + i], PR_I[ha*TN + i] - *PR_T2[3][i] )) ; 
				//Rprintf("%8.4f\n", ll_i) ;
			}
			//ll_i = log(1 - exp(ll_i)) + log(lambdao[hb][ha]) - log(lambdao_sum);
			}
			//printf("ll_i = %8.4f\n", ll_i) ; 
			//if (Ancsum > 1) ll+= ll_i ;
			ll += ll_i ;
			}
			else
			{
				//Rprintf("end \n") ;
			}
			
			
	//		Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[ha], Anc[hb]) ;
			
		
		}
		//Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[0], Anc[1], Anc[2], Anc[3], Anc[4], Ancsum, intNode) ;
			
		
}


/*******************************RESULT TO OUTPUT*****************************************/

	PROTECT(R_List = allocVector(VECSXP ,3)) ;
	PROTECT(R_bt=allocMatrix(REALSXP,tr->NNode,4)) ;
	PROTECT(R_ll = allocVector(REALSXP, 1)) ;
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
	PR_ll = REAL(R_ll) ;
	PR_ll[0] = ll ;
	
	SET_VECTOR_ELT(R_List, 0, R_tr) ;
	SET_VECTOR_ELT(R_List, 1, R_bt) ;
	SET_VECTOR_ELT(R_List, 2, R_ll) ;
	UNPROTECT(3) ;
	
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
		free(lambdao[j]) ;
	}
	free(Nodes) ; 
	free(Tend) ;  
	free(NNodes) ;
	free(minNodes) ;
	free(ch) ;
	free(PR_T2) ;
 	free(SNsum) ;
 	free(Anc) ;
 	free(lambdao) ;
 	free(vvec) ;
	return R_List;
}


SEXP tree_reconstruct(SEXP R_sir, SEXP R_NHosts, SEXP R_dat)
{
  int *PR_I, *PR_S, *PR_Iend, *SN, NHosts, i, j, k,st=0, st2=-1, TN, ha, hb, *Tend, *NNodes, *minNodes, mn, r1, r2, *ch, ntips=0, ei, intNode;
  int bn , *SNsum, m, ll_calc=1;
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
	SNsum = calloc(NHosts, sizeof(int)) ;
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
	
	
	
  	j = 0;
  	// maybe still some infected left at end:
  
	for (i=0 ; i<NHosts ; i++)
	{
		Nodes[i].t = NULL ;
		Nodes[i].n = NULL ;
   		NNodes[i] = 0 ;
    	minNodes[i] = 0 ;
    	SNsum[i] = SN[i] + j ;
    	j+= SN[i] ;
    	st2 = -1 ;
    	//Rprintf("Iend in host[%i] = %i\n", i, PR_I[TN * (i+1) - 1]) ;
    	for (k=0 ; k< PR_I[TN * (i+1) - 1] ; k++)
    	{
  	    	Nodes[i].n = llist_add_el_i(st2, Nodes[i].n) ;
			Nodes[i].t = llist_add_el_d(ST[NHosts - 1], Nodes[i].t) ;
			//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
        	NNodes[i]++ ;
			st2-- ;
		  	minNodes[i] -- ;         
		}
	}

	// Rprintf("Nodes 0 %i, ", (Nodes[0].n)->val) ;
	// Find last time point for each host
	for (i=0 ; i< TN ; i++)
	{
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		if (ha == hb && (*PR_T2[0][i] <= ST[ha]))
		{
			Tend[ha] = i ; 
		}
	}
	
  
	for (i=TN-1 ; i>-1 ; i--)
	{
    // for each event get hosts
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		// print out all info
		//llist_print(Nodes[0].n) ;
		//Rprintf("\n") ;
   		//Rprintf("%8.4f %i %i ", *PR_T2[0][i], ha+1, hb+1) ;
   		for (j=0 ; j<NHosts ; j++)
   		{
   			//Rprintf("N[%i]=%i, ",j+1,  NNodes[j]) ;
   		}
   		
   		for (j=0 ; j<NHosts ; j++)
   		{
   			//Rprintf("min[%i]=%i, ",j+1,  minNodes[j]) ;
   		}
   		
    // if terminal event
   		if (i==Tend[ha])
    	{
      		st = (ha>0 ? SNsum[ha-1] : 0) ;
      		//NNodes[ha] += -*PR_T2[3][i] ;
      		st2 = minNodes[ha] - 1;
		  	//Rprintf("Nodes in host %i are ", ha) ;
		  	//Rprintf("Host sampled\n") ;
      		for (j=0 ; j<SN[ha] ; j++)
		  	{
				Nodes[ha].n = llist_add_el_i(st, Nodes[ha].n) ;
			  	Nodes[ha].t = llist_add_el_d(ST[ha], Nodes[ha].t) ;
			  	NNodes[ha]++ ;
       			//Rprintf("%i, ",st) ;
			  	st++ ;	
		  	}
		  	for (j=SN[ha] ; j< -*PR_T2[3][i] ; j++)
		  	{
				//Rprintf("%i, ",st2) ;
        	  	Nodes[ha].n = llist_add_el_i(st2, Nodes[ha].n) ;
			  	Nodes[ha].t = llist_add_el_d(ST[ha], Nodes[ha].t) ;
			  	//Rprintf("Nodes i %i, ", (Nodes[i].n)->val) ;	
              	NNodes[ha]++ ;
			  	st2=st2-1 ;
		  	}
		  	//Rprintf("\n") ;
		  	if (SN[ha]< -*PR_T2[3][i])
		  	{
				minNodes[ha] = st2+1 ;  
		  	}   
      	}
    	else
		{
    		if ((int) *PR_T2[3][i] == -1 )
			{
				// death but birth in reverse.
				//Rprintf("death\n") ;
				// can calculate the probability of no coalescence has occurred
				
				// minnode
        		mn = (NNodes[ha]==0) ? 0 : minNodes[ha] ;
				if (mn<= -1) 
				{
					mn=mn-1 ;
				}
				else
				{
					mn= -1 ;
				}
				// add new node
				Nodes[ha].n = llist_add_el_i(mn, Nodes[ha].n) ;
				Nodes[ha].t = llist_add_el_d(*PR_T2[0][i], Nodes[ha].t) ;
				NNodes[ha]++ ;
				minNodes[ha] = mn ;
			} 
			if ((int) *PR_T2[3][i] > 0 )
			{
        		//Rprintf("birth\n") ;
        		for (bn=0 ; bn < (int)*PR_T2[3][i] ; bn++)
        		{
          			// birth (coalescent event)
          			// choose  children
  			  		if (ha == hb)
				  	{
            			//Rprintf("got here\n" ) ;  
					  	//Rprintf("\twithin, ") ;
					  	// sampling without replacement
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % (NNodes[ha] - 1) ;
					  	r2 = (r2>=r1) ? r2+1 : r2 ;
					 	// Rprintf("r1,r2 %i,%i",r1,r2) ;
					  	ch[0] = llist_get_el_i(r1, Nodes[ha].n) ;
					  	ch[1] = llist_get_el_i(r2, Nodes[hb].n) ;
					  	//Rprintf("r1,r2 %i,%i",ch[0],ch[1]) ;
        			}
				  	else{
					  	//Rprintf("\tbetween, ") ;
					  	r1 = rand() % NNodes[ha] ;
					  	r2 = rand() % NNodes[hb] ;
					  	ch[0]=(NNodes[ha]>1) ? llist_get_el_i(r1, Nodes[ha].n) :  llist_get_el_i(0, Nodes[ha].n);
					  	ch[1]=(NNodes[hb]>1) ? llist_get_el_i(r2, Nodes[hb].n) :  llist_get_el_i(0, Nodes[hb].n);
				  	}
          
				  	// check for extinct nodes.
				  	if (ch[1] < 0) // -
				  	{
				  		//Rprintf("one extinct 1\n") ;
					  	Nodes[hb].n = llist_delete_el_i(r2, Nodes[hb].n) ;
					  	Nodes[hb].t = llist_delete_el_d(r2, Nodes[hb].t) ;
					  	minNodes[hb] = (ch[1]==minNodes[hb]) ? llist_min(Nodes[hb].n) : minNodes[hb] ; 
					  	NNodes[hb]-- ;
				  	}
				  	else{
					  	if (ch[0] < 0) // +-
					  	{
					  		//Rprintf("one extinct 2\n") ; 
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
					  	else{ //++
							// coalescent event
							// create pair of edges with new internal node.
							tr->edge[ei][0] = intNode ;
							tr->edge[ei][1] = ch[0] ;
							tr->el[ei] = llist_get_el_d(r1, Nodes[ha].t) - *PR_T2[0][i] ;
							//Rprintf("event %i %i %8.4f %8.4f %8.4f\n",intNode, ch[0], llist_get_el_d(r1, Nodes[ha].t), *PR_T2[0][i], tr->el[ei]) ;
							//Rprintf("event %i %i %i %8.4f\n",ha+1, ha+1, hb+1, *PR_T2[0][i]) ;
						
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
 	free(SNsum) ;
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

void bnsizeupdate(int * bnsizes, int NHosts, int opt, int * opt_val, double opt_val2, int nS, int K)
{
  int i ;
  
  switch (opt)
  {
    case 0: 
    
      // bnsize of opt_val ;
      for (i=0 ; i<NHosts ; i++)
      {
        bnsizes[i] = K ;
        
      }
    	break ;
    case 1:
    
      for (i=0 ; i<NHosts ; i++)
      {
        bnsizes[i] = (int) rbinom((double) opt_val[i], opt_val2) ;
        
      }
    
    break ;
    case 2:
    
      for (i=0 ; i<NHosts ; i++)
      {
        bnsizes[i] = (int) rbinom((double) nS, opt_val2) ;
        
      } 
    
    break ;
     
  }
  
  
}

void cSIR_iters_ST(int *n, int I0, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN, double bnprob, int K)
{
	int **pp1,i,ot=0, *bnsizes, opt=0;
	double **pp2,tmax=-1 ;
  bnsizes = calloc(NHosts, sizeof(int)) ;
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
    	//bnsizeupdate(bnsizes, NHosts, opt, &pp1[1][(n[0]-1)*NHosts], bnprob, nS) ;
  		bnsizeupdate(bnsizes, NHosts, opt, &pp1[1][(n[0]-1)*NHosts], bnprob, nS, K) ;
  		//Rprintf("%i\n ", bnsizes[0]) ; 
		cSIR_iter_ST(n, nS, NHosts, B, dr, pp1, pp2, ST, SN, bnsizes) ;

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
  free(bnsizes) ;
}

/*
int getindofhostpair(int **bh, int ha, int hb, int N)
{
  int i, found=0, ind=-1;
  while ( i< N && found==0 )
  {
    if ((bh[i][0] == ha && bh[i][1] == hb) || (bh[i][1] == ha && bh[i][0] == hb))
    {
      found = 1 ;
      ind = i ;
    }
    i++ ;
  }
  return(ind) ;
}

void cSIR_time_samples(double *tpts, int *lo, int *lo_h, int *s, int ntips, int Ntpt)
{
	// find topology given new sample and current topology
  int i, ha, hb, ind ;
  int *l , **bh, *lo_new, *s_new, *h;
  
  // initialisation
  
  // ordered matrix of coalescent events (most recent first)
  bh = calloc((ntips-1), sizeof(int *)) ;
  for (i=0 ; i < (ntips-1) ; i++)
  {
    bh[i] = calloc(3, sizeof(int)) ;
    bh[i][0] = lo_h[s[i]] ;
    bh[i][1] = lo_h[s[i] + 1] ;
    bh[i][2] = s[i] ;
  }
  
  h = calloc(ntips, sizeof(int)) ;
  for (i=0 ; )
  
  lo_new = calloc(ntips, sizeof(int)) ;
  s_new = calloc(ntips - 1, sizeof(int)) ;
  
  
  for (i=Ntpt-1 ; i>=0 ; i--)
  {
    // find hosts of coalescent event
    ha = tpts[2*Ntpt + i] ;
    hb = tpts[3*Ntpt + i] ;
    
    // find first corresponding event in current structure
    ind = getindofhostpair(bh, ha, hb, ntips-1) ;
    if (ind >= 0)
    {
      // if exists
      s_new[Ntpt - 1 -i] = bh[ind,2]
    }else{
      
    }
    
  }
}
*/

void cSIR_iter(int n, int nS, int NHosts, double *B, double dr, int **p1, double **p2)
{
  /********
  Returns SIR trajectory
  
  
  *********/
  
	double *R, sumR=0, t_n, Rtot, *T, rand_e;	// rate matrix
	R = Calloc(NHosts * (NHosts+1),double) ;
	int i,j,e,ec,er, *S, *I, In, Sn, np;
	
	S=p1[0];
	I=p1[1];
	T=p2[0];
	
	for (i=0; i<NHosts; i++)
	{
    // death rates
		R[i + NHosts*NHosts]=I[(n-1)*NHosts + i]*dr ;
		sumR+= R[i + NHosts*NHosts] ;
		for (j=0; j<NHosts; j++)
		{
      In = I[(n-1)*NHosts + j] ;
      Sn = S[(n-1)*NHosts + i] ;
			
      R[i*NHosts + j] = In * Sn * B[i * NHosts + j]/(double)(nS) ;
			sumR+= R[i*NHosts + j] ;
		}
	}
	
	for (i=0; i<(NHosts*(NHosts+1));i++)
	{
		R[i]=R[i]/sumR ;
	}
	
	// Draw time of next event
	t_n = rexp(sumR) ;
	
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
    np = 1 ;
		I[n*NHosts + ec] = I[n*NHosts + ec] + np;
		S[n*NHosts + ec] = S[n*NHosts + ec] - np;
		T[n*4+3] = np ;
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

void cSIR_iter_ST(int *np, int nS, int NHosts, double *B, double dr, int **p1, double **p2, double *ST, int *SN, int * bnsizes)
{
	double *R, sumR=0, t_n, Rtot, *T, rand_e, *ST2, T_n, t_n1;	// rate matrix
	R = Calloc(NHosts * (NHosts+1),double) ;
	ST2 = Calloc(NHosts,double) ;
	int i, j,e,ec,er, *S, *I, *Iend, n, h_i, In, Sn, bn, bnsize, all_out=0;
	
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
      In = I[(n-1)*NHosts + j] ;
      Sn = S[(n-1)*NHosts + i] ;
      //Rprintf("%i\t",In) ;
  		if (i!=j)
      {
        bnsize = bnsizes[j] ;
        In = (In >= bnsize) ? In : 0 ;
        Sn = (Sn >= bnsize) ? Sn : 0 ;
      //	R[i*NHosts + j] = In * Sn * B[i * NHosts + j]/((double)(nS) * bnsize) ;
	  R[i*NHosts + j] = In * Sn * B[i * NHosts + j]/((double)(nS)) ;
	  
      }
      else
      R[i*NHosts + j] = In * Sn * B[i * NHosts + j]/((double)(nS)) ;
	  		
			sumR+= R[i*NHosts + j] ;
			
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
			//
			if (all_out==0)
      {
        //Iend[h_i] = I[h_i+(n*NHosts)] >= SN[h_i] ? SN[h_i] :  I[h_i+(n*NHosts)] ;
        Iend[h_i] = I[h_i+(n*NHosts)] ;
        I[h_i+(n*NHosts)] = I[h_i+(n*NHosts)] - SN[h_i] ;
        I[h_i+(n*NHosts)] = I[h_i+(n*NHosts)] > 0 ? I[h_i+(n*NHosts)] : 0 ;
        //S[h_i+(n*NHosts)] = 0 ;
        S[h_i+(n*NHosts)] = S[h_i+(n*NHosts)] ;
			  //T[n*4+3] = - I[h_i+((n-1)*NHosts)];
			  T[n*4+3] = - SN[h_i];
      }else{
        I[h_i+(n*NHosts)] = 0 ;
        S[h_i+(n*NHosts)] = 0 ;
        T[n*4+3] = - I[h_i+((n-1)*NHosts)];
			  Iend[h_i] = I[h_i+((n-1)*NHosts)]; 
      }
			
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
      In = I[(n-1)*NHosts + j] ;
      Sn = S[(n-1)*NHosts + i] ;
    	if (i!=j)
      {
        bnsize = bnsizes[j] ;
        In = (In >= bnsize) ? In : 0 ;
        Sn = (Sn >= bnsize) ? Sn : 0 ;
      }
      R[i*NHosts + j] = In * Sn * B[i * NHosts + j]/(double)(nS) ;
			
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
      bnsize = bnsizes[er] ;
      bn = (ec==er) ? 1 : bnsize ;
			I[n*NHosts + ec] = I[n*NHosts + ec] + bn;
			S[n*NHosts + ec] = S[n*NHosts + ec] - bn;
			T[n*4+3] = bn;
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
