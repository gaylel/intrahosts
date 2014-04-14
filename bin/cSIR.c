#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include "llists.h"
#include "phylo.h"
#include "leaf.h"
#include "beaglefuncs.h"


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
SEXP tree_reconstruct_with_partialll3(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS2, SEXP R_OBS, SEXP R_tr_old, SEXP R_mig2) ;
double ** calculateSIR_rates(double ** lambdao, double *B, int *PR_I, int *PR_S, int ha, int hb, int i, int kt, int NHosts, int NS, int TN, int *Anc, int *OBS2) ;
double calculatelambdaosum(double ** lambdao, int NHosts) ;
double ** Migmx_addnode(int mig_n, double **mig, double t, int ha, int hb, int Node) ;
SEXP draw_subtree_R(SEXP R_sirlist, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS2, SEXP R_OBS, SEXP R_tr_old_list, SEXP R_Ntips) ;
void smc_resample_R(SEXP R_ts, SEXP R_fs, SEXP R_p) ;
void smc_resample(int to_size, int from_size, double *p, int *s, int smode) ;  
double smc_logspace_add(double *logx, int xl) ;

struct smcparticle {
	phylo *tr ;		// current tree
	int sir_i ;	// index of sir trajectory to use
	double *mig;	// record of migratory nodes
	int mig_n ;		// number of migrations
} ;
typedef struct smcparticle smcparticle ;

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

struct smcinfo {
	smcparticle *sp ;
	int **sir1 ;
	double **sir2 ;
	int n ;
	
} ;

typedef struct smcinfo smcinfo ;

struct smcinfo2 {
	smcparticle *sp ;
	int ***sir1 ;
	double ***sir2 ;
	int *n ;
	double ll ;
	int smp ;
} ;

typedef struct smcinfo2 smcinfo2 ;

void smc_treereconstruct(int *PR_I, double *PR_T, int NHosts, int TN, int *SN, double *ST, int *OBS, int *OBS2, int has_old, phylo *tr_old3, double *mig3, int mig3_N, smcparticle* sp) ;
void smc_obs(int *SN, int NHosts, int tot, int* o) ;
SEXP smc_draw_prior_R(SEXP R_Np, SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN, SEXP R_bnprob, SEXP R_K) ;
hnode2  Nodes_addelement(int val, double t, hnode2 Nodes) ;
hnode2  Nodes_initialise(hnode2 Nodes) ;
int * Nodes_pick_pair(int *ch, hnode2 Nodesa, hnode2 Nodesb, int ha, int hb) ;
int * Nodes_pick_child(int *ch, hnode2 Nodesa) ;
hnode2 Nodes_deleteelement(int val, hnode2 Nodes) ;
int * Nodes_pick_sets(int *ov, hnode2 Nodesa1, hnode2 Nodesa2, hnode2 Nodesb1, hnode2 Nodesb2, int ha, int hb) ;
int * Nodes_pick_sets2(int *ov, hnode2 Nodesa1, hnode2 Nodesa2, hnode2 Nodesb1, hnode2 Nodesb2, int *NI, int ha, int hb, int has_old) ;
int * Nodes_pick_children(int *ch, int *ov, hnode2 Nodesa1, hnode2 Nodesa2, hnode2 Nodesb1, hnode2 Nodesb2, int ha, int hb) ;
double * Migmx_replacenode(int mig_n, int mig_N, double *mig, int oldnode, int newnode) ;
hnode2 Nodes_copy_all(hnode2 Nodes1, hnode2 Nodes2) ;
hnode2 Nodes_delete_all(hnode2 Nodes1) ;
hnode2 Nodes_deleteind(int ind, hnode2 Nodes) ;
void beagle_init_R(SEXP R_seqs, SEXP R_Nseqs, SEXP R_NSites, SEXP R_w, SEXP R_tr) ;
SEXP smc_draw_R(SEXP R_Np, SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN, SEXP R_bnprob, SEXP R_K, SEXP R_seqs, SEXP R_Nseqs, SEXP R_NSites, SEXP R_w, SEXP R_mu) ;
smcinfo2* smc_draw(int Np, int I0, int NS, int NHosts, double *B, double dr, double *ST, int *SN, double bnprob, int K, int NSeqs, int NSites, int** seqs, int* w, double mu) ;
int smc_trajcheck(int *SN, int *Iend, int NHosts) ;
void smc_freeparticle(smcparticle *sp) ;
SEXP smc_draw_prior_R(SEXP R_Np, SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN, SEXP R_bnprob, SEXP R_K) ;
smcinfo* smc_draw_prior(int Np, int I0, int NS, int NHosts, double *B, double dr, double *ST, int *SN, double bnprob, int K) ;

void beagle_init_R(SEXP R_seqs, SEXP R_Nseqs, SEXP R_NSites, SEXP R_w, SEXP R_tr)
{
	int NSeqs , NSites, *seq, i, j, **seqs, *w;
	double ll , mu=1;
	phylo * tr ;
	tr = R_to_phylo(R_tr) ;
	Rprintf("read in tree!\n") ;
	SEXP R_seq;
	NSeqs = INTEGER(coerceVector(R_Nseqs, INTSXP))[0] ;
	NSites = INTEGER(coerceVector(R_NSites, INTSXP))[0] ;
	seqs = calloc(NSeqs, sizeof(int*)) ;
	w = calloc(NSites, sizeof(int)) ;
	for (i=0 ; i<NSites ; i++)
	{
		w[i] =  INTEGER(coerceVector(R_w, INTSXP))[i] ;
	}
	for (i=0 ; i<NSeqs ; i++)
	{
		R_seq = VECTOR_ELT(R_seqs, i) ;
		seq = INTEGER(coerceVector(R_seq, INTSXP)) ;
	
		seqs[i] = calloc(NSites, sizeof(int)) ;
		for (j=0 ; j<NSites ; j++)
		{
			seqs[i][j] = seq[j] - 1;
			//Rprintf("%i ", seqs[i][j]) ;	
		}
		Rprintf("\n") ;
	}
	
	
	
	i=1 ;
	ll = beagle_init(NSeqs, NSites, seqs, w, tr->edge, tr->el, mu, tr->NNode) ;
	Rprintf("%8.4f\n", ll) ;	
	
	for (i=0 ; i<NSeqs ; i++)
	{
		free(seqs[i]) ;
	}
	free(seqs) ;
	free(w) ;
}
/*
smcparticle* smc_copyparticle(smcparticle *sp)
{
	smcparticle* sp2 = calloc(1, sizeof(smcparticle));
	int ntips = sp->NNode + 1 ;
	int i ;
	sp2 = phylo_create(ntips) ;
	sp2->Nedge=sp->Nedge ;
	for (i=0 ; i<
	tr->el=Calloc(Nedges,double) ;
	tr->NNode=ntips-1 ;
	tr->edge=Calloc(Nedges, int*) ;
	tr->nodelabel=Calloc(tr->NNode,int) ;
	for (m=0 ; m< Nedges ; m++)
	{
		tr->edge[m] = Calloc(2, int) ;
	}
	
	tr->tiplabel=Calloc(ntips,int) ;
	for (m=0 ; m< ntips ; m++)
	{
		tr->tiplabel[m]=m ;
	}
	
	sp2->sir_i = sp->sir_i ;
	sp2->mig_n = sp->mig_n ;
	for (i=0 ; i<(sp2->mig_n * 4) ; i++)
	{
		sp2->mig[i] = sp->mig[i] ;
	}
	return sp2 ;
}
*/
void smc_resample_R(SEXP R_ts, SEXP R_fs, SEXP R_p)
{
	int to_size = INTEGER(coerceVector(R_ts, INTSXP))[0] ;
	int from_size = INTEGER(coerceVector(R_fs, INTSXP))[0] ;
	double *p = REAL(coerceVector(R_p, REALSXP)) ;
	int *s = calloc(to_size, sizeof(int)) ;
	smc_resample(to_size, from_size, p, s, 1) ;
	
	int i ;
	for (i=0 ; i<to_size ; i++)
	{
		Rprintf("%i ", s[i]) ;
	}
	Rprintf("\n") ;
	free(s) ;
	
	 
}

void smc_resample(int to_size, int from_size, double *p, int *s, int smode)
{

	//1:  stratified sampling..
	//2:  systematic sampling
	int i, j ;
	double r , cp, r2 ; //cumsum of p;
	switch( smode ){
	case 1:
		j = 0 ;
		cp = p[0] ;
		for (i=0 ; i<to_size ; i++)
		{
			r = (unif_rand() + (double)i) / to_size ;
			//Rprintf("%8.4f %8.4f", r, cp) ;
			if (r <= cp)
			{
				s[i] = j ;
			}
			else
			{
				while (r>cp)
				{
				j++ ;
				cp+= p[j] ;
				
				}
				s[i] = j ;
			}
			//Rprintf(" %i\n", s[i]) ;
		}
	
		break ;
	case 2:
		j = 0 ;
		cp = p[0] ;
		r2 = unif_rand() ;
		for (i=0 ; i<to_size ; i++)
		{
			r = (r2 + (double)i) / to_size ;
			//Rprintf("%8.4f %8.4f", r, cp) ;
			if (r <= cp)
			{
				s[i] = j ;
			}
			else
			{
				while (r>cp)
				{
				j++ ;
				cp+= p[j] ;
				
				}
				s[i] = j ;
			}
			//Rprintf(" %i\n", s[i]) ;
		}
	
		break ;	
	}
	

}

double smc_logspace_add(double *logx, int xl)
{
	
	// returns log(sum(exp(logx)))
	double L ;
	int i ;
	if (xl == 1)
	{
		L = logx[0] ;
	}
	else
	{
		L = logspace_add(logx[0], logx[1]) ;
		for (i=2 ; i<xl ; i++)
		{
			L = logspace_add(L, logx[i]) ;
		}
	}
	return L ;
		
}

SEXP smc_draw_R(SEXP R_Np, SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN, SEXP R_bnprob, SEXP R_K, SEXP R_seqs, SEXP R_Nseqs, SEXP R_NSites, SEXP R_w, SEXP R_mu)
{
	// wrapper for smc_draw
	// draws a tree and trajectory based on an smc algorithm
	
	// input parameters
	
	int Np = INTEGER(coerceVector(R_Np, INTSXP))[0] ;		// number of SMC particles
	int I0 = INTEGER(coerceVector(R_I0, INTSXP))[0] ;		// initial number of infected particles
	int NS = INTEGER(coerceVector(R_NS, INTSXP))[0] ;		// number of particles in each host
	int NHosts = INTEGER(coerceVector(R_NHosts, INTSXP))[0] ; // number of hosts
	double *B = REAL(coerceVector(R_B, REALSXP)) ; 			// matrix of transmission coefficients
	double dr = REAL(coerceVector(R_dr, REALSXP))[0] ; 		// death rate
	double *ST = REAL(coerceVector(R_ST, REALSXP)) ;		// vector of sampling times
	int *SN = INTEGER(coerceVector(R_SN, INTSXP)) ;			// vector of viral population sizes
	double bn_prob = REAL(coerceVector(R_bnprob, REALSXP))[0] ;	//parameter for bottleneck size
	int K = INTEGER(coerceVector(R_K, INTSXP))[0] ;			// bottleneck size
	int NSeqs , NSites, *seq, i, j, **seqs, *w;
	double ll ;
	SEXP R_seq;
	NSeqs = INTEGER(coerceVector(R_Nseqs, INTSXP))[0] ;
	NSites = INTEGER(coerceVector(R_NSites, INTSXP))[0] ;
	seqs = calloc(NSeqs, sizeof(int*)) ;
	w = calloc(NSites, sizeof(int)) ;
	for (i=0 ; i<NSites ; i++)
	{
		w[i] =  INTEGER(coerceVector(R_w, INTSXP))[i] ;
	}
	for (i=0 ; i<NSeqs ; i++)
	{
		R_seq = VECTOR_ELT(R_seqs, i) ;
		seq = INTEGER(coerceVector(R_seq, INTSXP)) ;
	
		seqs[i] = calloc(NSites, sizeof(int)) ;
		for (j=0 ; j<NSites ; j++)
		{
			seqs[i][j] = seq[j] - 1;
			//Rprintf("%i ", seqs[i][j]) ;	
		}
		//Rprintf("\n") ;
	}
	double mu = REAL(coerceVector(R_mu, REALSXP))[0] ;
	
	//Rprintf("read in parameters\n") ;
	smcinfo2 * si ;
	si = smc_draw(Np, I0, NS, NHosts, B, dr, ST, SN, bn_prob, K, NSeqs, NSites, seqs, w, mu) ;
	
	
	
	
	/*int i ;
	SEXP R_sirlist, R_sir ;
	PROTECT(R_sirlist = allocVector(VECSXP ,Np)) ;
	for (i=0 ; i<Np ; i++)
	{
		R_sir = sample_cSIR_S_R(R_I0, R_NS, R_NHosts, R_B, R_dr, R_ST, R_SN, R_bnprob, R_K) ;
		SET_VECTOR_ELT(R_sirlist, i, R_sir) ;
	}
	UNPROTECT(1) ;*/
	// reconstruct trees
	
	SEXP R_list ;
	PROTECT(R_list = allocVector(VECSXP, 6)) ;
	SEXP R_tr , R_I, R_S, R_T, R_Iend, R_ll;
	int smp = si->smp, ti = si->sp[smp].sir_i, n_ti = si->n[ti];
	

	R_tr = phylo_to_R(si->sp[smp].tr) ;
	SET_VECTOR_ELT(R_list, 0, R_tr) ;
	PROTECT(R_I=allocMatrix(INTSXP,n_ti,NHosts)) ;
	PROTECT(R_S=allocMatrix(INTSXP,n_ti,NHosts)) ;
	PROTECT(R_T=allocMatrix(REALSXP,n_ti,4)) ;
	PROTECT(R_Iend=allocMatrix(INTSXP,1,NHosts)) ;
	PROTECT(R_ll = allocVector(REALSXP, 1)) ;
	REAL(R_ll)[0] = si->ll ;
	int * PR_I = INTEGER(R_I) ;
	int * PR_S = INTEGER(R_S) ;
	double * PR_T = REAL(R_T) ;
	int *PR_Iend = INTEGER(R_Iend) ;

	for (j=0; j< NHosts ; j++)
	{
		PR_Iend[j] = si->sir1[ti][2][j] ;
		for (i=0; i< n_ti; i++)
		{
		
			PR_I[i+j*n_ti]=si->sir1[ti][1][i*NHosts + j] ;
			PR_S[i+j*n_ti]=si->sir1[ti][0][i*NHosts + j] ;
		}
	}

	for (i=0; i<(n_ti) ; i++)
	{
		for (j=0; j<4; j++)
		{
			PR_T[i+j*n_ti] = si->sir2[ti][0][i*4+j] ;
		}
	}
	
	SET_VECTOR_ELT(R_list, 1, R_I) ;
	SET_VECTOR_ELT(R_list, 2, R_S) ;
	SET_VECTOR_ELT(R_list, 3, R_T) ;
	SET_VECTOR_ELT(R_list, 4, R_Iend) ;
	SET_VECTOR_ELT(R_list, 5, R_ll) ;
	UNPROTECT(6) ;	
	// free memory
	
	for (i=0 ; i<Np ; i++)
	{
		//Rprintf("free particle\n") ;
		//smc_freeparticle(&sp_old[i]) ;
		//smc_freeparticle(si[i].sp) ;
		Free(si->sir1[i][0]) ;
		Free(si->sir1[i][1]) ;
		Free(si->sir1[i][2]) ;
		free(si->sir1[i]) ;
		Free(si->sir2[i][0]) ;
		free(si->sir2[i]) ;
		if (i != smp)
		{
			smc_freeparticle(&si->sp[i]) ;
		}
		else
		{
		if (si->sp[i].mig)
		{	
			free(si->sp[i].mig) ;
		}
		}
	}
	free(si->sir1) ;
	free(si->sir2) ;
	free(si->sp) ;
	free(si->n) ;
	free(si) ;
	
	
	
	
	// free memory
	free(w) ;
	for (j=0 ; j<NSeqs ; j++)
	{
		free(seqs[j]) ;
	}
	free(seqs) ;
	
	return R_list ;
}


SEXP smc_draw_prior_R(SEXP R_Np, SEXP R_I0, SEXP R_NS, SEXP R_NHosts, SEXP R_B, SEXP R_dr, SEXP R_ST, SEXP R_SN, SEXP R_bnprob, SEXP R_K)
{
	// wrapper for smc_draw
	// draws a tree and trajectory based on an smc algorithm
	
	// input parameters
	
	int Np = INTEGER(coerceVector(R_Np, INTSXP))[0] ;		// number of SMC particles
	int I0 = INTEGER(coerceVector(R_I0, INTSXP))[0] ;		// initial number of infected particles
	int NS = INTEGER(coerceVector(R_NS, INTSXP))[0] ;		// number of particles in each host
	int NHosts = INTEGER(coerceVector(R_NHosts, INTSXP))[0] ; // number of hosts
	double *B = REAL(coerceVector(R_B, REALSXP)) ; 			// matrix of transmission coefficients
	double dr = REAL(coerceVector(R_dr, REALSXP))[0] ; 		// death rate
	double *ST = REAL(coerceVector(R_ST, REALSXP)) ;		// vector of sampling times
	int *SN = INTEGER(coerceVector(R_SN, INTSXP)) ;			// vector of viral population sizes
	double bn_prob = REAL(coerceVector(R_bnprob, REALSXP))[0] ;	//parameter for bottleneck size
	int K = INTEGER(coerceVector(R_K, INTSXP))[0] ;			// bottleneck size
	
	double ll ;
	
	Rprintf("read in parameters\n") ;
	smcinfo * si;
	
	si = smc_draw_prior(Np, I0, NS, NHosts, B, dr, ST, SN, bn_prob, K) ;
	
	
	
	
	/*int i ;
	SEXP R_sirlist, R_sir ;
	PROTECT(R_sirlist = allocVector(VECSXP ,Np)) ;
	for (i=0 ; i<Np ; i++)
	{
		R_sir = sample_cSIR_S_R(R_I0, R_NS, R_NHosts, R_B, R_dr, R_ST, R_SN, R_bnprob, R_K) ;
		SET_VECTOR_ELT(R_sirlist, i, R_sir) ;
	}
	UNPROTECT(1) ;*/
	// reconstruct trees
	
	
	SEXP R_list ;
	PROTECT(R_list = allocVector(VECSXP, 5)) ;
	SEXP R_tr , R_I, R_S, R_T, R_Iend;
	R_tr = phylo_to_R(si->sp[0].tr) ;
	SET_VECTOR_ELT(R_list, 0, R_tr) ;
	PROTECT(R_I=allocMatrix(INTSXP,si[0].n,NHosts)) ;
	PROTECT(R_S=allocMatrix(INTSXP,si[0].n,NHosts)) ;
	PROTECT(R_T=allocMatrix(REALSXP,si[0].n,4)) ;
	PROTECT(R_Iend=allocMatrix(INTSXP,1,NHosts)) ;
	int * PR_I = INTEGER(R_I) ;
	int * PR_S = INTEGER(R_S) ;
	double * PR_T = REAL(R_T) ;
	int *PR_Iend = INTEGER(R_Iend) ;
	int i, j ;
	
	for (j=0; j< NHosts ; j++)
	{
		PR_Iend[j] = si[0].sir1[2][j] ;
		for (i=0; i<si[0].n; i++)
		{
		
			PR_I[i+j*si[0].n]=si[0].sir1[1][i*NHosts + j] ;
			PR_S[i+j*si[0].n]=si[0].sir1[0][i*NHosts + j] ;
		}
	}
	
	for (i=0; i<(si[0].n) ; i++)
	{
		for (j=0; j<4; j++)
		{
			PR_T[i+j*si[0].n] = si[0].sir2[0][i*4+j] ;
		}
	}
	
	
	SET_VECTOR_ELT(R_list, 1, R_I) ;
	SET_VECTOR_ELT(R_list, 2, R_S) ;
	SET_VECTOR_ELT(R_list, 3, R_T) ;
	SET_VECTOR_ELT(R_list, 4, R_Iend) ;
	UNPROTECT(5) ;	
	// free memory
	
	
	
	for (i=0 ; i<Np ; i++)
	{
		//Rprintf("free particle\n") ;
		//smc_freeparticle(&sp_old[i]) ;
		//smc_freeparticle(si[i].sp) ;
		Free(si[i].sir1[0]) ;
		Free(si[i].sir1[1]) ;
		Free(si[i].sir1[2]) ;
		free(si[i].sir1) ;
		Free(si[i].sir2[0]) ;
		free(si[i].sir2) ;
		if (si[i].sp->mig)
		{
			free(si[i].sp->mig) ;
			si[i].sp->mig = NULL ;
		}
		free(si[i].sp) ;
		//free(si[i].n) ;
	      
	}
	
		
	free(si) ;	
	
	return R_list ;
	
	//SEXP R_A ;
	//PROTECT(R_A = allocVector(INTSXP, 1)) ;
	//R_A[0] = 1 ;	
	
	//return R_A ;
}

int smc_trajcheck(int *SN, int *Iend, int NHosts)
{
	int is_acc = 1 ;
	int i ;
	for (i=0 ; i<NHosts ; i++)
	{
		if (SN[i] > Iend[i])
		{
			is_acc = 0 ;
			break ;
		}
	}
	return is_acc ;

}


void smc_obs(int *SN, int NHosts, int tot, int* o)
{
	int i=0 ;
	while ((i < NHosts) & (tot>0))
	{
		o[i] = (tot > SN[i]) ? SN[i] : tot ;
		tot = tot - SN[i] ;
		i++ ;
	}
}

void smc_freeparticle(smcparticle *sp)
{
	//Rprintf("freeing tree at %p\n", sp->tr) ;
	
	if (sp->tr)
	phylo_free(sp->tr) ;
	
	sp->tr = 0 ;
	if (sp->mig)
	{
	//	Rprintf("freeing mig at %p\n", sp->mig) ;
		free(sp->mig) ;
		sp->mig = NULL ;
	}
	//free(sp) ;
}

smcinfo* smc_draw_prior(int Np, int I0, int NS, int NHosts, double *B, double dr, double *ST, int *SN, double bnprob, int K)
{
	// draws a tree and trajectory based on an smc algorithm
	smcinfo * si  = calloc(Np, sizeof(smcinfo)) ;
	
	// setup helper variables 
	int i, j ;
	double n_eff , n_eff_th = 1 * Np ;
	int *n = calloc(Np, sizeof(int)), Nt=Np ;		//vector of trajectory lengths
	for (i=0 ; i<Np ; i++)
	{
		n[i] = 0 ;
	}
	
	// set up particle systems
	smcparticle * sp = calloc(Np, sizeof(smcparticle)) ;
	
	for (i=0 ; i<Np ; i++)
	{
		sp[i].mig = NULL ;
		
	}
	
	int ***p_RVAL1 = calloc(Np, sizeof(int**)) ;
	double ***p_RVAL2 = calloc(Np, sizeof(double**)) ;
	
	double Ltot = 0;
	int *traj_acc = calloc(Np, sizeof(int)) ;	//indicator vector of whether trajectory is accepted
	int n_acc=0;
	int ti, k ;
	// draw Np trajectories
	for (i=0 ; i<Np ; i++)
	{
		p_RVAL1[i] = calloc(3, sizeof(int*)) ;
		p_RVAL2[i] = calloc(1, sizeof(double*)) ;
		traj_acc[i] = 0 ;
		while (traj_acc[i] == 0)
		{
			cSIR_iters_ST(&n[i], I0, NS, NHosts, B, dr, p_RVAL1[i], p_RVAL2[i], ST, SN, bnprob, K) ;
			traj_acc[i] = smc_trajcheck(SN, p_RVAL1[i][2], NHosts) ;
			if (traj_acc[i] == 0)
			{
				Free(p_RVAL1[i][0]) ;
				Free(p_RVAL1[i][1]) ;
				Free(p_RVAL1[i][2]) ;
				Free(p_RVAL2[i][0]) ;
			}
		}
		Rprintf("%i\t%i\t%i\n", n[i], traj_acc[i], i) ;
		//Rprintf("%8.4f\t%8.4f\t%8.4f\t%8.4f\n", p_RVAL2[i][0][4], p_RVAL2[i][0][5], p_RVAL2[i][0][6], p_RVAL2[i][0][7]) ;
		n_acc += traj_acc[i] ;
	}
	
	int *OBS = calloc(NHosts, sizeof(int)) ;
	int *OBS2 = calloc(NHosts, sizeof(int)) ;
	double SNsum = 0 ;
	
	for (i=0 ; i< NHosts ; i++)
	{
		SNsum += SN[i] ;
	}
	
	smc_obs(SN, NHosts, SNsum, OBS) ;
	smc_obs(SN, NHosts, 0, OBS2) ;
	
	for (i=0 ; i<Np ; i++)
	{
		smc_treereconstruct(p_RVAL1[i][1], p_RVAL2[i][0], NHosts, n[i], SN, ST, OBS, OBS2, 0, 0, 0, 0, &sp[i]) ;
		si[i].sp = &sp[i] ;
		si[i].sir1 = p_RVAL1[i] ;
		si[i].sir2 = p_RVAL2[i] ;
		si[i].n = n[i] ;
	}
	
	
	free(p_RVAL1) ;
	free(p_RVAL2) ;
	free(traj_acc) ;
	free(n) ;
	free(OBS) ;
	free(OBS2) ;
	Rprintf("got here\n") ;
	// free particles
	
	
	
	return(si) ;
}

smcinfo2 * smc_draw(int Np, int I0, int NS, int NHosts, double *B, double dr, double *ST, int *SN, double bnprob, int K, int NSeqs, int NSites, int** seqs, int* w, double mu)
{
	// draws a tree and trajectory based on an smc algorithm
	
	// setup helper variables 
	int i, j ;
	double n_eff , n_eff_th = 1 * Np ;
	int *n = calloc(Np, sizeof(int)), Nt=Np ;		//vector of trajectory lengths
	for (i=0 ; i<Np ; i++)
	{
		n[i] = 0 ;
	}
	
	// set up particle systems
	smcparticle * sp = calloc(Np, sizeof(smcparticle)) ;
	smcparticle * sp_old = calloc(Np, sizeof(smcparticle)) ;
	for (i=0 ; i<Np ; i++)
	{
		sp_old[i].mig = NULL ;
		sp[i].mig = NULL ;
	}
	
	int ***p_RVAL1 = calloc(Np, sizeof(int**)) ;
	double ***p_RVAL2 = calloc(Np, sizeof(double**)) ;
	double *p = calloc(Np, sizeof(double)) ;
	
	double Ltot = 0;
	int *s = calloc(Np, sizeof(int)) ;
	int *traj_acc = calloc(Np, sizeof(int)) ;	//indicator vector of whether trajectory is accepted
	int n_acc=0;
	int ti, k ;
	// draw Np trajectories
	for (i=0 ; i<Nt ; i++)
	{
		p_RVAL1[i] = calloc(3, sizeof(int*)) ;
		p_RVAL2[i] = calloc(1, sizeof(double*)) ;
		traj_acc[i] = 0 ;
		while (traj_acc[i] == 0)
		{
			cSIR_iters_ST(&n[i], I0, NS, NHosts, B, dr, p_RVAL1[i], p_RVAL2[i], ST, SN, bnprob, K) ;
			traj_acc[i] = smc_trajcheck(SN, p_RVAL1[i][2], NHosts) ;
			if (traj_acc[i] == 0)
			{
				Free(p_RVAL1[i][0]) ;
				Free(p_RVAL1[i][1]) ;
				Free(p_RVAL1[i][2]) ;
				Free(p_RVAL2[i][0]) ;
			}
		}
		//Rprintf("%i\t%i\t%i\n", n[i], traj_acc[i], i) ;
		//Rprintf("%8.4f\t%8.4f\t%8.4f\t%8.4f\n", p_RVAL2[i][0][4], p_RVAL2[i][0][5], p_RVAL2[i][0][6], p_RVAL2[i][0][7]) ;
		n_acc += traj_acc[i] ;
	}
	Rprintf("%i\n", n_acc) ;
	for (i=0 ; i<Np ; i++)
	{
		p[i] = (traj_acc[i]==1) ? (double)1 / n_acc : 0 ;
	}
	smc_resample(Np, Nt, p, s, 2) ;
	
	int *OBS = calloc(NHosts, sizeof(int)) ;
	int *OBS2 = calloc(NHosts, sizeof(int)) ;
	smc_obs(SN, NHosts, 2, OBS) ;
	smc_obs(SN, NHosts, 0, OBS2) ;
	
	for (i=0 ; i<Np ; i++)
	{
		ti = s[i] ;
		
		sp[i].sir_i = ti ;
		//Rprintf("parent = %i, %i  %i\n", s[i], ti, i) ;
		smc_treereconstruct(p_RVAL1[ti][1], p_RVAL2[ti][0], NHosts, n[ti], SN, ST, OBS, OBS2, 0, 0, 0, 0, &sp[i]) ;
		p[i] = beagle_init(2, NSites, seqs, w, sp[i].tr->edge, sp[i].tr->el, mu, 1);					
	}
	
	Ltot = smc_logspace_add(p, Np) ;
	
	n_eff = 0 ;
	for (i=0 ; i<Np ; i++)
	{
		p[i] = exp(p[i] -Ltot) ;
		//Rprintf("%8.4f ", p[i]) ;
		n_eff += p[i]*p[i] ;
	}
	n_eff = (double )1/n_eff ;
	//Rprintf("%8.4f\n", n_eff) ;
	if (n_eff < n_eff_th)
	{
		smc_resample(Np, Np, p, s, 2) ;
	}
	else
	{
		for (i=0 ; i<Np ; i++)
		{
			s[i] = i ;
		}
	}
	
	for (j=2 ; j<(NSeqs) ; j++)
	{
	//	Rprintf("j = %i\n", j) ;
		smc_obs(SN, NHosts, j+1, OBS) ;
		smc_obs(SN, NHosts, j, OBS2) ;
		for (i=0 ; i<NHosts ; i++)
		{
			OBS[i] = OBS[i] - OBS2[i] ; 
		}
	
		// update particles
		for (i=0 ; i<Np ; i++)
		{
			sp_old[i] = sp[i] ;
			
			sp[i].mig = NULL ;
			sp[i].tr = NULL ;
			//phylo_free(sp[i].tr) ;
			//sp[i].tr = 0 ;
		}	
	
		for (i=0 ; i<Np ; i++)
		{
			ti = sp_old[s[i]].sir_i ;
			
			//smc_freeparticle(&sp[i]) ;
			//Rprintf("freeing particle at %p\n", &sp[i]) ;
			sp[i].sir_i = ti ;
			//Rprintf("parent = %i, old tree at %p, old mig at   %p\n", s[i], sp_old[s[i]].tr, sp_old[s[i]].mig) ;
			smc_treereconstruct(p_RVAL1[ti][1], p_RVAL2[ti][0], NHosts, n[ti], SN, ST, OBS, OBS2, 1, sp_old[s[i]].tr, sp_old[s[i]].mig, sp_old[s[i]].mig_n, &sp[i]) ;
			//Rprintf("recon tree\n") ;
			for (k=0 ; k<sp[i].tr->Nedge ; k++)
			{
				//Rprintf("%i %i\n", sp[i].tr->edge[k][0], sp[i].tr->edge[k][1]) ;
			}
			
			p[i] = beagle_init(j+1, NSites, seqs, w, sp[i].tr->edge, sp[i].tr->el, mu, j);
		}
		
		for (i=0 ;  i<Np ;i++)
		{
			
			//Rprintf("freeing particle %i at %p, with tr at %p, mig at %p\n", i, &sp_old[i], sp_old[i].tr, sp_old[i].mig) ;
			smc_freeparticle(&sp_old[i]) ;
		//	Rprintf("new particle %i at %p, with tr at %p, mig at %p\n", i, &sp[i], sp[i].tr, sp[i].mig) ;
		}
		Ltot = smc_logspace_add(p, Np) ;
		n_eff = 0 ;
		for (i=0 ; i<Np ; i++)
		{
			p[i] = exp(p[i] -Ltot) ;
			//Rprintf("%i ", sp[i].sir_i) ;
			n_eff += p[i]*p[i] ;
		}
		//Rprintf("\n") ;
		n_eff = (double )1/n_eff ;
		//Rprintf("%8.4f\n", n_eff) ;
		if (n_eff < n_eff_th)
		{
			smc_resample(Np, Np, p, s, 1) ;
		}
		else
		{
			for (i=0 ; i<Np ; i++)
			{
				s[i] = i ;
			}
		}
	}
	
	
	// select particle and trajectory
	rmultinom(1, p, Np, s) ;
	
	
	for (j=0 ; j<Np ; j++)
	{
		if (s[j]==1)
		{
			break ;		
		}
	}
	
	
	
	smcinfo2 * si  = calloc(1, sizeof(smcinfo2)) ;
	si[0].sp = sp ;
	si[0].sir1 = p_RVAL1 ;
	si[0].sir2 = p_RVAL2 ;
	si[0].n = n ;
	si[0].ll = Ltot ;
	si[0].smp = j ;
	/*
	for (i=0 ; i<Nt ; i++)
	{
		//Rprintf("%i ", sp[i].sir_i) ;
		Free(p_RVAL1[i][0]) ;
		Free(p_RVAL1[i][1]) ;
		Free(p_RVAL1[i][2]) ;
		free(p_RVAL1[i]) ;
		Free(p_RVAL2[i][0]) ;
		free(p_RVAL2[i]) ;
	}
	//Rprintf("\n") ;
		
	free(p_RVAL1) ;
	free(p_RVAL2) ;
	free(n) ;
	for (i=0 ; i<Np ; i++)
	{
		//Rprintf("free particle\n") ;
		//smc_freeparticle(&sp_old[i]) ;
		
		smc_freeparticle(&sp[i]) ;
	}
	free(sp) ;
	
	*/
	free(traj_acc) ;
	free(OBS) ;
	free(OBS2) ;
	free(p) ;
	free(s) ;
	
	// free particles
	
	free(sp_old) ;
	return(si) ;
	
	
	
	
}

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

hnode2 Nodes_deleteind(int ind, hnode2 Nodes)
{
	int val;
	val = llist_get_el_i(ind, Nodes.n) ;
	Nodes.n = llist_delete_el_i(ind, Nodes.n) ;
	Nodes.t = llist_delete_el_d(ind, Nodes.t) ;
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
	return Nodes ;
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

int * Nodes_pick_child(int *ch, hnode2 Nodesa)
{
	// sampling without replacement
	int r1 ;
	r1 = (Nodesa.N > 1) ? rand() % Nodesa.N : 0;
	ch[0] = llist_get_el_i(r1, Nodesa.n) ;
	ch[1] = r1 ;
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
		if (ha!=hb)
		{
			if (rbinom(1, (double) Nodesb2.N/(Nodesb2.N + Nodesb1.N)) == 1)
			{
				ov[1] = 1 ; // 0,1
			}
		}
		else
		{
			if (rbinom(1, (double) Nodesa2.N/(Nodesa2.N + Nodesa1.N - 1)) == 1)
			{
				ov[1] = 1 ; // 0,1
			}
		}
	}
	// 0,0
	return ov ;
					  	
}


int * Nodes_pick_sets2(int *ov, hnode2 Nodesa1, hnode2 Nodesa2, hnode2 Nodesb1, hnode2 Nodesb2, int *NI, int ha, int hb, int has_old)
{
	ov[0] = 0;
	ov[1] = 0;
	double r1, r2 ;
	double pa[3], pb[3] ;
	int tota, totb ;
	int na[3], nb[3] ;
	int i, j ;
	tota = Nodesa1.N + Nodesa2.N + NI[ha] ;
	totb = (ha==hb) ? tota : Nodesb1.N + Nodesb2.N + NI[hb] ;
	
	
	// probability of node classes in hb
	
	pb[0] = (double) NI[hb]  ;	// non-ancestors
	pb[1] = (double) Nodesb1.N  ;	// new ancestors
	pb[2] = (double) Nodesb2.N  ; // old ancestors
	
	if ((ha!=hb) & (has_old==1))
	{
		totb -= pb[2] ;
		pb[2] = 0;
	}
	
	for (i=0 ; i<3 ; i++)
	{
		pb[i] = pb[i]/totb ;
	}
	
	rmultinom(1, pb, 3, nb) ;
	for (i=0 ; i<3 ; i++)
	{
		if (nb[i]==1)
		{
			ov[1] = i ;		
		}
	}
	
	// probability of node classes in ha
	pa[0] = (double) NI[ha] ;	// non-ancestors
	pa[1] = (double) Nodesa1.N  ;	// new ancestors
	pa[2] = (double) Nodesa2.N ; // old ancestors
	
	if (ha==hb)
	{
		pa[ov[1]]-- ;
		tota-- ;
	}
	
	if (ov[1]==2)
	{
		tota-=pa[2] ; 
		pa[2] = 0 ;
	}
	
	for (i=0 ; i<3 ; i++)
	{
		pa[i] = pa[i]/tota ;
	}
	
	rmultinom(1, pa, 3, na) ;
	for (i=0 ; i<3 ; i++)
	{
		if (na[i]==1)
		{
			ov[0] = i ;		
		}
	}
	
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
        	ch = Nodes_pick_pair(ch, Nodesa1, Nodesb2, ha, hb+1) ;
        
        }
    }else
        {
            ch = Nodes_pick_pair(ch, Nodesa2, Nodesb1, ha, hb+1) ;
    }
    
    return ch ;
}

hnode2 Nodes_copy_all(hnode2 Nodes1, hnode2 Nodes2)
{
	// copy nodes from Nodes2 to Nodes1
	int N_to_copy = Nodes2.N ;
	int i ;
	int val_i ;
	double val_d ;
	
	
	for (i=0 ; i<N_to_copy ; i++)
	{
		val_i = llist_get_el_i(i, Nodes2.n) ;
		val_d = llist_get_el_d(i, Nodes2.t) ;
		Nodes1 = Nodes_addelement(val_i, val_d, Nodes1) ;

	}
	return Nodes1 ;
}

hnode2 Nodes_delete_all(hnode2 Nodes1)
{
	int N_to_del = Nodes1.N ;
	int i ;
	int val_i ;
	for (i=0 ; i<N_to_del ; i++)
	{
		Nodes1 = Nodes_deleteind(0, Nodes1) ;
	}
	return Nodes1 ;
}

double ** Migmx_addnode(int mig_n, double **mig, double t, int ha, int hb, int Node)
{
	
		//Rprintf("%i\n", mig_n) ;
		if (mig_n == 1)				
		{
			mig = Calloc(mig_n, double *) ;
			mig[0] = Calloc(4, double) ;
					
		}
		else
		{
			mig = Realloc(mig, mig_n, double *) ;
			mig[mig_n - 1] = Calloc(4, double) ;
		}
		mig[mig_n - 1][0] = t ; 
		mig[mig_n - 1][1] = ha ; 
		mig[mig_n - 1][2] = hb ; 
		mig[mig_n - 1][3] =  (double) Node ;
		//Rprintf("%8.4f %8.4f %8.4f %8.4f\n", mig[mig_n - 1][0], mig[mig_n - 1][1], mig[mig_n - 1][2], mig[mig_n - 1][3]) ;
	
	return mig ;
}

double * Migmx_replacenode(int mig_n, int mig_N, double *mig, int oldnode, int newnode)
{
	// n : index of edge row
	int j ;
	
		for (j=mig_n ; j<mig_N ; j++)
		{
		if ((int) mig[j + 3*mig_N] == oldnode)
		{
			mig[j + 3*mig_N] = newnode ;
		}
		}
	return mig ;
}




SEXP draw_subtree_R(SEXP R_sirlist, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS2, SEXP R_OBS, SEXP R_tr_old_list, SEXP R_Ntips)
{
	// read list of sir trajectories
	int Np = length(R_sirlist) ;
	int Ntips = INTEGER(coerceVector(R_Ntips, INTSXP))[0] ;
	int i ;
	SEXP R_trlist, R_tr, R_tr_old ;
	PROTECT(R_trlist = allocVector(VECSXP ,Np)) ;
	
	if (Ntips == 2)
	{
		for (i=0 ; i<Np ; i++)
		{
			R_tr = tree_reconstruct_with_partialll3(VECTOR_ELT(R_sirlist, i), R_NHosts, R_dat, R_B, R_NS, R_BN, R_OBS2, R_OBS, 0, 0) ;
			SET_VECTOR_ELT(R_trlist, i, R_tr) ;
		}
	}
	else
	{
		for (i=0 ; i<Np ; i++)
		{
			R_tr_old = VECTOR_ELT(R_tr_old_list, i) ;
			R_tr = tree_reconstruct_with_partialll3(VECTOR_ELT(R_sirlist, i), R_NHosts, R_dat, R_B, R_NS, R_BN, R_OBS2, R_OBS, VECTOR_ELT(R_tr_old, 0), VECTOR_ELT(R_tr_old, 3)) ;
			SET_VECTOR_ELT(R_trlist, i, R_tr) ;
		}
	}
	UNPROTECT(1) ;
	return R_trlist ;
}




void smc_treereconstruct(int *PR_I, double *PR_T, int NHosts, int TN, int *SN, double *ST, int *OBS, int *OBS2, int has_old, phylo *tr_old3, double *mig3, int mig3_N, smcparticle* sp)
{
	int ce ;
	int ce1 ;
	double ll_i, small_val=1e-10 , ttmp2 ;
	int i, bn;
	
	// initialise helper variables
	
	double ** bt_old ;	// branching times of old tree
	int bt_n ;			// number of branching events of old tree
	int mig2_N = mig3_N;
	double *mig2 ; 
	phylo *tr_old ;
	if (has_old==1)
	{
		//Rprintf("has_old==1\n") ;
		// calculate branching times
		bt_old = phylo_bt(tr_old3, ST, OBS2, NHosts) ;
		bt_n = tr_old3->NNode - 1 ;
		mig2 = calloc(mig2_N * 4, sizeof(double)) ;
		for (i=0 ; i<(mig2_N * 4) ; i++)
		{
			mig2[i] = mig3[i] ;
		}
		// print out others
		
		//Rprintf("\nmig2_N=%i\n", mig2_N) ;
		tr_old = phylo_dup(tr_old3) ;
		//Rprintf("mig2 at %p\n", mig2) ;
		for (i=0 ; i<4*mig2_N ; i++)
		{
			//Rprintf("%8.4f ", mig2[i]) ;
		}
		//
		//Rprintf("\n\n") ;
		//Rprintf("\nmig2_N=%i\n", mig2_N) ;
		
		for (i=0 ; i<tr_old->Nedge ; i++)
		{
		 //Rprintf("%i %i\n", tr_old->edge[i][0], tr_old->edge[i][1]) ;
		}
		//Rprintf("read in old\n") ;
	}
	int mig2_n=0 ;
	
	
		
	/*************** INTERNAL VARIABLES FOR ALGORITHM ***********************************/
	
	int j, k=0, ha, hb, oa, ob, Ntrans ;
	double ***PR_T2, ttmp ;
	PR_T2 = Calloc(4,double**) ;	//reshape(PR_T)
	for (j=0 ; j<4 ; j++)
	{
		PR_T2[j] = Calloc(TN,double*) ;
		
	}
	
	for (i=0 ; i<TN ; i++)
	{
		for (j=0 ; j<4 ; j++)
		{
			PR_T2[j][i] = &PR_T[k++];
		}
	}
	
	int *Tend = Calloc(NHosts,int) ; // times at which hosts are sampled
	for (i=0 ; i< TN ; i++)
	{
		ha = (int) *PR_T2[1][i] ;
		hb = (int) *PR_T2[2][i] ;
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
		ntips+= OBS[i] + OBS2[i] ;
		lambdao[i] = calloc(NHosts, sizeof(double)) ;
	}
	hnode2 * Nodes = Calloc(NHosts,hnode2) ; 
	hnode2 * Nodes2 = Calloc(NHosts, hnode2) ;
	hnode2 * tNodes = Calloc(1,hnode2) ; 
	hnode2 * tNodes2 = Calloc(1, hnode2) ;
	int *NI = Calloc(NHosts, int) ;
	
	
	for (i=0 ; i<NHosts ; i++)
	{
		Nodes[i] = Nodes_initialise(Nodes[i]) ;
		Nodes2[i] = Nodes_initialise(Nodes2[i]) ;
	}
	
	// initialise temporary nodes
	*tNodes = Nodes_initialise(*tNodes) ;
	*tNodes2 = Nodes_initialise(*tNodes2) ;
	
	phylo *tr = phylo_create(ntips) ;
	
	int *ch = Calloc(4, int) ;		// vector of child nodes and indices
	int * ov = Calloc(2, int) ;
	
	int ei = tr->Nedge - 1 ;		// index of edges
	int intNode = 2*tr->NNode  ;	// index of internal node
	int *Anc = calloc(NHosts, sizeof(int)) ;
	
	// matrix to store the number of migrating lineages
	double ** mig ;
	int mig_n = 1 ;
	int NI_less ;
	int v1 ;
	int st ;
	int st2 ;
	int kt ;
	
  	j = 0;
  	// maybe still some infected left at end:
  
	for (i=0 ; i<NHosts ; i++)
	{
		//SNsum[i] = SN[i] + j ;
    	//j+= SN[i] ;
    	v1 = OBS[i] + OBS2[i] ;
    	SNsum[i] = v1 + j ;
    	j+= v1 ;
    	//NI[i] = PR_I[TN * (i+1) - 1] ;
    	NI[i] = PR_I[i + (TN - 1) * NHosts] ;
		Anc[i] = 0 ;
	}

	for (i=TN-1 ; i>-1 ; i--)
	{
		//Rprintf("%i\n", i) ;
		ll_i = 0 ;
		ce = 0 ;
		ce1 = 0 ;
    // for each event get hosts
		ha = (int) *PR_T2[1][i] ;
		hb = (int) *PR_T2[2][i] ;
		
    
    // if terminal event
   		if (i==Tend[ha])
    	{
    		//Rprintf("Terminal event %i\n", i) ;
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
		  	
		  	NI[ha] += -*PR_T2[3][i] - (OBS[ha] + OBS2[ha]) ;
		  	Anc[ha]+=OBS[ha] + OBS2[ha];  
		  	
		}
    	else
		{
			if ((i > 0) & (intNode >= ntips))
			{
				// calculate rates for each host
				kt = *PR_T2[3][i] ;	
				NI_less = 0 ;
					
				if ((int) *PR_T2[3][i] == -1 )
				{
					// death but birth in reverse.
					// can calculate the probability of no coalescence has occurred
					NI[ha]++ ;
				} 
				
				//if ((int) *PR_T2[3][i] == 1 )
				//{
					// 1 coalescent event within host
				//	ov = Nodes_pick_sets(ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], ha, hb) ;

				//}
				//Rprintf("%i\n", bt_n) ;
				if ((int) *PR_T2[3][i] > 0 )
				{
					Ntrans = (int)*PR_T2[3][i] ;
					// if coincides with coalescent interval among old observed tree
					if (has_old==1) //oa,ob=11
					{			
						if (bt_n >= 0)
						{
						while ((bt_old[bt_n][0]-small_val) > *PR_T2[0][i-1])
						{
							ch[0] = tr_old->edge[bt_n*2 + 1][1] ;
        					ch[1] = tr_old->edge[bt_n*2][1] ;
        					/*
							Rprintf("Existing coalescent event between %i from %i and  %i from %i at %8.4f, last bt = %8.4f, bt=%8.4f\n", ch[0], ha, ch[1], hb, *PR_T2[0][i], *PR_T2[0][i-1], bt_old[bt_n][0] ) ;
							llist_print(Nodes2[0].n) ;
							Rprintf(" ha2\n") ;
							llist_print(Nodes2[1].n) ;
							Rprintf(" hb2\n") ;
							Rprintf("Nodes in 0=\n") ;
        				    llist_print(Nodes[0].n) ;
        				    Rprintf("\nNodes in 1=\n") ;
        				    llist_print(Nodes[1].n) ;
							*/
							// create coalescent event in new tree
        					
        					ch[2] = llist_get_ind_i(ch[0], Nodes2[ha].n) ;
        					ch[3] = llist_get_ind_i(ch[1], Nodes2[hb].n) ;
        					//Rprintf("birth %i %i %i %i %i\n", ch[0], ch[1], ch[2], ch[3], ei) ;
        					// update new tree
        					tr = phylo_addnewedge(tr, ei, intNode, ch[0], llist_get_el_d(ch[2], Nodes2[ha].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr = phylo_addnewedge(tr, ei, intNode, ch[1], llist_get_el_d(ch[3], Nodes2[hb].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr->nodelabel[intNode-ntips] = ha ;
							//Rprintf("added edges\n") ;
							
							
							//Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
							//Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
							
							Nodes2[ha] = Nodes_deleteind(ch[2], Nodes2[ha]) ;
							Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
							
							
							// add later
							
							*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
							//Nodes2[ha] = Nodes_addelement(intNode, *PR_T2[0][i], Nodes2[ha]) ;
							
							if (bt_n > 0)
								tr_old = phylo_replacenode(tr_old, (bt_n*2)-1, tr_old->edge[bt_n*2][0], intNode) ;
							mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, tr_old->edge[bt_n*2][0], intNode) ;
							//Rprintf("Replacing %i with %i, mig2_n = %i, bt_n =%i\n", tr_old->edge[bt_n*2][0], intNode, mig2_n, bt_n ) ;
							
							intNode -- ;
							Anc[hb] -- ;
							//OBS2[hb]-- ;
							//ce++ ;
        					//Rprintf("%i\n", ei) ;
        					bt_n -- ;	
							Ntrans -- ;
						
							
						/*
							Rprintf("old tree\n") ;
							if (sumOBS > 0)
							{
								for (j=0 ; j<tr_old->Nedge ; j++)
								{
									//Rprintf("%i %i\n", tr_old->edge[j][0], tr_old->edge[j][1]) ;
								}
								
							}
							Rprintf("new tree\n") ;
							for (j=0 ; j<tr->Nedge ; j++)
							{
								Rprintf("%i %i\n", tr->edge[j][0], tr->edge[j][1]) ;
							}
							
							*/
							if (bt_n < 0) break ;
							
        				}
        				// look at migrations
        				}
        				if (mig2_n < mig2_N)
        				{
        					while ((mig2[mig2_n]-small_val) > *PR_T2[0][i-1])
        					{
        						//Rprintf("Old migration, node=%i, ha=%i, hb=%i\n", (int) mig2[mig2_n + 3*mig2_N], ha, hb) ;
        						
        						j = llist_get_ind_i(mig2[mig2_n + 3*mig2_N], Nodes2[hb].n) ;
        						ttmp = llist_get_el_d(j, Nodes2[hb].t) ;
						  		
								//Nodes2[hb] = Nodes_deleteelement(mig2[mig2_n + 3*mig2_N], Nodes2[hb]) ;
								Nodes2[hb] = Nodes_deleteind(j, Nodes2[hb]) ;
								
								//Nodes2[ha] = Nodes_addelement(mig2[mig2_n + 3*mig2_N], ttmp, Nodes2[ha] ) ;
								*tNodes2 = Nodes_addelement(mig2[mig2_n + 3*mig2_N], ttmp, *tNodes2) ;
								
								//Nodes[ha] = Nodes_deleteelement(Nodes[ha].minNode,  Nodes[ha] ) ;
								NI[ha]--;
        						
        						//update migration matrix
        						mig = Migmx_addnode(mig_n, mig, *PR_T2[0][i], ha, hb, mig2[mig2_n + 3*mig2_N]) ;
        						//Rprintf("Old migration mx at %8.4f\n", *PR_T2[0][i]) ;
        						for (j=0 ; j<mig2_N ; j++)
								{
								//	Rprintf("%8.4f %8.4f %8.4f %8.4f\n", mig2[j], mig2[j + mig2_N], mig2[j + 2*mig2_N], mig2[j + 3*mig2_N] ) ;
								}
								mig_n++ ;
								
								mig2_n++ ;
        						Ntrans-- ;
        						
        						if (mig2_n == mig2_N) 
        							break ;
        					}
        				}
        			
        			
        			
        		}
        		// separate two different kind of transmissions
        		//Ntrans = (sumOBS > 0 & ha!=hb) ? Ntrans - (int) mig2[mig2_n * 3 + (mig_n - 1)] : Ntrans ;
        		//Rprintf("%i %i\n", Nodes[ha].N, Nodes2[ha].N) ;
         		
         		for (bn=0 ; bn < Ntrans ; bn++)
        		{
        			ov = Nodes_pick_sets2(ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], NI, ha, hb, has_old) ;
   		  		  	//Rprintf("\n%i\t%i\t%i\t%i\n", ov[0], ov[1], ha, hb) ;
   		  		  	
   		  		  	switch(ov[1])
   		  		  	{
   		  		  		case 0:
   		  		  			NI[hb]-- ;
   		  		  			switch(ov[0])
   		  		  			{
   		  		  				case 0:
   		  		  					NI[ha]-- ;
   		  		  					NI_less++;
   		  		  				break ;
   		  		  				case 1:
   		  		  					ch = Nodes_pick_child(ch, Nodes[ha]) ;
   		  		  					ttmp = llist_get_el_d(ch[1], Nodes[ha].t) ;
									*tNodes = Nodes_addelement(ch[0], ttmp, *tNodes) ;
				  					//Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
									Nodes[ha] = Nodes_deleteind(ch[1], Nodes[ha]) ;
								
								break ;
								case 2:
									ch = Nodes_pick_child(ch, Nodes2[ha]) ;
   		  		  					ttmp = llist_get_el_d(ch[1], Nodes2[ha].t) ;
									*tNodes2 = Nodes_addelement(ch[0], ttmp, *tNodes2) ;
				   	  				//Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
				  					Nodes2[ha] = Nodes_deleteind(ch[1], Nodes2[ha]) ;
				  				
				  				break ;
   		  		  			}
   		  		  		break ;
   		  		  		case 1:
   		  		  			switch(ov[0])
   		  		  			{
   		  		  				case 0:
   		  		  					NI[ha]-- ;
   		  		  					ch = Nodes_pick_child(ch, Nodes[hb]) ;
   		  		  					ttmp = llist_get_el_d(ch[1], Nodes[hb].t) ;
   		  		  					//Nodes[hb] = Nodes_deleteelement(ch[0], Nodes[hb]) ;
									Nodes[hb] = Nodes_deleteind(ch[1], Nodes[hb]) ;
									
									*tNodes = Nodes_addelement(ch[0], ttmp, *tNodes) ;
									if (ha!=hb)
									{
										mig = Migmx_addnode(mig_n, mig, *PR_T2[0][i], ha, hb, ch[0]) ;
										mig_n++ ;
   		  		  					}
   		  		  					if (ha != hb) 
								 	{
								 		ce1++ ;
							 		}
						  			Anc[hb]-- ;
						  			Anc[ha]++ ;

   		  		  				break ;
   		  		  				case 1:
   		  		  					ch = Nodes_pick_pair(ch, Nodes[ha], Nodes[hb], ha, hb) ;
									ttmp = llist_get_el_d(ch[2], Nodes[ha].t) - *PR_T2[0][i] ;
									ttmp2 = llist_get_el_d(ch[3], Nodes[hb].t) - *PR_T2[0][i] ;
									//Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
									Nodes[ha] = Nodes_deleteind(ch[2], Nodes[ha]) ;
									Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
									*tNodes = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes) ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[0], ttmp) ;
									ei-- ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[1], ttmp2) ;
									ei-- ;
									tr->nodelabel[intNode-ntips] = ha ;
									//Rprintf("%i\t%i\t%8.4f\n", intNode, ch[0], ttmp) ;
									//Rprintf("%i\t%i\t%8.4f\n", intNode, ch[1], ttmp2) ;
									
									intNode -- ;
									Anc[hb] -- ;
									ce++ ;
								break ;
								case 2:
									
									ch = Nodes_pick_pair(ch, Nodes2[ha], Nodes[hb], ha, -1) ;
									/*Rprintf("\n%i\t%i\t%i\t%i\n", ov[0], ov[1], ch[0], ch[1]) ;
   		  		  					llist_print(Nodes[ha].n) ;
									Rprintf(" ha\n") ;
									llist_print(Nodes[hb].n) ;
									Rprintf(" hb\n") ;
        							llist_print(Nodes2[ha].n) ;
			Rprintf(" ha2\n") ;
			llist_print(Nodes2[hb].n) ;
			Rprintf(" hb2\n") ;*/
        	
        	
									ttmp = llist_get_el_d(ch[2], Nodes2[ha].t) - *PR_T2[0][i] ;
									ttmp2 = llist_get_el_d(ch[3], Nodes[hb].t) - *PR_T2[0][i] ;
									//Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
									Nodes2[ha] = Nodes_deleteind(ch[2], Nodes2[ha]) ;
									
									//Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
									Nodes[hb] = Nodes_deleteind(ch[3], Nodes[hb]) ;
									*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
									if (bt_n >= 0)
									{
										tr_old = phylo_replacenode(tr_old, ((bt_n+1)*2)-1, ch[0], intNode) ;
									}
									mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, ch[0], intNode) ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[0], ttmp) ;
									ei-- ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[1], ttmp2) ;
									ei-- ;
									tr->nodelabel[intNode-ntips] = ha ;
									
									intNode -- ;
									Anc[hb] -- ;
									ce++ ;
							
								break ;
   		  		  			}
   		  		  		break ;
   		  		  		case 2:
   		  		  			switch(ov[0])
   		  		  			{
   		  		  				case 0:
   		  		  					NI[ha]-- ;
   		  		  					ch = Nodes_pick_child(ch, Nodes2[hb]) ;
   		  		  					ttmp = llist_get_el_d(ch[1], Nodes2[hb].t) ;
						  		//	Nodes2[hb] = Nodes_deleteelement(ch[0], Nodes2[hb]) ;
									Nodes2[hb] = Nodes_deleteind(ch[1], Nodes2[hb]) ;
									*tNodes2 = Nodes_addelement(ch[0], ttmp, *tNodes2) ;
								break ;
   		  		  				case 1:
   		  		  					ch = Nodes_pick_pair(ch, Nodes[ha], Nodes2[hb], ha, -1) ;
   		  		  					ttmp = llist_get_el_d(ch[2], Nodes[ha].t) - *PR_T2[0][i] ;
									ttmp2 = llist_get_el_d(ch[3], Nodes2[hb].t) - *PR_T2[0][i] ;
									//Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
									
									//Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
									Nodes[ha] = Nodes_deleteind(ch[2], Nodes[ha]) ;
									
									Nodes2[hb] = Nodes_deleteind(ch[3], Nodes2[hb]) ;
									
									*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
									if (bt_n >= 0)
									{
										tr_old = phylo_replacenode(tr_old, ((bt_n+1)*2)-1, ch[1], intNode) ;
									}
									mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, ch[1], intNode) ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[0], ttmp) ;
									ei-- ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[1], ttmp2) ;
									ei-- ;
									tr->nodelabel[intNode-ntips] = ha ;
							
									intNode -- ;
									Anc[hb] -- ;
									ce++ ;
							
   		  		  				break ;
								case 2:
									Rprintf("not allowed\n") ;
				  				break ;
   		  		  			}
   		  		  		break ;
   		  		  	}
					
            	  	
        	}	
        	/*llist_print(Nodes[ha].n) ;
			Rprintf(" ha\n") ;
			llist_print(Nodes[hb].n) ;
			Rprintf(" hb\n") ;
        	llist_print(Nodes2[ha].n) ;
			Rprintf(" ha2\n") ;
			llist_print(Nodes2[hb].n) ;
			Rprintf(" hb2\n") ;
        	*/
        		// move in new nodes	
        		//Rprintf("Moving %i nodes\n", tNodes->N) ;
        		//Rprintf("%i nodes in ha, %i nodes in tNodes\n", Nodes[ha].N, tNodes->N) ;
        		
        		NI[ha] = NI[ha] + NI_less ;
        		Nodes[ha] =  Nodes_copy_all(Nodes[ha], *tNodes)	;			
        			
        		*tNodes =  Nodes_delete_all(*tNodes) ;
        		Nodes2[ha] =  Nodes_copy_all(Nodes2[ha], *tNodes2)	;				
        		*tNodes2 =  Nodes_delete_all(*tNodes2) ;
        		
        		
				
			}
			}
			else
			{
				//Rprintf("end \n") ;
			}
			
			
		//	Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[ha], Anc[hb]) ;
			
		
		}
		//Rprintf("ll=%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[0][i], *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, bt_n, Anc[0], Anc[1], Anc[2], Anc[3], Anc[4], Nodes[0].N, Nodes2[0].N, OBS[0], OBS2[0], Ancsum, intNode) ;
	//Rprintf("%i\t%i\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ha+1, hb+1, *PR_T2[0][i], *PR_T2[3][i], ce, ce1, bt_n, Nodes[0].N, Nodes2[0].N, Nodes[1].N, Nodes2[1].N, intNode, NI[0], NI[1]) ;
			
		
}

	for (j= 0; j< 4 ; j++)
	{
		Free(PR_T2[j]) ;
	}
	
	
	for (j=0 ; j<NHosts ; j++)
	{
		llist_destroy_i(Nodes[j].n)  ;
		llist_destroy_d(Nodes[j].t)  ;
		llist_destroy_i(Nodes2[j].n)  ;
		llist_destroy_d(Nodes2[j].t)  ;
		free(lambdao[j]) ;
	}
	llist_destroy_i(tNodes->n)  ;
	llist_destroy_d(tNodes->t)  ;
	llist_destroy_i(tNodes2->n)  ;
	llist_destroy_d(tNodes2->t)  ;
	Free(Nodes) ; 
	Free(Nodes2) ;
	Free(tNodes) ;
	Free(tNodes2) ;
	Free(Tend) ;  
	Free(ch) ;
	Free(ov) ;
	Free(PR_T2) ;
 	free(SNsum) ;
 	free(Anc) ;
 	free(lambdao) ;
 	if (has_old == 1)
 	{
 		free(mig2) ;
 		
 		for (i=0 ; i<tr_old->NNode ; i++)
 		{
 			Free(bt_old[i]) ;
 		}
 		
 		Free(bt_old) ;
 		phylo_free(tr_old) ;
 		
	}
	Free(NI) ;

	//Rprintf("end\n") ;
	sp->tr = tr ;
	
	mig_n-- ;
	if (mig_n > 0)
	{
		//if (sp->mig == NULL)
		//{
			sp->mig = (double *)calloc(4 * mig_n, sizeof(double)) ;	
			//Rprintf("%p\n", sp->mig) ;
		//}
		//else
		//{
		//	sp->mig = realloc(sp->mig, 4 * mig_n * sizeof(double)) ;
		//}
	i = 0 ;
	for (j=0 ; j<4 ; j++)
	{
		for (k=0 ; k<mig_n ; k++)
		{
			sp->mig[i++] = mig[k][j] ;
			
		}
		
	}
		//Rprintf("%8.4f\n", sp->mig) ;
		
	for (k=0 ; k<mig_n ; k++)
	{
		Free(mig[k]) ;
	}
	Free(mig) ;
	}
	sp->mig_n = mig_n ;
	//Rprintf("end\n") ;
	
}


SEXP tree_reconstruct_with_partialll3(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS2, SEXP R_OBS, SEXP R_tr_old, SEXP R_mig2)
{
   int  st=0, st2=-1, *minNodes, mn, r1, r2, *ch, ei, intNode;
  	int bn , m,  hi, hj, Ii, Sj, *Anc, kt, ktt, Ancsum=0, ce = 0, ce1= 0;
  	double *PR_bt, **bt, v1, v2, *vvec, ll=0.0, lambdao_sum=0.0, *PR_ll, ll_i, small_val=1e-10 ;

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
	int OBS3[NHosts] ;
	int sumOBS=0, i1 ;
	int ll_calc = 0 ;
	phylo *tr_old ;
	for (i1=0 ; i1 <NHosts ; i1++)
	{
		sumOBS += OBS2[i1] ;
		OBS3[i1] = OBS2[i1] + OBS[i1] ;
	}
	//Rprintf("%i\n", sumOBS) ;
	double **bt_old ;
	double *mig2 ;
	double ttmp, ttmp2 ;
	int mig2_N, mig2_n=0 ;
	int bt_n, i ;
	int has_old = (sumOBS>0) ? 1 : 0 ;
	if (has_old==1)
	{
		tr_old = R_to_phylo(R_tr_old) ;
		// calculate branching times
		bt_old = phylo_bt(tr_old, ST, OBS2, NHosts) ;
		bt_n = tr_old->NNode - 1 ;
		mig2_N = INTEGER(getAttrib(R_mig2, R_DimSymbol))[0] ;
		mig2 = calloc(mig2_N * 4, sizeof(double))  ;
		for (i=0 ; i< (mig2_N * 4); i++)
		{
			mig2[i] = REAL(coerceVector(R_mig2, REALSXP))[i] ; 
		}
	}
	
	
	/*************** INTERNAL VARIABLES FOR ALGORITHM ***********************************/
	
	int j, k=0, ha, hb, oa, ob, Ntrans ;
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
		ntips+= OBS3[i] ;
		lambdao[i] = calloc(NHosts, sizeof(double)) ;
	}
	hnode2 * Nodes = Calloc(NHosts,hnode2) ; 
	hnode2 * Nodes2 = Calloc(NHosts, hnode2) ;
	hnode2 * tNodes = Calloc(1,hnode2) ; 
	hnode2 * tNodes2 = Calloc(1, hnode2) ;
	int *NI = Calloc(NHosts, int) ;
	
	
	
	hnode2 hnode2v ;
	
	for (i=0 ; i<NHosts ; i++)
	{
		Nodes[i] = Nodes_initialise(Nodes[i]) ;
		Nodes2[i] = Nodes_initialise(Nodes2[i]) ;
	}
	
	// initialise temporary nodes
	*tNodes = Nodes_initialise(*tNodes) ;
	*tNodes2 = Nodes_initialise(*tNodes2) ;
	
	phylo *tr = phylo_create(ntips) ;
	
	ch = Calloc(4, int) ;
	int * ov = Calloc(2, int) ;
	
	ei = tr->Nedge - 1 ;
	intNode = 2*tr->NNode  ;
	Anc = calloc(NHosts, sizeof(int)) ;
	vvec = calloc(NS, sizeof(double)) ;
	
	// matrix to store the number of migrating lineages
	double ** mig ;
	int mig_n = 1 ;
	int NI_less ;
	
	
  	j = 0;
  	// maybe still some infected left at end:
  
	for (i=0 ; i<NHosts ; i++)
	{
		//SNsum[i] = SN[i] + j ;
    	//j+= SN[i] ;
    	v1 = OBS[i] + OBS2[i] ;
    	SNsum[i] = v1 + j ;
    	j+= v1 ;
    	NI[i] = PR_I[TN * (i+1) - 1] ;
		Anc[i] = 0 ;
	}


/*******************************ALGORITHM************************************************/

	for (i=TN-1 ; i>-1 ; i--)
	{
		//Rprintf("%i\n", i) ;
		ll_i = 0 ;
		ce = 0 ;
		ce1 = 0 ;
    // for each event get hosts
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		
   								//Rprintf("\nNodes in 0 at %i=\n", i) ;
        						//llist_print(Nodes[0].n) ;
        						
    // if terminal event
   		if (i==Tend[ha])
    	{
    		//Rprintf("Terminal event %i\n", i) ;
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
		  	
		  	NI[ha] += -*PR_T2[3][i] - (OBS[ha] + OBS2[ha]) ;
		  	//if ((OBS[ha] + OBS2[ha])< -*PR_T2[3][i])
		  	//{
			//	minNodes[ha] = st2+1 ;  
		  	//}
		  	Anc[ha]+=OBS[ha] + OBS2[ha];  
		  	/* 
		  	llist_print(Nodes[ha].n) ;
	 		Rprintf(" ha\n") ;	
	 		llist_print(Nodes2[ha].n) ;
			Rprintf(" ha2\n") ;
			*/
		}
    	else
		{
			/*
			Ancsum = 0 ;
			for (kt = 0 ; kt < NHosts ; kt++)
			{
				Ancsum += Anc[kt] ; 
			}
			*/
			//if (i > 0 & Ancsum > 1)
			
			if ((i > 0) & (intNode >= ntips))
			//if ((i > 0))
			{
				// calculate rates for each host
				kt = *PR_T2[3][i] ;	
				NI_less = 0 ;
					
				if (ll_calc == 1)
				{
				lambdao =  calculateSIR_rates(lambdao, B, PR_I, PR_S, ha, hb, i, kt, NHosts, NS, TN, Anc, OBS2) ;
				lambdao_sum = calculatelambdaosum(lambdao, NHosts) ;
				ll_i = -lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
				}
				if ((int) *PR_T2[3][i] == -1 )
				{
					// death but birth in reverse.
					// can calculate the probability of no coalescence has occurred
				
				// minnode
					
        			NI[ha]++ ;
					
					//Nodes[ha] = Nodes_addelement(Nodes[ha].minNode - 1, *PR_T2[0][i], Nodes[ha]) ;
					
				} 
				
				//if ((int) *PR_T2[3][i] == 1 )
				//{
					// 1 coalescent event within host
				//	ov = Nodes_pick_sets(ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], ha, hb) ;

				//}
				//Rprintf("%i\n", bt_n) ;
				if ((int) *PR_T2[3][i] > 0 )
				{
					Ntrans = (int)*PR_T2[3][i] ;
					// if coincides with coalescent interval among old observed tree
					if (has_old==1) //oa,ob=11
					{			
						if (bt_n >= 0)
						{
						while ((bt_old[bt_n][0]-small_val) > *PR_T2[0][i-1])
						{
							ch[0] = tr_old->edge[bt_n*2 + 1][1] ;
        					ch[1] = tr_old->edge[bt_n*2][1] ;
        					/*
							Rprintf("Existing coalescent event between %i from %i and  %i from %i at %8.4f, last bt = %8.4f, bt=%8.4f\n", ch[0], ha, ch[1], hb, *PR_T2[0][i], *PR_T2[0][i-1], bt_old[bt_n][0] ) ;
							llist_print(Nodes2[0].n) ;
							Rprintf(" ha2\n") ;
							llist_print(Nodes2[1].n) ;
							Rprintf(" hb2\n") ;
							Rprintf("Nodes in 0=\n") ;
        				    llist_print(Nodes[0].n) ;
        				    Rprintf("\nNodes in 1=\n") ;
        				    llist_print(Nodes[1].n) ;
							*/
							// create coalescent event in new tree
        					
        					ch[2] = llist_get_ind_i(ch[0], Nodes2[ha].n) ;
        					ch[3] = llist_get_ind_i(ch[1], Nodes2[hb].n) ;
        					//Rprintf("birth %i %i %i %i %i\n", ch[0], ch[1], ch[2], ch[3], ei) ;
        					// update new tree
        					tr = phylo_addnewedge(tr, ei, intNode, ch[0], llist_get_el_d(ch[2], Nodes2[ha].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr = phylo_addnewedge(tr, ei, intNode, ch[1], llist_get_el_d(ch[3], Nodes2[hb].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr->nodelabel[intNode-ntips] = ha ;
							//Rprintf("added edges\n") ;
							
							
							//Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
							//Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
							
							Nodes2[ha] = Nodes_deleteind(ch[2], Nodes2[ha]) ;
							Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
							
							
							// add later
							
							*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
							//Nodes2[ha] = Nodes_addelement(intNode, *PR_T2[0][i], Nodes2[ha]) ;
							
							if (bt_n > 0)
								tr_old = phylo_replacenode(tr_old, (bt_n*2)-1, tr_old->edge[bt_n*2][0], intNode) ;
							mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, tr_old->edge[bt_n*2][0], intNode) ;
							//Rprintf("Replacing %i with %i, mig2_n = %i, bt_n =%i\n", tr_old->edge[bt_n*2][0], intNode, mig2_n, bt_n ) ;
							
							intNode -- ;
							Anc[hb] -- ;
							//OBS2[hb]-- ;
							//ce++ ;
        					//Rprintf("%i\n", ei) ;
        					bt_n -- ;	
							Ntrans -- ;
						
							
						/*
							Rprintf("old tree\n") ;
							if (sumOBS > 0)
							{
								for (j=0 ; j<tr_old->Nedge ; j++)
								{
									//Rprintf("%i %i\n", tr_old->edge[j][0], tr_old->edge[j][1]) ;
								}
								
							}
							Rprintf("new tree\n") ;
							for (j=0 ; j<tr->Nedge ; j++)
							{
								Rprintf("%i %i\n", tr->edge[j][0], tr->edge[j][1]) ;
							}
							
							*/
							if (bt_n < 0) break ;
							
        				}
        				// look at migrations
        				}
        				if (mig2_n < mig2_N)
        				{
        					while ((mig2[mig2_n]-small_val) > *PR_T2[0][i-1])
        					{
        						//Rprintf("Old migration, node=%i, ha=%i, hb=%i\n", (int) mig2[mig2_n + 3*mig2_N], ha, hb) ;
        						
        						j = llist_get_ind_i(mig2[mig2_n + 3*mig2_N], Nodes2[hb].n) ;
        						ttmp = llist_get_el_d(j, Nodes2[hb].t) ;
						  		
								//Nodes2[hb] = Nodes_deleteelement(mig2[mig2_n + 3*mig2_N], Nodes2[hb]) ;
								Nodes2[hb] = Nodes_deleteind(j, Nodes2[hb]) ;
								
								//Nodes2[ha] = Nodes_addelement(mig2[mig2_n + 3*mig2_N], ttmp, Nodes2[ha] ) ;
								*tNodes2 = Nodes_addelement(mig2[mig2_n + 3*mig2_N], ttmp, *tNodes2) ;
								
								//Nodes[ha] = Nodes_deleteelement(Nodes[ha].minNode,  Nodes[ha] ) ;
								NI[ha]--;
        						
        						//update migration matrix
        						mig = Migmx_addnode(mig_n, mig, *PR_T2[0][i], ha, hb, mig2[mig2_n + 3*mig2_N]) ;
        						//Rprintf("Old migration mx at %8.4f\n", *PR_T2[0][i]) ;
        						for (j=0 ; j<mig2_N ; j++)
								{
								//	Rprintf("%8.4f %8.4f %8.4f %8.4f\n", mig2[j], mig2[j + mig2_N], mig2[j + 2*mig2_N], mig2[j + 3*mig2_N] ) ;
								}
								mig_n++ ;
								
								mig2_n++ ;
        						Ntrans-- ;
        						
        						if (mig2_n == mig2_N) 
        							break ;
        					}
        				}
        			
        			
        			
        		}
        		// separate two different kind of transmissions
        		//Ntrans = (sumOBS > 0 & ha!=hb) ? Ntrans - (int) mig2[mig2_n * 3 + (mig_n - 1)] : Ntrans ;
        		//Rprintf("%i %i\n", Nodes[ha].N, Nodes2[ha].N) ;
         		
         		for (bn=0 ; bn < Ntrans ; bn++)
        		{
        			ov = Nodes_pick_sets2(ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], NI, ha, hb, has_old) ;
   		  		  	//Rprintf("\n%i\t%i\t%i\t%i\n", ov[0], ov[1], ha, hb) ;
   		  		  	
   		  		  	switch(ov[1])
   		  		  	{
   		  		  		case 0:
   		  		  			NI[hb]-- ;
   		  		  			switch(ov[0])
   		  		  			{
   		  		  				case 0:
   		  		  					NI[ha]-- ;
   		  		  					NI_less++;
   		  		  				break ;
   		  		  				case 1:
   		  		  					ch = Nodes_pick_child(ch, Nodes[ha]) ;
   		  		  					ttmp = llist_get_el_d(ch[1], Nodes[ha].t) ;
									*tNodes = Nodes_addelement(ch[0], ttmp, *tNodes) ;
				  					//Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
									Nodes[ha] = Nodes_deleteind(ch[1], Nodes[ha]) ;
								
								break ;
								case 2:
									ch = Nodes_pick_child(ch, Nodes2[ha]) ;
   		  		  					ttmp = llist_get_el_d(ch[1], Nodes2[ha].t) ;
									*tNodes2 = Nodes_addelement(ch[0], ttmp, *tNodes2) ;
				   	  				//Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
				  					Nodes2[ha] = Nodes_deleteind(ch[1], Nodes2[ha]) ;
				  				
				  				break ;
   		  		  			}
   		  		  		break ;
   		  		  		case 1:
   		  		  			switch(ov[0])
   		  		  			{
   		  		  				case 0:
   		  		  					NI[ha]-- ;
   		  		  					ch = Nodes_pick_child(ch, Nodes[hb]) ;
   		  		  					ttmp = llist_get_el_d(ch[1], Nodes[hb].t) ;
   		  		  					//Nodes[hb] = Nodes_deleteelement(ch[0], Nodes[hb]) ;
									Nodes[hb] = Nodes_deleteind(ch[1], Nodes[hb]) ;
									
									*tNodes = Nodes_addelement(ch[0], ttmp, *tNodes) ;
									if (ha!=hb)
									{
										mig = Migmx_addnode(mig_n, mig, *PR_T2[0][i], ha, hb, ch[0]) ;
										mig_n++ ;
   		  		  					}
   		  		  					if (ha != hb) 
								 	{
								 		ce1++ ;
							 		}
						  			Anc[hb]-- ;
						  			Anc[ha]++ ;

   		  		  				break ;
   		  		  				case 1:
   		  		  					ch = Nodes_pick_pair(ch, Nodes[ha], Nodes[hb], ha, hb) ;
									ttmp = llist_get_el_d(ch[2], Nodes[ha].t) - *PR_T2[0][i] ;
									ttmp2 = llist_get_el_d(ch[3], Nodes[hb].t) - *PR_T2[0][i] ;
									//Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
									Nodes[ha] = Nodes_deleteind(ch[2], Nodes[ha]) ;
									Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
									*tNodes = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes) ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[0], ttmp) ;
									ei-- ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[1], ttmp2) ;
									ei-- ;
									tr->nodelabel[intNode-ntips] = ha ;
							
									intNode -- ;
									Anc[hb] -- ;
									ce++ ;
								break ;
								case 2:
									
									ch = Nodes_pick_pair(ch, Nodes2[ha], Nodes[hb], ha, -1) ;
									/*Rprintf("\n%i\t%i\t%i\t%i\n", ov[0], ov[1], ch[0], ch[1]) ;
   		  		  					llist_print(Nodes[ha].n) ;
									Rprintf(" ha\n") ;
									llist_print(Nodes[hb].n) ;
									Rprintf(" hb\n") ;
        							llist_print(Nodes2[ha].n) ;
			Rprintf(" ha2\n") ;
			llist_print(Nodes2[hb].n) ;
			Rprintf(" hb2\n") ;*/
        	
        	
									ttmp = llist_get_el_d(ch[2], Nodes2[ha].t) - *PR_T2[0][i] ;
									ttmp2 = llist_get_el_d(ch[3], Nodes[hb].t) - *PR_T2[0][i] ;
									//Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
									Nodes2[ha] = Nodes_deleteind(ch[2], Nodes2[ha]) ;
									
									//Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
									Nodes[hb] = Nodes_deleteind(ch[3], Nodes[hb]) ;
									*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
									if (bt_n >= 0)
									{
										tr_old = phylo_replacenode(tr_old, ((bt_n+1)*2)-1, ch[0], intNode) ;
									}
									mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, ch[0], intNode) ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[0], ttmp) ;
									ei-- ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[1], ttmp2) ;
									ei-- ;
									tr->nodelabel[intNode-ntips] = ha ;
									
									intNode -- ;
									Anc[hb] -- ;
									ce++ ;
							
								break ;
   		  		  			}
   		  		  		break ;
   		  		  		case 2:
   		  		  			switch(ov[0])
   		  		  			{
   		  		  				case 0:
   		  		  					NI[ha]-- ;
   		  		  					ch = Nodes_pick_child(ch, Nodes2[hb]) ;
   		  		  					ttmp = llist_get_el_d(ch[1], Nodes2[hb].t) ;
						  		//	Nodes2[hb] = Nodes_deleteelement(ch[0], Nodes2[hb]) ;
									Nodes2[hb] = Nodes_deleteind(ch[1], Nodes2[hb]) ;
									*tNodes2 = Nodes_addelement(ch[0], ttmp, *tNodes2) ;
								break ;
   		  		  				case 1:
   		  		  					ch = Nodes_pick_pair(ch, Nodes[ha], Nodes2[hb], ha, -1) ;
   		  		  					ttmp = llist_get_el_d(ch[2], Nodes[ha].t) - *PR_T2[0][i] ;
									ttmp2 = llist_get_el_d(ch[3], Nodes2[hb].t) - *PR_T2[0][i] ;
									//Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
									
									//Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
									Nodes[ha] = Nodes_deleteind(ch[2], Nodes[ha]) ;
									
									Nodes2[hb] = Nodes_deleteind(ch[3], Nodes2[hb]) ;
									
									*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
									if (bt_n >= 0)
									{
										tr_old = phylo_replacenode(tr_old, ((bt_n+1)*2)-1, ch[1], intNode) ;
									}
									mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, ch[1], intNode) ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[0], ttmp) ;
									ei-- ;
									tr = phylo_addnewedge(tr, ei, intNode, ch[1], ttmp2) ;
									ei-- ;
									tr->nodelabel[intNode-ntips] = ha ;
							
									intNode -- ;
									Anc[hb] -- ;
									ce++ ;
							
   		  		  				break ;
								case 2:
									Rprintf("not allowed\n") ;
				  				break ;
   		  		  			}
   		  		  		break ;
   		  		  	}
					
            	  	
        	}	
        	/*llist_print(Nodes[ha].n) ;
			Rprintf(" ha\n") ;
			llist_print(Nodes[hb].n) ;
			Rprintf(" hb\n") ;
        	llist_print(Nodes2[ha].n) ;
			Rprintf(" ha2\n") ;
			llist_print(Nodes2[hb].n) ;
			Rprintf(" hb2\n") ;
        	*/
        		// move in new nodes	
        		//Rprintf("Moving %i nodes\n", tNodes->N) ;
        		//Rprintf("%i nodes in ha, %i nodes in tNodes\n", Nodes[ha].N, tNodes->N) ;
        		
        		NI[ha] = NI[ha] + NI_less ;
        		Nodes[ha] =  Nodes_copy_all(Nodes[ha], *tNodes)	;			
        			
        		*tNodes =  Nodes_delete_all(*tNodes) ;
        		Nodes2[ha] =  Nodes_copy_all(Nodes2[ha], *tNodes2)	;				
        		*tNodes2 =  Nodes_delete_all(*tNodes2) ;
        		
        		
				
			}
			if (ll_calc == 1)
			{
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
			}
			
			ll += ll_i ;
			}
			else
			{
				//Rprintf("end \n") ;
			}
			
			
		//	Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[ha], Anc[hb]) ;
			
		
		}
		//Rprintf("ll=%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[0][i], *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, bt_n, Anc[0], Anc[1], Anc[2], Anc[3], Anc[4], Nodes[0].N, Nodes2[0].N, OBS[0], OBS2[0], Ancsum, intNode) ;
	//Rprintf("%i\t%i\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ha+1, hb+1, *PR_T2[0][i], *PR_T2[3][i], ce, ce1, bt_n, Nodes[0].N, Nodes2[0].N, Nodes[1].N, Nodes2[1].N, intNode, NI[0], NI[1]) ;
			
		
}


/*******************************RESULT TO OUTPUT*****************************************/

	PROTECT(R_List = allocVector(VECSXP ,4)) ;
	PROTECT(R_bt=allocMatrix(REALSXP,tr->NNode,4)) ;
	PROTECT(R_ll = allocVector(REALSXP, 1)) ;
	PROTECT(R_mig = allocMatrix(REALSXP, mig_n-1, 4)) ;
	PR_bt = REAL(R_bt) ;
	
	/*
	bt = phylo_bt(tr, ST, OBS3, NHosts) ;
	for (j=0 ; j<tr->NNode ; j++)
	{
		for (i=0;i<4;i++)
		{
			PR_bt[j+i*tr->NNode]=(i>0) ? (bt[j][i] + 1) : bt[j][i] ;
			
		}
	}*/
	
	for (j=0 ; j<(mig_n - 1) ; j++)
	{
		for (i=0;i<4;i++)
		{
			REAL(R_mig)[j+i*(mig_n - 1)]=  mig[j][i] ;
			
		}
	}
	
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
		//Free(bt[j]) ;
		
	}
	//Free(bt) ;
	
	for (j=0 ; j<NHosts ; j++)
	{
		llist_destroy_i(Nodes[j].n)  ;
		llist_destroy_d(Nodes[j].t)  ;
		llist_destroy_i(Nodes2[j].n)  ;
		llist_destroy_d(Nodes2[j].t)  ;
		free(lambdao[j]) ;
	}
	llist_destroy_i(tNodes->n)  ;
	llist_destroy_d(tNodes->t)  ;
	llist_destroy_i(tNodes2->n)  ;
	llist_destroy_d(tNodes2->t)  ;
	free(Nodes) ; 
	free(Tend) ;  
	free(ch) ;
	free(PR_T2) ;
 	free(SNsum) ;
 	free(Anc) ;
 	free(lambdao) ;
 	free(vvec) ;
 	if (has_old == 1)
 	{
 		free(mig2) ;
	}
	Free(NI) ;
	return R_List;
}


SEXP tree_reconstruct_with_partialll2(SEXP R_sir, SEXP R_NHosts, SEXP R_dat, SEXP R_B, SEXP R_NS, SEXP R_BN, SEXP R_OBS2, SEXP R_OBS, SEXP R_tr_old, SEXP R_mig2)
{
   int  st=0, st2=-1, *minNodes, mn, r1, r2, *ch, ei, intNode;
  	int bn , m,  hi, hj, Ii, Sj, *Anc, kt, ktt, Ancsum=0, ce = 0, ce1= 0;
  	double *PR_bt, **bt, v1, v2, *vvec, ll=0.0, lambdao_sum=0.0, *PR_ll, ll_i, small_val=1e-10 ;

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
	int OBS3[NHosts] ;
	int sumOBS=0, i1 ;
	int ll_calc = 0 ;
	phylo *tr_old ;
	for (i1=0 ; i1 <NHosts ; i1++)
	{
		sumOBS += OBS2[i1] ;
		OBS3[i1] = OBS2[i1] + OBS[i1] ;
	}
	//Rprintf("%i\n", sumOBS) ;
	double **bt_old ;
	double *mig2 ;
	double ttmp, ttmp2 ;
	int mig2_N, mig2_n=0 ;
	int bt_n, i ;
	if (sumOBS > 0)
	{
		tr_old = R_to_phylo(R_tr_old) ;
		// calculate branching times
		bt_old = phylo_bt(tr_old, ST, OBS2, NHosts) ;
		bt_n = tr_old->NNode - 1 ;
		mig2_N = INTEGER(getAttrib(R_mig2, R_DimSymbol))[0] ;
		mig2 = calloc(mig2_N * 4, sizeof(double))  ;
		for (i=0 ; i< (mig2_N * 4); i++)
		{
			mig2[i] = REAL(coerceVector(R_mig2, REALSXP))[i] ; 
		}
	}
	
	
	/*************** INTERNAL VARIABLES FOR ALGORITHM ***********************************/
	
	int j, k=0, ha, hb, oa, ob, Ntrans ;
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
	hnode2 * tNodes = Calloc(1,hnode2) ; 
	hnode2 * tNodes2 = Calloc(1, hnode2) ;
	
	hnode2 hnode2v ;
	
	for (i=0 ; i<NHosts ; i++)
	{
		Nodes[i] = Nodes_initialise(Nodes[i]) ;
		Nodes2[i] = Nodes_initialise(Nodes2[i]) ;
	}
	
	// initialise temporary nodes
	*tNodes = Nodes_initialise(*tNodes) ;
	*tNodes2 = Nodes_initialise(*tNodes2) ;
	
	phylo *tr = phylo_create(ntips) ;
	
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
    	v1 = OBS[i] + OBS2[i] ;
    	SNsum[i] = v1 + j ;
    	j+= v1 ;
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
		//Rprintf("%i\n", i) ;
		ll_i = 0 ;
		ce = 0 ;
		ce1 = 0 ;
    // for each event get hosts
		ha = (int) *PR_T2[1][i] - 1;
		hb = (int) *PR_T2[2][i] - 1;
		/*
   		Rprintf("\nNodes in 0 at %i=\n", i) ;
        						llist_print(Nodes2[0].n) ;
        						Rprintf("\nNodes in 1 at %i=\n", i) ;
        						llist_print(Nodes2[1].n) ;
        */						
    // if terminal event
   		if (i==Tend[ha])
    	{
    		//Rprintf("Terminal event %i\n", i) ;
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
		  	//llist_print(Nodes[ha].n) ;
	 		//Rprintf(" ha\n") ;	
	 		//llist_print(Nodes2[ha].n) ;
			//Rprintf(" ha2\n") ;
		}
    	else
		{
			/*
			Ancsum = 0 ;
			for (kt = 0 ; kt < NHosts ; kt++)
			{
				Ancsum += Anc[kt] ; 
			}
			*/
			//if (i > 0 & Ancsum > 1)
			
			if ((i > 0) & (intNode >= ntips))
			//if ((i > 0))
			{
				// calculate rates for each host
				kt = *PR_T2[3][i] ;	
				
					
				if (ll_calc == 1)
				{
				lambdao =  calculateSIR_rates(lambdao, B, PR_I, PR_S, ha, hb, i, kt, NHosts, NS, TN, Anc, OBS2) ;
				lambdao_sum = calculatelambdaosum(lambdao, NHosts) ;
				ll_i = -lambdao_sum * (*PR_T2[0][i] - *PR_T2[0][i-1]); 
				}
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
					
					//Nodes[ha] = Nodes_addelement(Nodes[ha].minNode - 1, *PR_T2[0][i], Nodes[ha]) ;
					
				} 
				
				//if ((int) *PR_T2[3][i] == 1 )
				//{
					// 1 coalescent event within host
				//	ov = Nodes_pick_sets(ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], ha, hb) ;

				//}
				//Rprintf("%i\n", bt_n) ;
				if ((int) *PR_T2[3][i] > 0 )
				{
					Ntrans = (int)*PR_T2[3][i] ;
					// if coincides with coalescent interval among old observed tree
					if (sumOBS > 0) //oa,ob=11
					{			
						if (bt_n >= 0)
						{
						while ((bt_old[bt_n][0]-small_val) > *PR_T2[0][i-1])
						{
							ch[0] = tr_old->edge[bt_n*2 + 1][1] ;
        					ch[1] = tr_old->edge[bt_n*2][1] ;
        					/*
							Rprintf("Existing coalescent event between %i from %i and  %i from %i at %8.4f, last bt = %8.4f, bt=%8.4f\n", ch[0], ha, ch[1], hb, *PR_T2[0][i], *PR_T2[0][i-1], bt_old[bt_n][0] ) ;
							llist_print(Nodes2[0].n) ;
							Rprintf(" ha2\n") ;
							llist_print(Nodes2[1].n) ;
							Rprintf(" hb2\n") ;
							Rprintf("Nodes in 0=\n") ;
        				    llist_print(Nodes[0].n) ;
        				    Rprintf("\nNodes in 1=\n") ;
        				    llist_print(Nodes[1].n) ;
							*/
							// create coalescent event in new tree
        					
        					ch[2] = llist_get_ind_i(ch[0], Nodes2[ha].n) ;
        					ch[3] = llist_get_ind_i(ch[1], Nodes2[hb].n) ;
        					//Rprintf("birth %i %i %i %i %i\n", ch[0], ch[1], ch[2], ch[3], ei) ;
        					// update new tree
        					tr = phylo_addnewedge(tr, ei, intNode, ch[0], llist_get_el_d(ch[2], Nodes2[ha].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr = phylo_addnewedge(tr, ei, intNode, ch[1], llist_get_el_d(ch[3], Nodes2[hb].t) - *PR_T2[0][i]) ;
							ei-- ;
							tr->nodelabel[intNode-ntips] = ha ;
							//Rprintf("added edges\n") ;
							Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
							Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
							// add later
							
							*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
							//Nodes2[ha] = Nodes_addelement(intNode, *PR_T2[0][i], Nodes2[ha]) ;
							
							if (bt_n > 0)
								tr_old = phylo_replacenode(tr_old, (bt_n*2)-1, tr_old->edge[bt_n*2][0], intNode) ;
							mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, tr_old->edge[bt_n*2][0], intNode) ;
							//Rprintf("Replacing %i with %i, mig2_n = %i, bt_n =%i\n", tr_old->edge[bt_n*2][0], intNode, mig2_n, bt_n ) ;
							
							intNode -- ;
							Anc[hb] -- ;
							OBS2[hb]-- ;
							//ce++ ;
        					//Rprintf("%i\n", ei) ;
        					bt_n -- ;	
							Ntrans -- ;
						
							
						/*
							Rprintf("old tree\n") ;
							if (sumOBS > 0)
							{
								for (j=0 ; j<tr_old->Nedge ; j++)
								{
									//Rprintf("%i %i\n", tr_old->edge[j][0], tr_old->edge[j][1]) ;
								}
								
							}
							Rprintf("new tree\n") ;
							for (j=0 ; j<tr->Nedge ; j++)
							{
								Rprintf("%i %i\n", tr->edge[j][0], tr->edge[j][1]) ;
							}
							
							*/
							if (bt_n < 0) break ;
							
        				}
        				// look at migrations
        				}
        				if (mig2_n < mig2_N)
        				{
        					while ((mig2[mig2_n]-small_val) > *PR_T2[0][i-1])
        					{
        						//Rprintf("Old migration, node=%i, ha=%i, hb=%i\n", (int) mig2[mig2_n + 3*mig2_N], ha, hb) ;
        						
        						
        						ttmp = llist_get_el_d(llist_get_ind_i(mig2[mig2_n + 3*mig2_N], Nodes2[hb].n), Nodes2[hb].t) ;
						  		
								Nodes2[hb] = Nodes_deleteelement(mig2[mig2_n + 3*mig2_N], Nodes2[hb]) ;
								//Nodes2[ha] = Nodes_addelement(mig2[mig2_n + 3*mig2_N], ttmp, Nodes2[ha] ) ;
								*tNodes2 = Nodes_addelement(mig2[mig2_n + 3*mig2_N], ttmp, *tNodes2) ;
								
								Nodes[ha] = Nodes_deleteelement(Nodes[ha].minNode,  Nodes[ha] ) ;
								
        						
        						//update migration matrix
        						mig = Migmx_addnode(mig_n, mig, *PR_T2[0][i], ha, hb, mig2[mig2_n + 3*mig2_N]) ;
        						//Rprintf("Old migration mx at %8.4f\n", *PR_T2[0][i]) ;
        						for (j=0 ; j<mig2_N ; j++)
								{
								//	Rprintf("%8.4f %8.4f %8.4f %8.4f\n", mig2[j], mig2[j + mig2_N], mig2[j + 2*mig2_N], mig2[j + 3*mig2_N] ) ;
								}
								mig_n++ ;
								
								mig2_n++ ;
        						Ntrans-- ;
        						
        						if (mig2_n == mig2_N) 
        							break ;
        					}
        				}
        			
        			
        			
        		}
        		// separate two different kind of transmissions
        		//Ntrans = (sumOBS > 0 & ha!=hb) ? Ntrans - (int) mig2[mig2_n * 3 + (mig_n - 1)] : Ntrans ;
        		//Rprintf("%i %i\n", Nodes[ha].N, Nodes2[ha].N) ;
         		
         		for (bn=0 ; bn < Ntrans ; bn++)
        		{
        			ov = Nodes_pick_sets(ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], ha, hb) ;
   		  		  	if ((ha!=hb) & (sumOBS > 0))
        			{
        				ov[1] = 0 ;
   		  		  	}
        			
   		  		  	//Rprintf("%i\t%i\n", ov[0], ov[1]) ;
					
            	  	ch = Nodes_pick_children(ch, ov, Nodes[ha], Nodes2[ha], Nodes[hb], Nodes2[hb], ha, hb) ;
					//Rprintf("%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ov[0], ov[1], ch[0], ch[1], Nodes[ha].N, Nodes[hb].N, tNodes->N) ;
					/*
					llist_print(Nodes[ha].n) ;
							Rprintf(" ha\n") ;
							llist_print(Nodes[hb].n) ;
							Rprintf(" hb\n") ;
							llist_print(tNodes->n) ;
							Rprintf(" t\n") ;
							llist_print(tNodes2->n) ;
							Rprintf(" t2\n") ;
							llist_print(Nodes2[ha].n) ;
							Rprintf(" ha2\n") ;
							llist_print(Nodes2[hb].n) ;
							Rprintf(" hb2\n") ;
							*/
					if (ch[1] < 0) // -
				  	{
				  		// oa,ob 10 00
				  		Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
				  		if (ov[0] == 0)
				  		{  		
				  			ttmp = llist_get_el_d(llist_get_ind_i(ch[0], Nodes[ha].n), Nodes[ha].t) ;
							*tNodes = Nodes_addelement(ch[0], ttmp, *tNodes) ;
				  			
				  			Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
				  		
				   	  	}
				   	  	else
				   	  	{
				   	  		/*
				   	  		llist_print(Nodes[ha].n) ;
							Rprintf(" ha\n") ;
							llist_print(Nodes[hb].n) ;
							Rprintf(" hb\n") ;
							llist_print(tNodes->n) ;
							Rprintf(" t\n") ;
							llist_print(tNodes2->n) ;
							Rprintf(" t2\n") ;
							llist_print(Nodes2[ha].n) ;
							Rprintf(" ha2\n") ;
							llist_print(Nodes2[hb].n) ;
							Rprintf(" hb2\n") ;
							*/
        	
				   	  		ttmp = llist_get_el_d(llist_get_ind_i(ch[0], Nodes2[ha].n), Nodes2[ha].t) ;
							*tNodes2 = Nodes_addelement(ch[0], ttmp, *tNodes2) ;
				   	  		
				  			Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
				  		
				   	  	}
				  		
				  	}
				  	else{
				  		// node from hb is observed
					  	if (ch[0] < 0) // +-
					  	{
					  		//oa,ob 00 01
					  		// node from ha is unobserved
					  		// move node from hb into ha.
						  	Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
							if (ov[1]==0)
							{
								ttmp = llist_get_el_d(llist_get_ind_i(ch[1], Nodes[hb].n), Nodes[hb].t) ;
						  		
								Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
								//Nodes[ha] = Nodes_addelement(ch[1], ttmp, Nodes[ha] ) ;
								*tNodes = Nodes_addelement(ch[1], ttmp, *tNodes) ;
								if (ha!=hb)
								{
								mig = Migmx_addnode(mig_n, mig, *PR_T2[0][i], ha, hb, ch[1]) ;
								//Rprintf("New migration : Moved %i from host %i to %i\n", ch[1], hb, ha) ;
								
								mig_n++ ;
								//Rprintf("mig_n = %i\n", mig_n) ;
								}
							}
							else
							{
								ttmp = llist_get_el_d(llist_get_ind_i(ch[1], Nodes2[hb].n), Nodes2[hb].t) ;
						  		//Rprintf("moving %i from %i to %i\n", ch[1], hb, ha) ;
								Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
								*tNodes2 = Nodes_addelement(ch[1], ttmp, *tNodes2) ;
								//Nodes2[ha] = Nodes_addelement(ch[1], ttmp, Nodes2[ha] ) ;
								//Rprintf("Nodes in 0=\n") ;
        						//llist_print(Nodes2[0].n) ;
        						//Rprintf("\nNodes in 1=\n") ;
        						//llist_print(Nodes2[1].n) ;
							}
							
						 	
						 	
						 	if (ha != hb) 
						 	{
						 		ce1++ ;
						 		// record migrations between different hosts.
						 		
						 	}
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
							tr->nodelabel[intNode-ntips] = ha ;
							if (ov[0] == 0)
							{
								
								ttmp = llist_get_el_d(ch[2], Nodes[ha].t) - *PR_T2[0][i] ;
								if (ov[1] == 0)
								{
									//00
									ttmp2 = llist_get_el_d(ch[3], Nodes[hb].t) - *PR_T2[0][i] ;
									Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
									Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
									*tNodes = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes) ;
								
									//Nodes[ha] = Nodes_addelement(intNode, *PR_T2[0][i], Nodes[ha]) ;
									//OBS[hb]-- ;
									
								}
								else
								{
									//01
									ttmp2 = llist_get_el_d(ch[3], Nodes2[hb].t) - *PR_T2[0][i] ;
									Nodes[ha] = Nodes_deleteelement(ch[0], Nodes[ha]) ;
									Nodes2[hb] = Nodes_deleteelement(ch[1], Nodes2[hb]) ;
									*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
									
									//Nodes2[ha] = Nodes_addelement(intNode, *PR_T2[0][i], Nodes2[ha]) ;
									//Rprintf("New coalescent event between %i and %i\n", ch[0], ch[1]) ;
									//Rprintf("bt_n=%i", bt_n) ;
									if (bt_n >= 0)
									{
										tr_old = phylo_replacenode(tr_old, ((bt_n+1)*2)-1, ch[1], intNode) ;
									//Rprintf("replacing %i with %i\n", ch[1], intNode) ;	
									}
									//if (ha != hb)
									//{
									mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, ch[1], intNode) ;
								for (j=0 ; j<mig2_N ; j++)
								{
									//Rprintf("%8.4f %8.4f %8.4f %8.4f %i\n", mig2[j], mig2[j + mig2_N], mig2[j + 2*mig2_N], mig2[j + 3*mig2_N], mig2_n ) ;
								}
									//}
									//OBS2[hb]-- ;
									
								}
								
								
							}
							else
							{
							//10
								ttmp = llist_get_el_d(ch[2], Nodes2[ha].t) - *PR_T2[0][i] ;
								ttmp2 = llist_get_el_d(ch[3], Nodes[hb].t) - *PR_T2[0][i] ;
								Nodes2[ha] = Nodes_deleteelement(ch[0], Nodes2[ha]) ;
								Nodes[hb] = Nodes_deleteelement(ch[1], Nodes[hb]) ;
								//Nodes2[ha] = Nodes_addelement(intNode, *PR_T2[0][i], Nodes2[ha]) ;
								*tNodes2 = Nodes_addelement(intNode, *PR_T2[0][i], *tNodes2) ;
								//Rprintf("New coalescent event between %i and %i\n", ch[0], ch[1]) ;
								//Rprintf("bt_n=%i", bt_n) ;
									
								if (bt_n >= 0)
								{
									tr_old = phylo_replacenode(tr_old, ((bt_n+1)*2)-1, ch[0], intNode) ;
									//Rprintf("replacing %i with %i\n", ch[0], intNode) ;
								}
								mig2 = Migmx_replacenode(mig2_n, mig2_N, mig2, ch[0], intNode) ;
								//Rprintf("old migration mx after replacing %i with %i\n", ch[0], intNode) ;
								for (j=0 ; j<mig2_N ; j++)
								{
									//Rprintf("%8.4f %8.4f %8.4f %8.4f %i\n", mig2[j], mig2[j + mig2_N], mig2[j + 2*mig2_N], mig2[j + 3*mig2_N] , mig2_n) ;
								}
								//OBS2[ha]-- ;
							}
							tr = phylo_addnewedge(tr, ei, intNode, ch[0], ttmp) ;
							ei-- ;
							tr = phylo_addnewedge(tr, ei, intNode, ch[1], ttmp2) ;
							ei-- ;
							
							intNode -- ;
							Anc[hb] -- ;
							ce++ ;
							//OBS[hb] -- ;
							
							/*
							Rprintf("old tree\n") ;
							if (sumOBS > 0)
							{
								for (j=0 ; j<tr_old->Nedge ; j++)
								{
									Rprintf("%i %i\n", tr_old->edge[j][0], tr_old->edge[j][1]) ;
								}
								
							}
							Rprintf("new tree\n") ;
							for (j=0 ; j<tr->Nedge ; j++)
							{
								Rprintf("%i %i\n", tr->edge[j][0], tr->edge[j][1]) ;
							}*/
							
						}
				  }  
        	}	
        	//llist_print(Nodes[ha].n) ;
			//Rprintf(" ha\n") ;
			//llist_print(Nodes[hb].n) ;
			//Rprintf(" hb\n") ;
        		// move in new nodes	
        		//Rprintf("Moving %i nodes\n", tNodes->N) ;
        		//Rprintf("%i nodes in ha, %i nodes in tNodes\n", Nodes[ha].N, tNodes->N) ;
        		
        		Nodes[ha] =  Nodes_copy_all(Nodes[ha], *tNodes)	;			
        			
        		*tNodes =  Nodes_delete_all(*tNodes) ;
        		Nodes2[ha] =  Nodes_copy_all(Nodes2[ha], *tNodes2)	;				
        		*tNodes2 =  Nodes_delete_all(*tNodes2) ;
        		
        		
				
			}
			if (ll_calc == 1)
			{
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
			}
			
			ll += ll_i ;
			}
			else
			{
				//Rprintf("end \n") ;
			}
			
			
		//	Rprintf("ll=%8.4f\t%8.4f\t%i\t%i\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, i, Anc[ha], Anc[hb]) ;
			
		
		}
		//Rprintf("ll=%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ll_i, ll, ha+1, hb+1, *PR_T2[0][i], *PR_T2[3][i], lambdao[hb][ha], lambdao_sum, ce, ce1, bt_n, Anc[0], Anc[1], Anc[2], Anc[3], Anc[4], Nodes[0].N, Nodes2[0].N, OBS[0], OBS2[0], Ancsum, intNode) ;
		//Rprintf("%i\t%i\t%8.4f\t%8.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", ha+1, hb+1, *PR_T2[0][i], *PR_T2[3][i], ce, ce1, bt_n, Nodes[0].N, Nodes2[0].N, Nodes[1].N, Nodes2[1].N, intNode) ;
			
		
}


/*******************************RESULT TO OUTPUT*****************************************/

	PROTECT(R_List = allocVector(VECSXP ,4)) ;
	PROTECT(R_bt=allocMatrix(REALSXP,tr->NNode,4)) ;
	PROTECT(R_ll = allocVector(REALSXP, 1)) ;
	PROTECT(R_mig = allocMatrix(REALSXP, mig_n-1, 4)) ;
	PR_bt = REAL(R_bt) ;
	
	/*
	bt = phylo_bt(tr, ST, OBS3, NHosts) ;
	for (j=0 ; j<tr->NNode ; j++)
	{
		for (i=0;i<4;i++)
		{
			PR_bt[j+i*tr->NNode]=(i>0) ? (bt[j][i] + 1) : bt[j][i] ;
			
		}
	}
	*/
	for (j=0 ; j<(mig_n - 1) ; j++)
	{
		for (i=0;i<4;i++)
		{
			REAL(R_mig)[j+i*(mig_n - 1)]=  mig[j][i] ;
			
		}
	}
	
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
		//Free(bt[j]) ;
		
	}
	//Free(bt) ;
	
	for (j=0 ; j<NHosts ; j++)
	{
		llist_destroy_i(Nodes[j].n)  ;
		llist_destroy_d(Nodes[j].t)  ;
		llist_destroy_i(Nodes2[j].n)  ;
		llist_destroy_d(Nodes2[j].t)  ;
		free(lambdao[j]) ;
	}
	llist_destroy_i(tNodes->n)  ;
	llist_destroy_d(tNodes->t)  ;
	llist_destroy_i(tNodes2->n)  ;
	llist_destroy_d(tNodes2->t)  ;
	free(Nodes) ; 
	free(Tend) ;  
	free(ch) ;
	free(PR_T2) ;
 	free(SNsum) ;
 	free(Anc) ;
 	free(lambdao) ;
 	free(vvec) ;
 	if (sumOBS > 0)
 	{
 		free(mig2) ;
	}
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
	Free(pp1) ;
	Free(pp2) ;
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
