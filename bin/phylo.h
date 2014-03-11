#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

struct phylo {
	int **edge ;
	double *el ;
	int *nodelabel ;
	int NNode ;
	int *tiplabel ;	
	int Nedge ;
	
} ;

typedef struct phylo phylo ;
phylo* phylo_create(int ntips) ;
SEXP phylo_to_R(phylo* tr) ;
double** phylo_bt(phylo* tr, double *tipinfo, int *hostinfo, int NHosts) ;
phylo * R_to_phylo(SEXP R_tr) ;
phylo * phylo_addnewedge(phylo * tr, int ei, int from_node, int to_node, double elength) ;
phylo * phylo_replacenode(phylo* tr, int n, int oldnode, int newnode) ;
