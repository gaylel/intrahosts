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
<<<<<<< HEAD


=======
>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714
phylo* phylo_create(int ntips) ;
SEXP phylo_to_R(phylo* tr) ;
double** phylo_bt(phylo* tr, double *tipinfo, int *hostinfo, int NHosts) ;
