#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "llists.h"

// leaf structure for ordering tips of the trees
struct leaf {
	int NLeaves ;
	int *lo ; // the permutation of the tips
	int *h ; // the host of each tip
	int *s ; // the order of coalescence events between neighbouring tips
} ;

typedef struct leaf leaf ;
struct leaf* leaf_create(int NLeaves) ;
struct leaf* R_to_leaf(SEXP R_lo, SEXP R_s) ;
