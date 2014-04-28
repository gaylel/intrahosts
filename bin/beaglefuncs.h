#include <libhmsbeagle/beagle.h>
#include <stdlib.h>
#include <time.h>


struct smc_instance {
	int key ;		// instance number
	
	
} ;

typedef struct smc_instance smc_instance ;

int * beagle_nodeIndices_init(int NInd) ;
double * beagle_edgelengths_init(int NInd, double *el, double mu, int **edge) ;
double beagle_init(int NSeqs, int NSites, int** seqs, int* w, int **edge, double *el, double mu, int NNode) ;
int* beagle_create_instances(int NSeqs, int NSites, int NNode, int Np) ;
int beagle_initialize_instance(int instance, int NSeqs, int NSites, int** seqs, int* w) ;
double beagle_update_instance(int instance, int NTips, int **edge, double *el, double mu, int NNode) ;
void beagle_free_instances(int Np, int *instances) ;
BeagleOperation* beagle_readtree(int **edge, double *el, int NNode) ;