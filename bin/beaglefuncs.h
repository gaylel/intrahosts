#include <libhmsbeagle/beagle.h>
#include <stdlib.h>
#include <time.h>


int * beagle_nodeIndices_init(int NInd) ;
double * beagle_edgelengths_init(int NInd, double *el, double mu, int **edge) ;
double beagle_init(int NSeqs, int NSites, int** seqs, int* w, int **edge, double *el, double mu, int NNode) ;
int* beagle_init_instances(int NSeqs, int NSites, int NNode, int Np) ;
void beagle_free_instances(int Np, int *instances) ;
double beagle_update_instance(int instance, int NSeqs, int NSites, int** seqs, int* w, int **edge, double *el, double mu, int NNode) ;
BeagleOperation* beagle_readtree(int **edge, double *el, int NNode) ;