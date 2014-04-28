#include "beaglefuncs.h"

BeagleOperation* beagle_readtree(int **edge, double *el, int NNode)
{
	// convert tree to set of beagle operations
	// tree edges are ordered according to branching time
	int i, j=0 ;
	int NEdges = NNode * 2 ;
	int N = NEdges + NNode + 1;
	int c1, c2 ;
	// create Operation here?
	BeagleOperation * b_op = calloc(NNode, sizeof(BeagleOperation)) ;
	
	for (i=NEdges-1 ; i>0 ; i=i-2)
	{
		b_op[j].destinationPartials = N - edge[i][0] ;
		b_op[j].destinationScaleWrite = BEAGLE_OP_NONE ;
		b_op[j].destinationScaleRead = BEAGLE_OP_NONE ;
		c1 = (edge[i][1] > NNode) ? N-edge[i][1] : edge[i][1] ;
		c2 = (edge[i-1][1] > NNode) ? N-edge[i-1][1] : edge[i-1][1] ;
		b_op[j].child1Partials = c1 ;
		b_op[j].child1TransitionMatrix = c1 ;
		b_op[j].child2Partials = c2 ;
		b_op[j].child2TransitionMatrix = c2 ;
		//Rprintf("%i %i %i %i %i\n", b_op[j].destinationPartials, b_op[j].child1Partials, b_op[j].child1TransitionMatrix, b_op[j].child2Partials, b_op[j].child2TransitionMatrix) ;
		j++ ;
	}

	return b_op ;
}

int * beagle_nodeIndices_init(int NInd)
{
	int *nodeInd = calloc(NInd, sizeof(int)) ;
	int i ;
	for (i=0 ; i<NInd ; i++)
	{
		nodeInd[i] = i ;
	}
	return nodeInd ;
}

double * beagle_edgelengths_init(int NInd, double *el, double mu, int **edge)
{
	int NNode = NInd / 2; 
	int i ;
	int N = 3 * NNode + 1 ;
	int NEdges = NInd ;
	int c ;
	double * el_out = calloc(NInd, sizeof(double)) ;
	for (i=0 ; i<NEdges ; i++)
	{
		c = (edge[i][1] > NNode) ? N-edge[i][1] : edge[i][1] ;
		el_out[c] = el[i] * mu; 
	}
	return el_out ;
	
}

int * beagle_create_instances(int NSeqs, int NSites, int NNode, int Np)
{
	BeagleInstanceDetails* returnInfo ;
	returnInfo = (BeagleInstanceDetails*)calloc(1, sizeof(BeagleInstanceDetails)) ;
	// create an instance of the BEAGLE library
	int i, *instance;
	instance = calloc(Np, sizeof(int)) ;
	for (i=0 ; i<Np ; i++)
	{
		instance[i]	 = beagleCreateInstance(
				NSeqs,		/**< Number of tip data elements (input) */
				NSeqs-1,	        /**< Number of partials buffers to create (input) -- internal node count */
				NSeqs,		/**< Number of compact state representation buffers to create -- for use with setTipStates (input) */
				4,		/**< Number of states in the continuous-time Markov chain (input) -- DNA */
				NSites,	/**< Number of site patterns to be handled by the instance (input) -- not compressed in this case */
				1,		/**< Number of eigen-decomposition buffers to allocate (input) */
				NNode * 2,		/**< Number of transition matrix buffers (input) -- one per edge */
				1,		/**< Number of rate categories */
				0,		/**< Number of scaling buffers -- can be zero if scaling is not needed*/
				NULL,		/**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
				0,		/**< Length of resourceList list (input) -- not needed to use the default hardware config */
				0,		/**< Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input) */
				0,		/**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
			    returnInfo
				);
	}
	return instance ;
}

void beagle_free_instances(int Np, int *instances)
{
	int i ;
	for (i=0 ; i<Np ; i++)
	{
		beagleFinalizeInstance(instances[i]);
	}
}

int beagle_initialize_instance(int instance, int NSeqs, int NSites, int** seqs, int* w)
{
	int i , er;
	
	for (i=0 ; i< NSeqs ; i++)
	{
		er = beagleSetTipStates(instance, i, seqs[i] ) ;
	}
	
	// equal weighting for each site.
	double * pattern_weights=calloc(NSites, sizeof(double)) ;
	for (i=0 ; i< NSites ; i++)
	{
		pattern_weights[i] = w[i] ;
	//	printf("%i ", w[i]) ;
	}
	er =beagleSetPatternWeights(instance, pattern_weights) ;
	
	// set state background frequencies
	double freqs[4] = {0.25, 0.25, 0.25, 0.25} ;
	er = beagleSetStateFrequencies(instance, 0, freqs) ;
	
	
	//create an array containing site category weights and rates
	const double weights[1] = { 1.0 } ;
	const double rates[1] = { 1.0 } ;
	er = beagleSetCategoryWeights(instance, 0, weights) ;
	er = beagleSetCategoryRates(instance, rates) ;
	
	 // an eigen decomposition for the JC69 model
    double evec[4 * 4] = {
                 1.0,  2.0,  0.0,  0.5,
                 1.0,  -2.0,  0.5,  0.0,
                 1.0,  2.0, 0.0,  -0.5,
                 1.0,  -2.0,  -0.5,  0.0
        };

    double ivec[4 * 4] = {
                 0.25,  0.25,  0.25,  0.25,
                 0.125,  -0.125,  0.125,  -0.125,
                 0.0,  1.0,  0.0,  -1.0,
                 1.0,  0.0,  -1.0,  0.0
        };

    double eval[4] = { 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333 };

	// set the Eigen decomposition
    er = beagleSetEigenDecomposition(instance, 0, evec, ivec, eval);
	
	free(pattern_weights) ;
	return 1 ; 
}

double beagle_update_instance(int instance, int NTips, int **edge, double *el, double mu, int NNode)
{
	// NNode : number of internal nodes of (sub) tree
	// NTips : number of tips of full tree

	int i , j, k, er;
	
        // a list of indices and edge lengths
        // these get used to tell beagle which edge length goes with which node
        //int  nodeIndices[4] = { 0, 1, 2, 3 };
        int *nodeIndices = beagle_nodeIndices_init(NNode * 2) ;
        
        for (i=NNode+1 ; i<NNode*2 ; i++)
        {
        	nodeIndices[i] += NTips - NNode - 1 ; 
        }
        
        double *edgeLengths = beagle_edgelengths_init(NNode * 2, el, mu, edge) ;
	//	printf("%8.4f %8.4f %8.4f %8.4f \n", edgeLengths[0], edgeLengths[1], edgeLengths[2], edgeLengths[3] ) ;
        // tell BEAGLE to populate the transition matrices for the above edge lengthss
        er = beagleUpdateTransitionMatrices(instance,     // instance
                                 0,             // eigenIndex
                                 nodeIndices,   // probabilityIndices
                                 NULL,          // firstDerivativeIndices
                                 NULL,          // secondDervativeIndices
                                 edgeLengths,   // edgeLengths
                                 NNode * 2);            // count
         for (i=0 ; i<NNode * 2 ; i++)
        {
        	//Rprintf("node index %i, %8.4f\n", nodeIndices[i], edgeLengths[i]) ;
        }
        
        /*double *tm =calloc(16, sizeof(double)) ;
        
        for (k=0 ; k<NNode * 2 ; k++)
        {
        er = beagleGetTransitionMatrix(instance, k,  tm) ;
        
        for (i=0 ; i<4 ; i++)
        {
        	for (j=0 ; j<4 ; j++)
        	Rprintf("%8.4f ",tm[i*4 + j]) ;
        	Rprintf("\n") ; 
        }         
        Rprintf("\n") ;
        }
        free(tm) ;
        */
        // create a list of partial likelihood update operations
        // the order is [dest, sourceScaling, destScaling, source1, matrix1, source2, matrix2]
        // these operations say: first peel node 0 and 1 to calculate the per-site partial likelihoods, and store them
        // in buffer 3.  Then peel node 2 and buffer 3 and store the per-site partial likelihoods in buffer 4.
        /*BeagleOperation operations[2] = {
                {3, BEAGLE_OP_NONE, BEAGLE_OP_NONE, 0, 0, 1, 1},
                {4, BEAGLE_OP_NONE, BEAGLE_OP_NONE, 2, 2, 3, 3}
        };*/
        BeagleOperation  * operations = beagle_readtree(edge, el, NNode) ;
        int offset = NTips - NNode - 1; 
        for (i=0 ; i<NNode ; i++)
        {
        	if (operations[i].destinationPartials > NNode)
        			operations[i].destinationPartials += offset ;
        	if (operations[i].child1Partials > NNode)
        	{
        		operations[i].child1Partials += offset ;
        		operations[i].child1TransitionMatrix += offset ;
        	}
        	if (operations[i].child2Partials > NNode)
        	{
        		operations[i].child2Partials += offset ;
        		operations[i].child2TransitionMatrix += offset ;
     		}
     		//Rprintf("Operations %i %i %i %i %i\n", operations[i].destinationPartials, operations[i].child1Partials, operations[i].child1TransitionMatrix, operations[i].child2Partials, operations[i].child2TransitionMatrix) ;
        }
		//Rprintf("read in tree %i\n", instance) ;
        // this invokes all the math to carry out the likelihood calculation
        
         for (j=1 ; j<2 ; j++)
         {
        	// beagleGetPartials(instance, j, BEAGLE_OP_NONE, tm) ;
        
        	//Rprintf("%8.4f ",tm[j]) ;
        }	
        //Rprintf("\n") ;
        er = beagleUpdatePartials( instance,      // instance
                        operations,     // eigenIndex
                        NNode,              // operationCount
                        BEAGLE_OP_NONE);             // cumulative scale index
		//Rprintf("Updated partials %i\n", instance) ;
        double logL = 0;
        int rootIndex[1] = {NTips + NNode - 1};
        //Rprintf("root index = %i\n", rootIndex[0]) ;
        int categoryWeightIndex[1] = {0};
        int stateFrequencyIndex[1] = {0};
        int cumulativeScaleIndex[1] = {BEAGLE_OP_NONE};

		
		// calculate the site likelihoods at the root node
        // this integrates the per-site root partial likelihoods across sites, background state frequencies, and rate categories
        // results in a single log likelihood, output here into logL
        beagleCalculateRootLogLikelihoods(instance,               // instance
                                    rootIndex,// bufferIndices
                                    categoryWeightIndex,                // weights
                                    stateFrequencyIndex,                 // stateFrequencies
                                    cumulativeScaleIndex,          // scaleBuffer to use
                                    1,                      // count
                                    &logL);         // outLogLikelihoods
		//Rprintf("calculated likelihoods %8.4f\n", logL) ;
	free(nodeIndices) ;
	free(edgeLengths) ;
	free(operations) ;
	 
	return logL ;
}



double beagle_init(int NSeqs, int NSites, int** seqs, int* w, int **edge, double *el, double mu, int NNode)
{
	int i , er;
	int NTips = NNode + 1 ;
	
	BeagleInstanceDetails* returnInfo ;
	returnInfo = (BeagleInstanceDetails*)calloc(1, sizeof(BeagleInstanceDetails)) ;
	// create an instance of the BEAGLE library
	int instance = beagleCreateInstance(
				NSeqs,		/**< Number of tip data elements (input) */
				NSeqs-1,	        /**< Number of partials buffers to create (input) -- internal node count */
				NSeqs,		/**< Number of compact state representation buffers to create -- for use with setTipStates (input) */
				4,		/**< Number of states in the continuous-time Markov chain (input) -- DNA */
				NSites,	/**< Number of site patterns to be handled by the instance (input) -- not compressed in this case */
				1,		/**< Number of eigen-decomposition buffers to allocate (input) */
				NNode * 2,		/**< Number of transition matrix buffers (input) -- one per edge */
				1,		/**< Number of rate categories */
				0,		/**< Number of scaling buffers -- can be zero if scaling is not needed*/
				NULL,		/**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
				1,		/**< Length of resourceList list (input) -- not needed to use the default hardware config */
				0,		/**< Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input) */
				0,		/**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
				returnInfo
				);

	if (instance < 0) {
		//Rprintf("Failed to obtain BEAGLE instance\n\n");
	}
	for (i=0 ; i< NSeqs ; i++)
	{
		er = beagleSetTipStates(instance, i, seqs[i] ) ;
	}
	
	// equal weighting for each site.
	double * pattern_weights=calloc(NSites, sizeof(double)) ;
	for (i=0 ; i< NSites ; i++)
	{
		pattern_weights[i] = w[i] ;
	//	printf("%i ", w[i]) ;
	}
	er =beagleSetPatternWeights(instance, pattern_weights) ;
	//printf("weights set %i\n", er) ;
	// set state background frequencies
	double freqs[4] = {0.25, 0.25, 0.25, 0.25} ;
	er = beagleSetStateFrequencies(instance, 0, freqs) ;
	//printf("Frequencies set %i\n", er) ;
	
	//create an array containing site category weights and rates
	const double weights[1] = { 1.0 } ;
	const double rates[1] = { 1.0 } ;
	er = beagleSetCategoryWeights(instance, 0, weights) ;
	//printf("category weights set %i\n", er) ;
	
	
	er = beagleSetCategoryRates(instance, rates) ;
	//printf("category rates set %i\n", er) ;
	
	 // an eigen decomposition for the JC69 model
        double evec[4 * 4] = {
                 1.0,  2.0,  0.0,  0.5,
                 1.0,  -2.0,  0.5,  0.0,
                 1.0,  2.0, 0.0,  -0.5,
                 1.0,  -2.0,  -0.5,  0.0
        };

        double ivec[4 * 4] = {
                 0.25,  0.25,  0.25,  0.25,
                 0.125,  -0.125,  0.125,  -0.125,
                 0.0,  1.0,  0.0,  -1.0,
                 1.0,  0.0,  -1.0,  0.0
        };

        double eval[4] = { 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333 };



        // set the Eigen decomposition
        er = beagleSetEigenDecomposition(instance, 0, evec, ivec, eval);
	//	printf("buffer set %i\n", er) ;
	
        // a list of indices and edge lengths
        // these get used to tell beagle which edge length goes with which node
        //int  nodeIndices[4] = { 0, 1, 2, 3 };
        int *nodeIndices = beagle_nodeIndices_init(NNode * 2) ;
        //double edgeLengths[4] = { 0.1, 0.1, 0.1, 0.1 };
		double *edgeLengths = beagle_edgelengths_init(NNode * 2, el, mu, edge) ;
	//	printf("%8.4f %8.4f %8.4f %8.4f \n", edgeLengths[0], edgeLengths[1], edgeLengths[2], edgeLengths[3] ) ;
        // tell BEAGLE to populate the transition matrices for the above edge lengthss
        er = beagleUpdateTransitionMatrices(instance,     // instance
                                 0,             // eigenIndex
                                 nodeIndices,   // probabilityIndices
                                 NULL,          // firstDerivativeIndices
                                 NULL,          // secondDervativeIndices
                                 edgeLengths,   // edgeLengths
                                 NNode * 2);            // count
                                 
      //  printf("tm set %i\n", er) ;
	                         
	//	printf("got here\n") ;
        // create a list of partial likelihood update operations
        // the order is [dest, sourceScaling, destScaling, source1, matrix1, source2, matrix2]
        // these operations say: first peel node 0 and 1 to calculate the per-site partial likelihoods, and store them
        // in buffer 3.  Then peel node 2 and buffer 3 and store the per-site partial likelihoods in buffer 4.
        /*BeagleOperation operations[2] = {
                {3, BEAGLE_OP_NONE, BEAGLE_OP_NONE, 0, 0, 1, 1},
                {4, BEAGLE_OP_NONE, BEAGLE_OP_NONE, 2, 2, 3, 3}
        };*/
        BeagleOperation  * operations = beagle_readtree(edge, el, NNode) ;
         for (i=0 ; i<NNode ; i++)
        {
      //  	printf("%i %i %i %i %i %i %i\n", operations[i].destinationPartials, operations[i].destinationScaleWrite, operations[i].destinationScaleRead , operations[i].child1Partials, operations[i].child1TransitionMatrix,  operations[i].child2Partials, operations[i].child2TransitionMatrix ) ;
        }
		//printf("got here %i\n", BEAGLE_OP_NONE) ;
        // this invokes all the math to carry out the likelihood calculation
        er = beagleUpdatePartials( instance,      // instance
                        operations,     // eigenIndex
                        NNode,              // operationCount
                        BEAGLE_OP_NONE);             // cumulative scale index

		//printf("got here %i\n" ,er) ;
        double logL = 0;
        int rootIndex[1] = {NNode * 2};
        int categoryWeightIndex[1] = {0};
        int stateFrequencyIndex[1] = {0};
        int cumulativeScaleIndex[1] = {BEAGLE_OP_NONE};

		
		// calculate the site likelihoods at the root node
        // this integrates the per-site root partial likelihoods across sites, background state frequencies, and rate categories
        // results in a single log likelihood, output here into logL
        beagleCalculateRootLogLikelihoods(instance,               // instance
                                    rootIndex,// bufferIndices
                                    categoryWeightIndex,                // weights
                                    stateFrequencyIndex,                 // stateFrequencies
                                    cumulativeScaleIndex,          // scaleBuffer to use
                                    1,                      // count
                                    &logL);         // outLogLikelihoods

		/*double outLogLikelihoods[NSites] ;
		double oll = 0 ;
		beagleGetSiteLogLikelihoods(instance, outLogLikelihoods);
		for (i = 0 ; i<NSites ; i++)
		{
			oll += w[i]*outLogLikelihoods[i];
			printf("%8.4f\n", outLogLikelihoods[i]) ;
		}*/
	//printf("%8.4f\n", logL) ;
	int err_bfi ;
	free(pattern_weights) ;
	free(nodeIndices) ;
	free(edgeLengths) ;
	free(operations) ;
	//free(rsrcList->list) ;
	//free(rsrcList) ; 
	//beagleFinalize() ;
	//BeagleResourceList* br_list = beagleGetResourceList() ;
	//Rprintf("resource list, length %i, at %p, %p\n", br_list->length, br_list[0], br_list->list) ;
	err_bfi = beagleFinalizeInstance(instance);
	//instance = -1 ;
	//Rprintf("%i\n", err_bfi) ;
	//for (i=0 ; i<br_list->length ; i++)
	//free(br_list[i].list) ;
	//free(br_list) ;
	
	//Rprintf("resource list, length %i, at %p, %p\n", br_list->length, br_list[0], br_list->list) ;
  	free(returnInfo) ;
  	return logL ;
		
}