#include "phylo.h"
#include <stdlib.h>

phylo* phylo_create(int ntips)
{
	int Nedges = (ntips-1)*2, m, n ;
	phylo *tr ;
	tr=Calloc(1,phylo);
	tr->Nedge=Nedges ;
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
	return tr ;
}


double** phylo_bt(phylo* tr, double *tipinfo, int *hostinfo, int NHosts)
{
	// returning branching times for a phylo structure where the edges are ordered by intNode indices
	// tipinfo[i] : time at host i
	// hostinfo[i] : number of tips in host i
	int NTips=tr->NNode+1, i, j, k, Ni, *hosts, ch1, ch2 ;
	double **bt, *tips ;
	// create structures
	hosts=Calloc(NTips+tr->NNode,int) ;
	tips=Calloc(NTips+tr->NNode,double) ;
	k=0 ;
	for (i=0 ; i<NHosts ; i++)
	{
		for (j=0 ; j<hostinfo[i] ; j++)
		{
			tips[k]=tipinfo[i] ;
			hosts[k++]=i ;
		}
	}
	
	for (i=0 ; i<tr->NNode ; i++)
	{
		tips[k]=-1 ;
		hosts[k++]=tr->nodelabel[i] ;
	}

	bt=Calloc(tr->NNode, double*) ;
	for (i=0 ; i< tr->NNode ; i++)
	{
		bt[i]=Calloc(4, double) ;
	}
	
	for (i=tr->NNode-1 ; i>=0 ; i--)
	{
		Ni=NTips+i ;
		j=2*i ;
		ch1=tr->edge[j][1] ;
		ch2=tr->edge[j+1][1] ;
		//Rprintf("%i %i %i\n",Ni, ch1,ch2) ; 
		
		bt[i][0]=tips[ch1] - tr->el[j] ;
		tips[Ni]=bt[i][0] ;
		bt[i][1]=hosts[Ni] ;
		bt[i][2]=hosts[ch1] ;
		bt[i][3]=hosts[ch2] ;
		//Rprintf("%8.4f %i %i %i\n", bt[i][0], hosts[Ni], hosts[ch1], hosts[ch2]) ; 
		
	}
	Free(hosts) ;
	Free(tips) ;
	return bt ;
}

phylo * R_to_phylo(SEXP R_tr)
{
	int ntips = INTEGER(VECTOR_ELT(R_tr, 4))[0] + 1;
	SEXP nl;
	//char **ch;
	//Rprintf("ntips = %i\n", ntips) ;
	phylo* tr = phylo_create(ntips) ;
	int nedges = 2 * (ntips - 1) ;
	int i ;
	for (i=0 ; i<(nedges) ; i++)
	{
		tr->edge[i][0] = INTEGER(VECTOR_ELT(R_tr, 0))[i] - 1 ;
		tr->edge[i][1] = INTEGER(VECTOR_ELT(R_tr, 0))[nedges + i] - 1 ;
		//Rprintf("edge = %i %i\n", tr->edge[i][0], tr->edge[i][1]) ;
	
		tr->el[i] = REAL(VECTOR_ELT(R_tr, 1))[i] ;
	}
	
	nl = VECTOR_ELT(R_tr, 3) ;
	for (i=0 ; i<tr->NNode ; i++)
	{
	//	Pmychar[0] = R_alloc(strlen(CHAR(STRING_ELT(mychar, 0))), sizeof(char));
    
		
		tr->nodelabel[i] = atoi(CHAR(STRING_ELT(nl, i)) ) ;
	}
	return tr ;
}


SEXP phylo_to_R(phylo* tr)
{
	SEXP R_edge, R_edgelength, R_NNode, R_nodelabel, R_tiplabel, R_list, class, R_names ;
	int *PR_edge, *PR_nodelabel, *PR_tiplabel, m, n, c ;
	double *PR_edgelength ;
	char ch[20] ;
	PROTECT(R_edge=allocMatrix(INTSXP,tr->Nedge,2)) ;
	PROTECT(R_edgelength=allocVector(REALSXP, tr->Nedge)) ;
	PROTECT(R_nodelabel=allocVector(STRSXP,tr->NNode)) ;
	PROTECT(R_NNode = allocVector(INTSXP,1)) ;
	PROTECT(R_tiplabel=allocVector(STRSXP, tr->NNode+1)) ;
	PROTECT(R_list = allocVector(VECSXP ,5)) ;
	PR_edge = INTEGER(R_edge) ;
	PR_edgelength = REAL(R_edgelength) ;
	//PR_tiplabel = STRING_PTR(R_tiplabel) ;
	//PR_nodelabel = CHAR(R_nodelabel) ;
	for (n=0 ; n<2 ; n++)
	{
		c=n*tr->Nedge ;
		for (m=0 ; m<tr->Nedge ; m++)
		{
			PR_edge[m+c]=tr->edge[m][n] + 1;
		}
	}
	
	for (m=0 ; m<tr->Nedge ; m++)
	{
		PR_edgelength[m]=tr->el[m] ;
	}
	
	for (m=0 ; m<tr->NNode ; m++)
	{
		//PR_nodelabel[m]=mkChar(tr->nodelabel[m] + 1);
		sprintf(ch, "%d",tr->nodelabel[m]+1 ) ;
		SET_STRING_ELT(R_nodelabel, m, mkChar(ch));
	}
	
	for (m=0 ; m<tr->NNode+1 ; m++)
	{
		//PR_tiplabel[m]=mkChar(tr->tiplabel[m]+1 );
		sprintf(ch, "%d",tr->tiplabel[m]+1 ) ;
		SET_STRING_ELT(R_tiplabel,m,mkChar(ch));
	}
	//SET_VECTOR_ELT(R_NNode, 0, tr->NNode) ;
	INTEGER(R_NNode)[0] = tr->NNode ;
	SET_VECTOR_ELT(R_list, 0, R_edge) ;
	SET_VECTOR_ELT(R_list, 1, R_edgelength) ;
	SET_VECTOR_ELT(R_list, 2, R_tiplabel) ;
	SET_VECTOR_ELT(R_list, 3, R_nodelabel) ;
	SET_VECTOR_ELT(R_list, 4, R_NNode) ;
	PROTECT(R_names=allocVector(STRSXP,5)) ;
	SET_STRING_ELT(R_names, 0, mkChar("edge")) ;
	SET_STRING_ELT(R_names, 1, mkChar("edge.length")) ;
	SET_STRING_ELT(R_names, 2, mkChar("tip.label")) ;
	SET_STRING_ELT(R_names, 3, mkChar("node.label")) ;
	SET_STRING_ELT(R_names, 4, mkChar("Nnode")) ;
	setAttrib(R_list, R_NamesSymbol, R_names);
	PROTECT(class = allocVector(STRSXP,1)) ;
	SET_STRING_ELT(class, 0, mkChar("phylo")) ;
	
	classgets(R_list,class) ;
	UNPROTECT(8) ;
	free(tr->el) ;
	free(tr->edge) ;
	free(tr->tiplabel) ;
	free(tr->nodelabel) ;
	free(tr) ;
	
	
	return R_list ;
}

phylo * phylo_addnewedge(phylo * tr, int ei, int from_node, int to_node, double elength)
{
	//Rprintf("%i %i\n", from_node, to_node) ;
	tr->edge[ei][0] = from_node ;
	tr->edge[ei][1] = to_node ;
	tr->el[ei] = elength ;
							
	return tr ;
}

phylo * phylo_replacenode(phylo* tr, int n, int oldnode, int newnode)
{
	// n : index of edge row
	int j ;
	
	if (n > 0)
	{
		// replace intNode in tr_old
		for (j = n ; j>=0 ; j--)
		{
			if (tr->edge[j][0]==oldnode)
			{
				tr->edge[j][0] = newnode ;			
			}
			if (tr->edge[j][1]==oldnode)
			{
				tr->edge[j][1] = newnode ;
			}
			
		}
	}
	return tr ;
}