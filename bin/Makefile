datadir=../RData 
testdir=../test

cSIR.so:	cSIR.c llists.c phylo.c leaf.c
	R CMD SHLIB cSIR.c llists.c phylo.c leaf.c

runmodel: cSIR.so
	./cSIRrun.sh 
	
runmodelplot:
	./mcmcplot.sh $(testdir)

#memcheck:	
#	R -d  "valgrind --tool=memcheck -v" --slave --vanilla --args "tree.1" < cSIRtest.R
	
		