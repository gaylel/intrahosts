testdir=../test
debugdir=../debug 
datadir=../../RData 


cSIR.so:	cSIR.c llists.c phylo.c leaf.c
	R CMD SHLIB cSIR.c llists.c phylo.c leaf.c

debugdir:
	mkdir $(debugdir) 
	
debug:	cSIR.so debugdir
	R -d "valgrind --leak-check=full --show-reachable=yes" --vanilla --args $(debugdir) $(datadir) < cSIRrun.R	

testdir:
	if [ ! -d "$(testdir)" ]; then \
	mkdir $(testdir) ;\
	fi 
	
test1:	testdir cSIR.so
	#R -d  gdb --vanilla --args "tree.1" < cSIRtest.R
	R --slave --vanilla --args "tree.1" $(testdir) < cSIRtest.R
	
test2:	testdir cSIR.so
	#R -d  gdb --vanilla --args "tree.1" < cSIRtest.R
	R --slave --vanilla --args "tree.2" $(testdir) < cSIRtest.R	
	
test3:	testdir cSIR.so
	#R -d  gdb --vanilla --args "tree.1" < cSIRtest.R
	R --slave --vanilla --args "tree.3" $(testdir) < cSIRtest.R		

test4:	testdir cSIR.so
	#R -d  gdb --vanilla --args "tree.1" < cSIRtest.R
	R --slave --vanilla --args "tree.4" $(testdir) < cSIRtest.R		
	
test5:	testdir cSIR.so
	#R -d  gdb --vanilla --args "tree.1" < cSIRtest.R
	R --slave --vanilla --args "tree.5" $(testdir) < cSIRtest.R		
	
testB:	testdir cSIR.so
	R --slave --vanilla --args "tree.B" $(testdir) < cSIRtest.R		
	
testBPlot: 
	./mcmcplot.sh $(testdir)	
	
simtest:	testdir cSIR.so
	./cSIRsimdata.sh "simtest"	
	
runsimtest:	testdir cSIR.so
	./cSIRrunsim.sh "simtest"	

runsimtestplot:
	./mcmcplot.sh ../simtest/20131014

runmodel:
	./cSIRrun.sh 2
	
runmodelplot:
	./mcmcplot.sh ../run_2

memcheck:	
	R -d  "valgrind --tool=memcheck -v" --slave --vanilla --args "tree.1" < cSIRtest.R
	
		