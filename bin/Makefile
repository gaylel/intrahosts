datadir=../RData 
testdir=../test


cSIR.so:	cSIR.c llists.c phylo.c leaf.c
	R CMD SHLIB cSIR.c llists.c phylo.c leaf.c

#runmodel: cSIR.so
#	./cSIRrun.sh 
	
runmodelplot:
	./mcmcplot.sh $(testdir)

.PHONY : simdir
simdir: 
	if [ ! -f ../sim ] ; then \
		mkdir ../sim ;\
		cp ../params/sim.params ../sim/sim.params ;\
	fi
	
simdata:	simdir
	./cSIRsimdata.sh ../sim ;	
	R --slave --vanilla --args ../sim < simData.R	

simtest.%:	
	mkdir ../sim/$@
	cat  ../params/test.params | \
	sed 's/firstN.*/firstN=seq(1,$*),/' | \
	sed 's/Niters.*/Niters=1000,/'  > ../sim/$@/test.params

simtime.%:	simtest.%
	R --slave --vanilla --args ../sim/simtest.$* ../sim/data ../sim/simtest.$*/test.params < cSIRtime.R    \
	> ../sim/simtest.$*/simtime.$*.log 
	#sed 's/real \([0-9].*\)/$* \1/' 
	 
simtimeall:
	for i in $$(seq 1 6) ; do \
		make simtime.$$i ; \
	done 
	 
timetest:	simtimeall
	R --slave --vanilla --args ../sim < cSIRtimeplot.R    
	
testdir:
	mkdir ../test
	cat  ../params/sim.params | \
	sed 's/h <-.*/h <- 4/' > ../test/sim.params
	
testdata:	testdir
	./cSIRsimdata.sh ../test ;	
	R --slave --vanilla --args ../test < simData.R	

testrun.%:
	#mkdir ../test/$@
	#cat ../params/test.params | \
	#sed 's/Niters.*/Niters=500000,/'  | \
	#sed 's/firstN.*/firstN=seq(1,4),/' >../test/$@/test.params
	R --slave --vanilla --args ../test/$@ ../test/data ../test/$@/test.params < cSIRrun.R	
				
testrunplot.%:
	R --slave --vanilla --args ../test/testrun.$* ../test/testrun.$* ../test/testrun.$*/test.params  < mcmcplot.R
#memcheck:	
#	R -d  "valgrind --tool=memcheck -v" --slave --vanilla --args "tree.1" < cSIRtest.R

rundir:
	mkdir ../run
	
runparams.%:	
	mkdir ../run/runmodel.$*
	cat ../params/test.params | \
	sed 's/firstN.*/firstN=c(4, 11, 16, 18, 24),/' | \
	sed 's/NS=.*/NS=1000,/' | \
	sed 's/Niters.*/Niters=500000,/' > ../run/runmodel.$*/test.params
	
runmodel.%:	
	R --slave --vanilla --args ../run/$@ ../RData ../run/$@/test.params < cSIRrun.R	
	
plotrunmodel.%:
	R --slave --vanilla --args ../run/runmodel.$* ../run/runmodel.$* ../run/runmodel.$*/test.params  < mcmcplot.R

debugdir:
	mkdir ../debug
	cat ../params/sim.params | \
	sed 's/h <-.*/h <- 5/' | \
	sed 's/nseq <-.*/nseq <- 10/' > ../debug/sim.params 
	
debugdata:
	./cSIRsimdata.sh ../debug ;	
	R --slave --vanilla --args ../debug < simData.R	

debugrun.%:
	R --slave --vanilla --args $* ../debug  < testmodel.R			