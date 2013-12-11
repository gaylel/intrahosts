datadir=../RData 
testdir=../test


cSIR.so:	cSIR.c llists.c phylo.c leaf.c
	R CMD SHLIB cSIR.c llists.c phylo.c leaf.c

runmodel: cSIR.so
	./cSIRrun.sh 
	
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
	
				
#memcheck:	
#	R -d  "valgrind --tool=memcheck -v" --slave --vanilla --args "tree.1" < cSIRtest.R
	
		