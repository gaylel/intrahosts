#!/bin/bash

fname="test_2609"
<<<<<<< HEAD
outdir=${fname}/$(date +%Y%m%d) ;
=======
outdir=${fname}/$(date +%Y%m%d)_2 ;
>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714
mkdir $outdir ;
datadir="${fname}/RData" ;
R --slave --vanilla --args $outdir $datadir < cSIRrunsim.R
