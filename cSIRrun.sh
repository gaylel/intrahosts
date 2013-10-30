#!/bin/bash


outdir=../run_$1
#outdir=../$(date +%Y%m%d) ;
mkdir $outdir ;
datadir="../../RData" ;
params="test.params"
R --slave --vanilla --args $outdir $datadir $outdir/$params < cSIRrun.R