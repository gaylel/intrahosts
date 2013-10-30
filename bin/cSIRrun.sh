#!/bin/bash


outdir=../test
#outdir=../$(date +%Y%m%d) ;
mkdir $outdir ;
datadir="../RData" ;
params="test.params"
R --slave --vanilla --args $outdir $datadir $outdir/$params < cSIRrun.R