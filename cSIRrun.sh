#!/bin/bash


outdir=../run_$1
#outdir=../$(date +%Y%m%d) ;
mkdir $outdir ;
datadir="../../RData" ;
R --slave --vanilla --args $outdir $datadir < cSIRrun.R