#!/bin/bash

fname=$1

outdir=../${fname}/$(date +%Y%m%d) ;
mkdir $outdir ;
datadir="../${fname}/RData" ;
R --slave --vanilla --args $outdir $datadir < cSIRrunsim.R
