#!/bin/bash

fname="test_2609"
outdir=${fname}/$(date +%Y%m%d) ;
mkdir $outdir ;
datadir="${fname}/RData" ;
R --slave --vanilla --args $outdir $datadir < cSIRrunsim.R
