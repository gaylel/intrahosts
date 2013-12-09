#!/bin/bash

# plot the traces of the variables 

smpdir=$1
outdir=$smpdir
paramfile=$smpdir/test.params
R --slave --vanilla --args $smpdir $outdir $paramfile  < mcmcplot.R