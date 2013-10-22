#!/bin/bash

# plot the traces of the variables 

smpdir=$1
outdir=$smpdir
R --slave --vanilla --args $smpdir $outdir  < mcmcplot.R