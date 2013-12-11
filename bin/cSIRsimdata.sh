#!/bin/bash

fname=$1


rundir=${fname}
datadir=$rundir/data
mkdir $datadir

R --slave --vanilla --args $rundir < simTree.R
seq-gen -mGTR -s1.83e-4 -l903 < ${datadir}/test.nwk > ${datadir}/test.dat
