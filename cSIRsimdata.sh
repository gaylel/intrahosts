#!/bin/bash

#fname=$1

fname="test_2609" 
mkdir $fname
cd $fname
datadir="Rdata"
mkdir $datadir ;
cd ..
R --slave --vanilla --args ${fname}/$datadir < simTree.R
seq-gen -mGTR -s1.83e-4 -l903 < ${fname}/${datadir}/test.nwk > ${fname}/${datadir}/test.dat
