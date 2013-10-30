#!/bin/bash

datafile=../data/newmarket/HE.fasta ;
less $datafile | sed "s/.*isolate \([0-9]*\).*clone \([0-9]*\)/>I\1S\2/" > ../data/newmarket/IS.fasta ;

# ape on mac can't handle capital G - 
less ../data/newmarket/IS.fasta | sed 's/G/g/g' > ../data/newmarket/IS_mac.fasta 
