#!/bin/bash

datafile=../data/newmarket/HE.fasta ;
less $datafile | sed "s/.*isolate \([0-9]*\).*clone \([0-9]*\)/>I\1S\2/" > ../data/newmarket/IS.fasta ;
