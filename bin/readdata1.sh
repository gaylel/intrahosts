#!/bin/bash

fname="/Users/gayle/Documents/PHE/data/newmarket/journal.ppat.1003081.s001.txt" 
datadir="../Rdata"
mkdir $datadir ;
R --slave --vanilla --args $fname $datadir < readindata1.R