args <- commandArgs()
source("~/Documents/Projects/R/SIR/plots.R")
# read in skyride plot
sr<-read.table(args[6],skip=2) ;
library("ape") ;
ph<-read.nexus(args[5]) ;
skyride_plots(ph,sr,args[7]) ;

