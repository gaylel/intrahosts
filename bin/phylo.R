s_lo_from_phylo<-function(tr)
{
    # given a phylo object return the ordering of the tips and the coalescent ordering
    
    t<-reorder(tr,"cladewise")
    ntips<-length(t$tip.label)
    tips<-t$edge[,2]<=length(t$tip.label)
    lo<-t$edge[tips,2]
    s<-NULL
    lo2<-lo
    for (j in seq((ntips*2-1),(ntips+1)))
    {
	k<-which(t$edge[,1]==j) ;
	ch<-t$edge[k,2] ;
	k1<-which(lo2==ch[1]) ;
	k2<-which(lo2==ch[2]) ;
	if (max(k1)>max(k2))
	{
		tmp<-k1 ;
		k1<-k2 ;
		k2<-tmp ;
	}
	lo2[k1]<-j
	lo2[k2]<-j
	#print(lo2)
	s<-c(s,max(k1))
    }
    return(list(s=s,lo=lo))
}

edges_from_lo_s<-function(s,lo)
{
    ntips<-length(lo)
    NEdges<-(ntips-1)*2
    lo2<-lo
    e<-matrix(0,nrow=NEdges,ncol=2)
    for (j in seq((ntips*2-1),(ntips+1)))
    {
    	s_i<-s[(ntips*2)-j]
    	ch1<-lo2[s_i]
    	ch2<-lo2[s_i+1]
    	k1<-which(lo2==ch1)
    	k2<-which(lo2==ch2)
    	k<-2*(j-ntips)
    	print(k)
    	e[k-1,]<-c(j,ch1)
    	e[k,]<-c(j,ch2)
    	lo2[k1]<-j
    	lo2[k2]<-j
    }
    return(e)
}

phylo_ne<-function(tr,h=-1)
{
	# get parent, child, grandchild, sibling 
	fam<-ne_findauntniece(tr)
	
	# constrain by host if necessary
	if (length(h)>1)
	{
		nlmx<-matrix(h[fam],nrow(fam),ncol(fam)) ;
		i<-which(nlmx[3,]==nlmx[4,]) ;
		fam<-fam[,i] ;
	}
	# constrain by fam 
	
	fam<-ne_filterauntniecebylength(tr,fam)
	# pick group
	if(ncol(fam)>0)
	{
		r<-ifelse(ncol(fam)>1,sample(1:ncol(fam),1),1) ;
		aunt<-fam[4,r] ;
		niece<-fam[3,r] ;
		#print(aunt)
		#print(niece)
		tr<-ne_swapnodes(tr,aunt,niece)
	}
	return(tr)
}

phylo_lo_s_ne<-function(tr,lo,s)
{
	# get parent, child, grandchild, sibling 
	fam<-ne_findauntniece(tr)
	
	# constrain by host if necessary
	
	# constrain by fam 
	
	fam<-ne_filterauntniecebylength(tr,fam)
	# pick group
	if(ncol(fam)>0)
	{
		r<-ifelse(ncol(fam)>1,sample(1:ncol(fam),1),1) ;
		aunt<-fam[4,r] ;
		niece<-fam[3,r] ;
		print(aunt)
		print(niece)
		
		
		tr<-ne_swapnodes(tr,aunt,niece)
	}
}

ne_swapnodes_lo_s<-function(tr,lo,s,a,b)
{
	# a is the aunt
	# b is the niece
	# parent of niece
	p_i<-which(tr$edge[,2]==b)
	# grandparent of niece
	gp_i<-which(tr$edge[,2]==tr$edge[p_i,1])
	# sister of niece
	n2_i<-setdiff(which(tr$edge[,1]==tr$edge[p_i,1]),p_i)
	n2<-tr$edge[n2_i,2]
	
	# get aunt's tips
	tip_a<-unlist(Descendants(tr,a,type="tips")) ;
	
	# get niece's tips
	tip_n<-unlist(Descendants(tr,b,type="tips")) ;
	
	# get niece2's tips
	tip_n2<-unlist(Descendants(tr,n2,type="tips")) ;
	
}

ne_swapnodes<-function(tr,a,b)
{
	# a is the aunt,
	# b is the niece
	ia<-which(tr$edge[,2]==a)
	ib<-which(tr$edge[,2]==b)
	ic<-which(tr$edge[,2]==tr$edge[ib,1])
	tr$edge[ia,2]<-b
	tr$edge[ib,2]<-a
	ea<-tr$edge.length[ia]
	eb<-tr$edge.length[ib]
	ec<-tr$edge.length[ic]
	tr$edge.length[ib]<-ea-ec
	tr$edge.length[ia]<-eb+ec
	return(tr)
	
}

ne_filterauntniecebylength<-function(tr,fam)
{
	# only retain valid swaps i.e. parent-aunt has to be longer than parent-sister
	# parent aunt
	e1<-tr$edge.length[match(fam[4,],tr$edge[,2])]
	# parent child
	e2<-tr$edge.length[match(fam[2,],tr$edge[,2])]
	fam<-fam[,e1>e2]
	return(fam)
}

ne_findauntniece<-function(tr)
{
	ntips<-length(tr$tip.label)
	ino<-order(tr$edge[,1]) ;
	inode<-seq((ntips+1),(2*ntips - 1) );
	ch<-matrix(tr$edge[ino,2],nrow=tr$Nnode,ncol=2,byrow=TRUE) ;
	i1<-tr$edge[,1] %in% ch[,1] ;
	i2<-tr$edge[,1] %in% ch[,2] ;
	
	# children  :
	l1_c<-tr$edge[i1,1] ;
	l2_c<-tr$edge[i2,1] ;
	
	# with parents :
	l1_p<-inode[match(l1_c,ch[,1])]; # parents of ch1
	l2_p<-inode[match(l2_c,ch[,2])]; # parents of ch2
	
	# grandchildren :
	l1_gc<-tr$edge[i1,2] ;
	l2_gc<-tr$edge[i2,2] ;
	
	# and siblings :
	l1_s<-ch[l1_p-ntips,2] ;
	l2_s<-ch[l2_p-ntips,1] ;
	
	mx <- cbind(rbind(l1_p,l1_c,l1_gc,l1_s),rbind(l2_p,l2_c,l2_gc,l2_s)) ;
	return(mx)
}