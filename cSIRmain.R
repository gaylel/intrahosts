cSIR<-function(I0=1,nS=100,br=2,dr=1,maxiters=50)
{
	S<-nS-I0 
	I<-I0 
	T<-0 

	n<-1 
	while(I[n]>0 & n<maxiters) {
		r<-c((br/nS)*I[n]*S[n],dr*I[n]) 
		
		# draw time interval
		T[n+1]<-rexp(1,sum(r)) + T[n]
		
		# draw event
		e<- sample(c(1,2),1,prob=r) 
		
		I[n+1]<-I[n] 
		S[n+1]<-S[n] 
		if (e==1)
		{
			I[n+1]<-I[n+1] + 1
			S[n+1]<-S[n+1] - 1
		}
		else
		{
			I[n+1]<-I[n+1] - 1
		}
		
		
		n<-n+1 
		
	}
	return(list(S=S,I=I,T=T)) 
}

cSIR_multi<-function(nHosts=2,I0=1,nS=100,br=2,br2=1,dr=1,maxiters=50)
{
	# between host transmission rates
	B<-matrix(br2,nrow=nHosts,ncol=nHosts) 
	for (i in seq(1,nHosts))
	{
		B[i,i]<-br 
	}
	
	S<-matrix(nS,nrow=1,ncol=nHosts) 
	S[1,1]<-nS-I0 
	
	I<-matrix(0,nrow=1,ncol=nHosts) 
	I[1,1]<-I0 
	
	T<-matrix(0,nrow=1,ncol=4) 
	T[1,1]<-0 

	R<-matrix(0,nrow=nHosts,ncol=nHosts+1) 
	n<-1 
	
	while(any(I[n,]>0) & n<maxiters) {
		
		R[,nHosts+1]<-I[n,]*dr 
		R[1:nHosts,1:nHosts]<-matrix(I[n,],nHosts,nHosts,byrow=TRUE)*matrix(S[n,],nHosts,nHosts)*B/nS 
		#r<-c((br/nS)*I[n]*S[n],dr*I[n]) 
		
		# draw time interval
		T<-rbind(T,rep(0,4)) 
		T[n+1,1]<-rexp(1,sum(R)) + T[n]
		
		# draw event
		e<- sample(1:length(R),1,prob=as.vector(R)) 
		ec<-ceiling(e/nHosts) 
		er<-e-(ec-1)*nHosts 
		
		I<-rbind(I,I[n,])
		#I[n+1,]<-I[n,] 
		
		S<-rbind(S,S[n,]) 
		#S[n+1,]<-S[n,] 
		
		if (ec==(nHosts+1))
		{
			I[n+1,er]<-I[n+1,er] - 1
			T[n+1,4]<- -1 
			T[n+1,2:3]<-er 
		}
		else
		{
		
			I[n+1,er]<-I[n+1,er] + 1
			S[n+1,er]<-S[n+1,er] - 1
			T[n+1,4]<-1 
			T[n+1,2]<-ec 
			T[n+1,3]<-er 
		}
		
		
		
		
		n<-n+1 
		
	}
	return(list(S=S,I=I,T=T)) 
}

cSIR_eval<-function(sir,SN)
{

	p<-1;
	
	Iend<-sir$Iend ;
	NHosts<-length(Iend) ;
	i=1 ;
	while (p==1 & i<=NHosts)
	{
	 	#if (Iend[i]<SN[i] || Iend[i]>(5*SN[i])) p<-0 ;
	 	if (Iend[i]<SN[i] ) p<-0 ;
	 	i<-i+1 ;
	} 
	return(p) ;
}

cSIR_reconstructedtree<-function(tre, T)
{

	# find lineages that are extinct before time T
	nnode<-t$Nnode ;
	ntips<-length(t$tip.label) ;
	d<-node.depth.edgelength(t) ;
	dtips<-d[1:ntips] ;
	dnodes<-d[ntips+1:(ntips+nnode)] ;	
	stips<-which(dtips>=T) ;
	
	# look for edges spanning the cutoff
	e<-cbind(d[t$edge[,1]],d[t$edge[,2]]) ;
	
	# new tips
	sp<-intersect(which(e[,1]<T),which(e[,2]>T)) ;
	tips2<-t$edge[sp,2] ;
	
	
	edges<-t$edge[sp,] ;
	el<-t$edge.length[sp] ;
	
	ti<-NULL ;
	while (nrow(edges)>0)
	{
	# get parents
	pa <- which(tabulate(edges[,1])==2) ;
	
	# most recent coalescence event
	tii<- which.max(d[pa]) ;
	ti<-cbind(ti,d[pa[tii]]) ;
	ch<-which(edges[,1]==pa[tii]) ;
	edges<-edges[-ch,] ;
	edges<-rbind(edges,t$edge[which(t$edge[,2]==pa[tii]),] ) ;
	edges
	}
	
	nodes<-unique(t$edge[sp,1]) ;
	while (length(nodes)>0)
	{
		sp<-match(nodes,t$edge[,2],nomatch=0) ;
		edges<-rbind(edges,t$edge[sp[sp>0],]) ;
		el<-c(el,t$edge.length[sp[sp>0]]) ;
		nodes<-unique(t$edge[sp[sp>0],1]) ;
	}
	
	dup<-duplicated(edges) ;
	edges<-edges[!dup,] ;
	el<-el[!dup] ;
	
	# tidy and relabel tree
	
	nodes<-tips2 ;
	
	for (i in seq(1,length(nodes)))
	{
		j<-nodes[i] ;
		pa_i<-match(j,edges[,2]) ;
		pa<-edges[pa_i,1] ;
		nc<-length(which(edges[,1]==pa)) ;
		while (nc==1)
		{
			pa_i1<-match(pa,edges[,2]) ;
			pa<-edges[pa_i1,1] ;
			el[pa_i]=el[pa_i]+el[pa_i1] ;
			edges[pa_i,1]<-edges[pa_i1,1] ;
			el<-el[-pa_i1] ;
			edges<-edges[-pa_i1,] ;	
			nc<-length(which(edges[,1]==pa)) ;
		}
	
	

	}
}

cSIRtree_reconstruct<-function(sir,N=100,SN,ST)
{
	I<-sir[[1]] ;
	S<-sir[[2]] ;
	T<-sir[[3]] ;
	Iend<-sir[[4]] ;
	nT<-nrow(T) ;
	nHosts<-ncol(I) ;
	unborn=rep(1,nHosts) ;
	
	# graph edges
	g<-NULL ;
	
	# edge lengths
	el<-NULL ;

	# node labels
	nl<-NULL;
	
	NNodes<-sum(SN) ;
	Intnode<-NNodes+NNodes-1 ;
	Nodes<-list() ;
	st=0 ;
	# save node indices and corresponding times
	
	for (i in seq(1,nHosts))
	{
		Nodes[[i]]<-rbind(seq(st+1,SN[i]+st),rep(ST[i],SN[i])) ;
		nl<-c(nl,rep(i,SN[i])) ;
		
		# add extra nodes
		if (Iend[i]>SN[i])
		{
			enodes<-rbind(seq(-1, -(Iend[i]-SN[i])),rep(ST[i],Iend[i]-SN[i])) ;
			Nodes[[i]]<-cbind(Nodes[[i]],enodes) ;
		}
		st=st+SN[i] ;
	}
	
	Tend<-vector();
	for (i in seq(1,nHosts))
	{
		Tend[i]<-max(which(T[,2]==i & T[,3]==i)) ;
	}
	
	
	for (t in seq(nT-1,1,by=-1))
	{
		ha<-T[t,2] ;
		hb<-T[t,3] ;
		
		# skip terminal event
		if (!(ha==hb & t>=Tend[ha]))
		{
			if (T[t,4]==-1)
			{
				# death - but birth in reverse
				mn<-ifelse(ncol(Nodes[[ha]])==0,0,min(Nodes[[ha]])) ;
				if (mn<= -1) 
				{
					new_node<-mn-1 ;
				}
				else
				{
					new_node<- -1 ;
				}
				Nodes[[ha]]<-cbind(Nodes[[ha]],matrix(c(new_node,T[t,1]),2,1)) ;
			}	
			if (T[t,4]==1)
			{
				# birth - but coalescent event in reverse in host ha
				if (ha==hb)
				{
					ch<-sample(Nodes[[ha]][1,],2) ;
				}
				else
				{
					nha<-ncol(Nodes[[ha]]) ;
					nhb<-ncol(Nodes[[hb]]) ;
					ch1<-ifelse(nha>1,sample(Nodes[[ha]][1,],1),Nodes[[ha]][1,1]) ;
					ch2<-ifelse(nhb>1,sample(Nodes[[hb]][1,],1),Nodes[[hb]][1,1]) ;
					ch<-c(ch1,ch2);
				}
				
				# if any extinct nodes are sampled
				#ia<-match(ch[1],Nodes[[ha]][1,]) ;
				#ib<-match(ch[2],Nodes[[hb]][1,]) ;
				ia<-which(Nodes[[ha]][1,]==ch[1]) ;
				ib<-which(Nodes[[hb]][1,]==ch[2]) ;
				
				if (ch[2]<0)
				{
	
					Nodes[[hb]]<-matrix(nrow=2,Nodes[[hb]][,-ib] );
				}
				else
				{
					if (ch[1]<0)
					{
						# move node from Hb into Ha
						
						Nodes[[ha]]<-matrix(Nodes[[ha]][,-ia],nrow=2) ;
						#ib<-match(ch[2],Nodes[[hb]][1,]) ;
						
						ib<-which(Nodes[[hb]][1,]==ch[2] );
						Nodes[[ha]]<-cbind(Nodes[[ha]],Nodes[[hb]][,ib]) ;
						Nodes[[hb]]<-matrix(Nodes[[hb]][,-ib],nrow=2 );
					}
					else{
						# coalescent event
						
						# create pair of edges with new internal node
						g<-rbind(g,c(Intnode,ch[1])) ;
						g<-rbind(g,c(Intnode,ch[2])) ;
						el<-c(el,Nodes[[ha]][2,ia] - T[t,1]) ;
						el<-c(el,Nodes[[hb]][2,ib] - T[t,1]) ;
						nl<-c(nl,ha) ;
						
						Nodes[[ha]]<-matrix(Nodes[[ha]][,-ia],nrow=2) ;
						
						ib<-which(Nodes[[hb]][1,]==ch[2]) ;
						#ib<-match(ch[2],Nodes[[hb]][1,]) ;
						Nodes[[hb]]<-matrix(Nodes[[hb]][,-ib],nrow=2) ;
						
						Nodes[[ha]]<-cbind(Nodes[[ha]],matrix(c(Intnode,T[t,1]),2,1)) ;
						Intnode<-Intnode-1 ;
						
					}
				}
				
				
				
				
				
			}
		
		
		
		}	
	}
	tr<-list() ;
	tr$edge <- g[nrow(g):1,] ; #g[order(g[,1]),] ;
	tr$Nnode<-NNodes-1 ;
	#t2$node.label<-nodelabel2 ;
	tr$edge.length<-rev(el) ;
	tlab<-NULL
	for (i in seq(1, nHosts))
	{
		hname<-paste("H",i,sep="") ;
		tlab<-c(tlab,paste(hname,paste("S",1:SN[i],sep=""),sep=""))
		
		
	}
	tr$node.label=as.character(nl[(NNodes+1):(NNodes+(NNodes-1))]) ;
	tr$node.label<-rev(tr$node.label) ;
	tr$tip.label<-tlab ;
	class(tr)<-"phylo" ;
		
	return(tr) ;
}


cSIR_p_recover<-function(nS,br,dr,s)
{
	p<-dr/(dr + ((br/nS)*s)) ;
	return(p) ;
}


sample_cSIR_C <-function(I0, NS, NHosts, B, dr)
{
	sir1<-.Call("sample_cSIR_R",I0, NS, NHosts, B, dr) ;
	sir<-list(I=sir1[[1]],S=sir1[[2]],T=sir1[[3]]) ;
	sir$T[,2:3]<-sir$T[,2:3]+1 ;
	return(sir) ;
}

sample_cSIR_S_C <-function(I0, NS, NHosts, B, dr, SN, ST)
{
	sir1<-.Call("sample_cSIR_S_R",I0, NS, NHosts, B, dr, ST, SN) ;
	sir<-list(I=sir1[[1]],S=sir1[[2]],T=sir1[[3]],Iend=sir1[[4]]) ;
	sir$T[,2:3]<-sir$T[,2:3]+1 ;
	sir$B<-B ;
	return(sir) ;
}

cSIR<-function(nHosts=2,I0=1,nS=100,br=2,br2=1,dr=1)
{
	
	#set.seed(runif(1)) ;
	# between host transmission rates
	B<-matrix(br2,nrow=nHosts,ncol=nHosts) ;
	for (i in seq(1,nHosts))
	{
		B[i,i]<-br ;
	}
	sir<-sample_cSIR_C(I0=I0, NS=nS, NHosts = nHosts, B=B, dr=dr) ;
	return(sir) ;
}


cSIR_S<-function(nHosts=2,I0=1,nS=100,br=2,br2=1,dr=1,ST=c(1,1),SN=c(30,30))
{
	
	#set.seed(runif(1)) ;
	# between host transmission rates
	B<-matrix(br2,nrow=nHosts,ncol=nHosts) ;
	for (i in seq(1,nHosts))
	{
		B[i,i]<-br ;
	}
	sir<-sample_cSIR_S_C(I0=I0, NS=nS, NHosts = nHosts, B=B, dr=dr, SN=SN, ST=ST) ;
	return(sir) ;
}

cSIR_SB<-function(nHosts=2,I0=1,nS=100,B,dr=1,ST=c(1,1),SN=c(30,30))
{
	
	B<-matrix(0.001*runif(nHosts*nHosts),nrow=nHosts,ncol=nHosts) ;
	for (i in seq(1,nHosts))
	{
		B[i,i]=1 ;
	}
	sir<-sample_cSIR_S_C(I0=I0, NS=nS, NHosts = nHosts, B=B, dr=dr, SN=SN, ST=ST) ;
	return(sir) ;
}

cSIR_Bproposal<-function(oldB)
{
	nHosts<-nrow(oldB) ;
	B<-exp(log(oldB)+rnorm(nHosts*nHosts,0,0.5)) ;
	return(B) ;
}

cSIR_Bijproposal<-function(oldB,Bcon)
{
	B1<-which(Bcon==1);
	B<-oldB ; 
	ind<-ifelse(length(B1)>1,sample(B1,1),1) ;
	B[ind]<-exp(log(oldB[ind])+rnorm(1,0,0.2)) ;
	return(B) ;
}

cSIR_drproposal<-function(olddr)
{
	dr<-exp(log(olddr)+rnorm(1,0,0.5)) ;
	return(dr) ;
}

cSIR_mrproposal<-function(oldmr)
{
	mr<-exp(log(oldmr)+rnorm(1,0,0.1)) ;
	return(mr) ;
}

cSIR_upgma<-function(seqs,SN,ST)
{
	nHosts<-length(SN) ;
	cSN<-c(0,cumsum(SN)) ;
	ho<-rev(order(ST)) ;
	i1<-0 ;
	lo<-NULL ;
	offset<- 0;
	s<-NULL ;
	for (i in ho)
	{
		i2<-cSN[i+1] ;
		i1<-cSN[i]+1 ;
		print(sprintf("%i %i",i1,i2)) ;
		lo2<-1 ;
		if (i2>i1)
		{
		tr<-upgma(dist.hamming(phyDat(seqs[i1:i2]))) ;
		lo2<-tr$edge[tr$edge[,2]<=SN[i],2] ;
		}
		lo<-cbind(lo,rbind(lo2+i1-1,rep(i,SN[i]))) ;
		
	
		# get s from phylo.
		
		if (i2>i1)
		{
			ntips<-SN[i];
		
			for (j in seq((ntips*2-1),(ntips+1)))
			{
				k<-which(tr$edge[,1]==j) ;
				ch<-tr$edge[k,2] ;
						
				k1<-which(lo2==ch[1]) ;
				k2<-which(lo2==ch[2]) ;
				if (max(k1)>max(k2))
				{
					tmp<-k1 ;
					k1<-k2 ;
					k2<-tmp ;
				}
				s<-c(s,max(k1)+offset) ;
				lo2[k1]<-j ;
				lo2[k2]<-j ;	
				
			}
		}
		
		s<-c(s,length(s)+1) ;
		offset<-offset+SN[i] ;
	}
	s<-s[-length(s)] ;
	
	return(list(lo=lo,s=s)) ;
}


cSIR_Bdrprior<-function(B,Bcon,dr,l1,l2,l3) 
{
	p<-0 ;
	nHosts<-nrow(B) ;
	ll1=log(l1);
	ll2=log(l2);
	ll3=log(l3);
	# use exponential distributions
	for (i in seq(1,nHosts))
	{
		for (j in seq(1,nHosts))
		{
			if (Bcon[i,j]==1)
			{
			if (i==j)
			{
				#p=p*(l1*exp(-l1*B[i,j]))
				p=p+ ll1 - l1*B[i,j] ;			
			}
			else
			{
				#p=p*(l2*exp(-l2*B[i,j])) ;
				p=p+ ll2 - l2*B[i,j] ;
			}
			}
		}
	}
	p=p+ ll3 - l3*dr  ;
	#p=p*(l3*exp(-l3*dr))
	return(p) ;
}

cSIR_Bstruct<-function(ST)
{
	# constrain B so that there are less parameters
	
	# birth rates within host > birth / transmission rate between hosts
	NHosts<-length(ST) ;
	Bcon<-matrix(0,nrow=NHosts,ncol=NHosts) ;
	for (i in seq(1,NHosts))
	{
		for (j in seq(1,NHosts))
		{
			if (ST[j]-ST[i] >= 0 & ST[j]-ST[i] <=14)
			{
			Bcon[i,j]=1 ;
			}
		}
	}

	return(Bcon) ;
}

<<<<<<< HEAD
=======
cSIR_ne<-function(tr,SN,lo,s)
{
	nHosts<-length(SN)
	nl<-NULL
	for (i in seq(1,nHosts))
	{
		nl<-c(nl,rep(i,SN[i])) ;
	}
	nl<-c(nl,as.integer(tr$node.label)) ;
	tr<-phylo_ne(tr,nl)
	res<-s_lo_from_phylo(tr)
	h<-match(res$lo,lo[1,])
	res$lo<-rbind(res$lo,lo[2,h])
	return(list(lo=res$lo,s=res$s,tr=tr))
}
>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714


cSIR_narrowexchange<-function(tr,SN, lo, s)
{
	ntips<-length(tr$tip.label)
	nHosts<-length(SN) ;
	nl<-NULL ;
	
	# host labels for each node
	for (i in seq(1,nHosts))
	{
		nl<-c(nl,rep(i,SN[i])) ;
	}
	nl<-c(nl,as.integer(tr$node.label)) ;
	
	
	inode<-seq((ntips+1),(2*ntips - 1) );
	ne=0 ;
	ino<-order(tr$edge[,1]) ;
	
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
	nlmx<-matrix(nl[mx],nrow(mx),ncol(mx)) ;

	i<-which(nlmx[2,]==nlmx[4,]) ;
	mx<-mx[,i] ;
	nlmx<-nlmx[,i] ;
	el_ps<-tr$edge.length[match(mx[4,],tr$edge[,2])] ;
	el_pc<-tr$edge.length[match(mx[2,],tr$edge[,2])] ;
	mx<-mx[,el_ps>el_pc] ;
	
	
	# sample aunt, niece
	r<-ifelse(ncol(mx)>1,sample(1:ncol(mx),1),1) ;
	aunt<-mx[4,r] ;
	niece<-mx[3,r] ;
	niece2<-setdiff(mx[3,which(mx[4,]==aunt)],niece) ;
	
	# get aunt's tips
	tip_a<-unlist(Descendants(tr,aunt,type="tips")) ;
	intn_a<-unlist(Descendants(tr,aunt,type="all")) ;
	intn_a<-setdiff(intn_a,tip_a) ;

	# get niece's tips
	tip_n<-unlist(Descendants(tr,niece,type="tips")) ;
	intn_n<-unlist(Descendants(tr,niece,type="all")) ;
	intn_n<-setdiff(intn_n,tip_n) ;


	# get niece2's tips
	tip_n2<-unlist(Descendants(tr,niece2,type="tips")) ;
	intn_n2<-unlist(Descendants(tr,niece2,type="all")) ;
	intn_n2<-setdiff(intn_n2,tip_n2) ;

	
	# look for tips in lo
	tip_ai<-match(tip_a,lo[1,]) ;
	tip_ni<-match(tip_n,lo[1,]) ;
	tip_n2i<-match(tip_n2,lo[1,]) ;
	lo_i<-seq(1,length(lo[1,])) ;
	
	# look for internal nodes in s
	allnodes<-ntips+ntips ;
	s_ai<-allnodes-c(intn_a,aunt[aunt>ntips]) ;
	s_ni<-allnodes-c(intn_n,niece[niece>ntips]) ;
	s_n2i<-allnodes -c(intn_n2,niece2[niece2>ntips]) ;
	s_p<-allnodes-mx[1,r] ;
	s_s<-allnodes-mx[2,r] ;
	
	
	alen<-length(tip_ai) ; 
	nlen<-length(tip_ni) ;
	n2len <-length(tip_n2i) ;
	if (max(tip_ni) > max(tip_ai))
	{ 
		
		if (max(tip_n2i) > max(tip_ni))
		{
			# a,n,n2
			lo_i[tip_ai]<- lo_i[tip_ai] + nlen ;
			lo_i[tip_ni]<- lo_i[tip_ni] - alen ;
			s[s_ai]<-s[s_ai] + nlen ;
			s[s_ni]<-s[s_ni] - alen ;
			s[s_p]<-s[s_p] - alen + nlen ;

		}
		else
		{
			# a n2 n
			lo_i[tip_ai]<- lo_i[tip_ai] + nlen + n2len ;
			lo_i[tip_ni]<- lo_i[tip_ni] - n2len - alen ;
			lo_i[tip_n2i]<- lo_i[tip_n2i] + nlen - alen ;
			s[s_ai]<-s[s_ai] + nlen + n2len;
			s[s_ni]<-s[s_ni] - n2len - alen;
			s[s_n2i]<-s[s_n2i] + nlen - alen;
			s[s_p]<-s[s_p] - alen + nlen ;
			s[s_s]<-s[s_s] - alen + nlen ;
		}
	}
	else
	{
		if (max(tip_n2i) > max(tip_ni))
		{
			# n n2 a
			lo_i[tip_ai]<- lo_i[tip_ai] - nlen - n2len ;
			lo_i[tip_ni]<- lo_i[tip_ni] + n2len + alen ;
			lo_i[tip_n2i]<- lo_i[tip_n2i] - nlen + alen ;
			s[s_ai]<-s[s_ai] - nlen - n2len;
			s[s_ni]<-s[s_ni] + n2len + alen;
			s[s_n2i]<-s[s_n2i] - nlen + alen;
			s[s_p]<-s[s_p] + alen - nlen ;
			s[s_s]<-s[s_s] + alen - nlen ;


		}
		else
		{
			# n2 n a
			lo_i[tip_ai]<- lo_i[tip_ai] - nlen ;
			lo_i[tip_ni]<- lo_i[tip_ni] + alen ;
			s[s_ai]<-s[s_ai] - nlen ;
			s[s_ni]<-s[s_ni] + alen ;
			s[s_p]<-s[s_p] + alen - nlen ;
		}
	}

	lo[,lo_i] <- lo
	# delete s-n
	e_n<-which(tr$edge[,2]==niece) ;
	e_s<-which(tr$edge[,2]==mx[2,r]) ;
	e_a<-which(tr$edge[,2]==aunt) ;
	
	
	# add p-n
	tr$edge[e_n,1]<-mx[1,r] ;
	# change edge length to p-s + s-n
	tr$edge.length[e_n] <- tr$edge.length[e_n] + tr$edge.length[e_s] ;
	
	# delete p-a
	# add s-a
	tr$edge[e_a,1]<-mx[2,r] ;
	# change edge length to p-a - p-s
	tr$edge.length[e_a] <- tr$edge.length[e_a] - tr$edge.length[e_s] ;
	return(list(tr=tr,lo=lo,s=s))
		
}	
	
	
	
	
	
	
	
cSIR_SB_metrop<-function(nHosts=2,I0=1,nS=100,dr=1,dat,x, N=1000)
{
	
	# initialise the parameters
	dr=0.1 ;
	mr=1e-3 ;
	B<-matrix(0.001*runif(nHosts*nHosts),nrow=nHosts,ncol=nHosts) ;
	SN<-dat$SN ;
	ST<-dat$ST ;
	for (i in seq(1,nHosts))
	{
		B[i,i]=1 ;
	}
	Bcon<-cSIR_Bstruct(ST=ST) ;
	
	B<-B*Bcon ;
	nB<-length(Bcon>0) ;
	l1=1 ;
	l2=1/0.1;
	l3=1/0.1;
	Nacc=0 ;
	Nrej=0 ;
	iend<-NULL ;
	Bchain<-list() ;
	drchain<-list() ;
	sirchain<-list() ;
	trchain<-list() ;
	llchain<-list() ;
	mrchain<-list() ;
	pold<-cSIR_Bdrprior(B,Bcon,dr,l1,l2,l3)  ;
	M=3;
	
	# initialise branching order and leaf order here.
	#s<-(sum(SN)-1):1 ;
	#lo<- NULL ;
	#for (i in seq(1,nHosts))
	#{
	#	lo<-c(lo,rep(i,SN[i])) ;
	#	
	#}
	#lo<-rbind(seq(1,sum(SN)),lo) ;
	
	ini<-cSIR_upgma(x,SN,ST)
	lo<-ini$lo ;
	s<-ini$s ;
	print(lo)
	print(s)
	x<-phyDat(x) ;
	llcur<- log(0) ;
	mvec<-c(rep(1,nB),2,3) ;
	while (Nacc<N)
	{
		# make new proposal
		newB<-B ;
		newmr<-mr ;
		newlo<-lo ;
		news<-s ;
		pacc=0 ;
		newdr<-dr ;
		
		m<-sample(seq(1,M),1) ;
		#m<-sample(mvec,1) ;
		#m<-1 ;
		if (Nacc==0)
		{
			m<-1 ;
		}
		else
		{
			newtr<-tr ;
		}
		switch(as.character(m),
		"1"={ 
			NBacc=0 ;
			NBi =0;
			NBmax=500 ;
			# sample from B
			newB<-cSIR_Bproposal(B) ;
			newB<-newB*Bcon ;
			
			#newdr<-0.1 ;
			ct<-0 ;
			newll_min = log(0) ;
			#while (ct<1)
			while (NBacc==0)
			{
				sir<-sample_cSIR_S_C(I0=I0, NS=nS, NHosts = nHosts, B=newB, dr=newdr, SN=SN, ST=ST) ;
				
				if (cSIR_eval(sir,SN)>0)
				{
					NBacc<-1 ;
					pnew<-cSIR_Bdrprior(newB,Bcon,newdr,l1,l2,l3) ;
				
					# sample phylogeny and calculate likelihood
					tre<-cSIR_phylosample(sir=sir,nS=nS,lo=lo,s=s,dat)
					ll<-pml(tre$tr,x, rate=newmr, model="GTR") ; 
				
					if (ll$logLik > newll_min)
					{
						newlo<-tre$lo ;
						news<-tre$s ;
						newtr<-tre$tr ;
					
						newll<-ll$logLik ;
						pold<-cSIR_Bdrprior(B,Bcon,dr,l1,l2,l3)  ;
						#print(ll$logLik-llcur + log(pnew) - log(pold)) ;
						#pacc=min(ll$logLik-llcur + log(pnew) - log(pold),0) ;
						pacc=min(ll$logLik-llcur + pnew - pold,0) ;
						#pacc=(min(exp(ll$logLik-llcur),1)) ;
						newll_min<- ll$logLik ;	
					}	
					
					
					#print(sprintf("%i   %8.4f",ct, newll)) ;
					
				}
				ct<-ct+1 ;
				NBi<-NBi+1 ;
			}
			#print(NBacc)
		},
		"4"={
				newdr<-cSIR_drproposal(dr) ;
				NBacc<-0 ;
				ct=0 ;
				while (NBacc==0 && ct<10)
				{
					sir<-sample_cSIR_S_C(I0=I0, NS=nS, NHosts = nHosts, B=newB, dr=newdr, SN=SN, ST=ST) ;
					if (cSIR_eval(sir,SN)>0)
					{
						NBacc<-1 ;
						# sample phylogeny and calculate likelihood
						tre<-cSIR_phylosample(sir=sir,nS=nS,lo=lo,s=s,dat)
						ll<-pml(tre$tr,x, rate=newmr, model="GTR") ; 
				
						newlo<-tre$lo ;
						news<-tre$s ;
						newtr<-tre$tr ;
					
						newll<-ll$logLik ;
						pacc=min(ll$logLik-llcur,0) ;
						ct=ct+1 ;	
					}	
					
					
				}
		},
		"3"={
			# sample mutation rate
			newmr<-cSIR_mrproposal(mr) ;
			ll<-pml(newtr,x, rate=newmr, model="GTR") ; 
			newll<-ll$logLik ;
			pacc=min(ll$logLik-llcur,0) ;

		},
		"2"={
			# narrow exchange
			res<-cSIR_narrowexchange(newtr,SN, lo, s)
			newlo<-res$lo ;
			news<-res$s ;
			newtr<-res$tr ;
			ll<-pml(newtr,x, rate=newmr, model="GTR") ; 
			newll<-ll$logLik ;
			pacc=min(ll$logLik-llcur,0) ;	
		}
		)
		#print(sprintf("%8.4f",pacc)) ;
		if (log(runif(1))<=pacc)
		{
			Nacc<-Nacc+1 ;
			lo<-newlo ;
			s<- news ;
			dr<-newdr ;
			B<-newB ;
			mr<-newmr ;
			tr<-newtr ;
			trchain[[Nacc]]<-newtr ;
			llchain[[Nacc]]<-newll ;
			
			llcur<-newll ;
			drchain[[Nacc]]<-newdr ;
			Bchain[[Nacc]]<-newB ;
			mrchain[[Nacc]]<-newmr ;
			#sirchain[[Nacc]]<-sir ;
			dr<-newdr ;
			B<-newB ;
			print(sprintf("%i loglikelihood   %8.4f   Mutation rate   %8.4f  Death rate %8.4f Move %i",Nacc,newll,newmr,newdr, m))
			if ((Nacc %% 10)==0)
			{
					plot.phylo(newtr, cex = 0.2) ;
			}
		}
		else
			{Nrej<-Nrej+1}
		
	}
	rList<-list(Nacc=Nacc,Nrej=Nrej,B=Bchain,dr=drchain,mr=mrchain,tr=trchain,ll=llchain) ;
	
}



cSIR_SB_metrop2<-function(nHosts=2,I0=1,nS=100,dr=1,dat,x, N=1000)
<<<<<<< HEAD
{
	
	# initialise the parameters
	dr=0.1 ;
	mr=1e-3 ;
	B<-matrix(0.001*runif(nHosts*nHosts),nrow=nHosts,ncol=nHosts) ;
	SN<-dat$SN ;
	ST<-dat$ST ;
	for (i in seq(1,nHosts))
	{
		B[i,i]=1 ;
	}
	Bcon<-cSIR_Bstruct(ST=ST) ;
	
	B<-B*Bcon ;
	nB<-length(Bcon>0) ;
	l1=1 ;
	l2=1/0.1;
	l3=1/0.1;
	
	iend<-NULL ;
	Bchain<-list() ;
	drchain<-list() ;
	sirchain<-list() ;
	trchain<-list() ;
	llchain<-list() ;
	mrchain<-list() ;
	dvec<-NULL
	drvec<-NULL
	pold<-cSIR_Bdrprior(B,Bcon,dr,l1,l2,l3)  ;
	M=4
	;
	
	# initialise 
	ini<-cSIR_upgma(x,SN,ST)
	lo<-ini$lo ;
	s<-ini$s ;
	x<-phyDat(x) ;
	llcur<- log(0) ;
	params<-list(B=B, dr=dr, mr= mr, ll=llcur, lo=lo, s=s, Bcon=Bcon)
	
	
	#params$B<-matrix(0,nrow=NHosts,ncol=NHosts)
	#for (i in seq(1,NHosts))
	#{
	#	params$B[i,i] <- 1
		
	#}
	#params$B[1,2]<-0.4
	#params$B[2,3]<-0.2
	#params$dr<-0.1
	#params$mr<-1.83e-4
	
	
	hyperparams<-list(l1=l1,l2=l2,l3=l3) 
	dat.params<-list(I0=I0, NS=nS, NHosts=NHosts,x=x)
	# initial tree ;
	
	tre<-cSIR_drawTre(params, dat.params, dat, Ntries=1000) 
	tr<-tre$tr
	params$tr<-tr
	params$T<-tre$T
	mvec<-c(rep(1,nB),2,3) ;
	
	Nacc=1;
	Nrej=0 ;
	while (Nacc<N)
	{
		# make new proposal
		m<-sample(seq(1,M),1) ;
		#m<-sample(mvec,1) ;
		#m<-4;
		m<-sample(c(2,4),1)
		switch(as.character(m),
		"1"={ 
				params<-cSIR_Bupdate(params=params, dat.params=dat.params, dat=dat, hyperparams=hyperparams)
			},
		"5"={
				params<-cSIR_drupdate(params, dat.params, dat)
			
		},
		"3"={
				params<-cSIR_mrupdate(params,dat.params,dat) 
			},
		"2"={
			# narrow exchange
			res<-cSIR_narrowexchange(params$tr,SN, params$lo, params$s)
			lo<-res$lo ;
			s<-res$s ;
			tr<-res$tr ;
			ll<-pml(tr,x, rate=params$mr, model="GTR") ; 
			newll<-ll$logLik ;
			pacc=min(ll$logLik-params$ll,0) ;	
			if (log(runif(1))<=pacc)
			{
				params$lo<-lo
				params$s<-s
				params$tr<-tr
				params$ll<-ll$logLik
			}
		},
		"4"={
				params<-cSIR_Tupdate(params, dat.params, dat)
		}
		)
		#print(params)
		#print(sprintf("%8.4f",pacc)) ;
			trchain[[Nacc]]<-params$tr ;
			llchain[[Nacc]]<-params$ll ;
			drchain[[Nacc]]<-params$dr ;
			Bchain[[Nacc]]<-params$B ;
			mrchain[[Nacc]]<-params$mr ;
			#sirchain[[Nacc]]<-sir ;
			
			print(sprintf("%i loglikelihood   %8.4f   Mutation rate   %8.4f  Death rate %8.4f Move %i",Nacc,params$ll,params$mr,params$dr, m))
			#print(params$T)
			if ((Nacc %% 10)==0)
			{
					#plot.phylo(params$tr, cex = 0.2) ;
			}
			
			dvec<-c(params$dr,dvec)
			#print(params)
			Nacc<-Nacc+1 ;
		
	}
	hist(dvec,50)
	rList<-list(Nacc=Nacc,Nrej=Nrej,B=Bchain,dr=drchain,mr=mrchain,tr=trchain,ll=llchain) ;
	return(rList)
}

cSIR_drawTre<-function(params, dat.params, dat, Ntries=1000)
{
	B<-params$B 
	dr<-params$dr
	lo<-params$lo
	s<-params$s
	I0<-dat.params$I0
	NS<-dat.params$NS
	NHosts<-dat.params$NHosts
	
	
	Nacc<-0
	
	T<-params$T
	tre<-list(tr=params$tr,T=T,s=s,lo=lo)
	nt<-1
	while (Nacc==0 & nt<Ntries)
	{
		sir<-sample_cSIR_S_C(I0=I0, NS=NS, NHosts = NHosts, B=B, dr=dr, SN=dat$SN, ST=dat$ST) 
		#print(sir)
		if (cSIR_eval(sir,dat$SN)>0)
		{
			Nacc<-1 ;
			
			# sample tree (branching times)
			tre<-cSIR_phylosample(sir=sir,nS=NS,lo=lo,s=s,dat)
			tre$T<-diff(tre$T[,1])
			
		}
		nt<-nt+1
	}
	print(Nacc)
	tre$Nacc<-Nacc
	return(tre)
}

cSIR_Tupdate<-function(params, dat.params, dat)
{
	res<-cSIR_drawTre(params=params,dat.params=dat.params,dat=dat)
	T<-res$T
	tr<-res$tr
	lo<-res$lo
	s<-res$s
	llcur<-params$ll
	ll<-pml(tr,dat.params$x, rate=params$mr, model="GTR")  
	pacc=min(ll$logLik-llcur,0)			
	if (log(runif(1))<=pacc)
	{
		params$lo<-lo
		params$s<-s
		params$tr<-tr
		params$ll<-ll$logLik
		params$T<-T
	}
	return(params)
}

cSIR_Bupdate<-function(params, dat.params, dat, hyperparams, th=2)
{
	newB<-cSIR_Bproposal(params$B) 
	newB<-newB*params$Bcon
	oldB<-params$B
	params$B<-newB
	res<-cSIR_drawTre(params=params,dat.params=dat.params,dat=dat)
	params$B<-oldB
	T<-res$T ;
	d<-sum((T-params$T)^2)/(length(T)-1) 
	print(d)
	
	if (d<th & res$Nacc==1)
	{
		pold<-cSIR_Bdrprior(oldB,params$Bcon,params$dr,hyperparams$l1,hyperparams$l2,hyperparams$l3)  ;
		pnew<-cSIR_Bdrprior(newB,params$Bcon,params$dr,hyperparams$l1,hyperparams$l2,hyperparams$l3)  ;
		h<-min(1,exp(pnew - pold)) ;
		print(h)
		
		if (runif(1)<=h)
		{ 
			params$B<-newB
		}
	}
	
	return(params)			
}

cSIR_drupdate<-function(params,dat.params,dat, hyperparams, th=0.1)
{
	newdr<-cSIR_drproposal(params$dr) 
	olddr<-params$dr
	params$dr<-newdr
	res<-cSIR_drawTre(params,dat.params,dat)
	T<-res$T ;
	d<-sum((T-params$T)^2)/(length(T)-1) 
	print(d)
	
	if (d>th | res$Nacc==0)
	{
		params$dr<-olddr
	}
	
=======
{
	source("phylo.R")
	# initialise the parameters
	dr=0.1 ;
	mr=1e-3 ;
	B<-matrix(0.001*runif(nHosts*nHosts),nrow=nHosts,ncol=nHosts) ;
	SN<-dat$SN ;
	ST<-dat$ST ;
	for (i in seq(1,nHosts))
	{
		B[i,i]=1 ;
	}
	Bcon<-cSIR_Bstruct(ST=ST) ;
	
	B<-B*Bcon ;
	nB<-length(Bcon>0) ;
	l1=1 ;
	l2=1/0.1;
	l3=1/0.1;
	
	iend<-NULL ;
	Bchain<-list() ;
	drchain<-list() ;
	sirchain<-list() ;
	trchain<-list() ;
	llchain<-list() ;
	mrchain<-list() ;
	dvec<-NULL
	drvec<-NULL
	pold<-cSIR_Bdrprior(B,Bcon,dr,l1,l2,l3)  ;
	M=4
	
	
	# initialise 
	ini<-cSIR_upgma(x,SN,ST)
	lo<-ini$lo ;
	s<-ini$s ;
	x<-phyDat(x) ;
	llcur<- log(0) ;
	params<-list(B=B, dr=dr, mr= mr, ll=llcur, lo=lo, s=s, Bcon=Bcon)
	print(params$tr$edge)
			print(lo)
			print(s)
	
	#params$B<-matrix(0,nrow=NHosts,ncol=NHosts)
	#for (i in seq(1,NHosts))
	#{
	#	params$B[i,i] <- 1
		
	#}
	#params$B[1,2]<-0.4
	#params$B[2,3]<-0.2
	params$dr<-0.1
	params$mr<-mr
	
	
	hyperparams<-list(l1=l1,l2=l2,l3=l3) 
	dat.params<-list(I0=I0, NS=nS, NHosts=NHosts,x=x)
	# initial tree ;
	
	tre<-cSIR_drawTre(params, dat.params, dat, Ntries=1000) 
	tr<-tre$tr
	params$tr<-tr
	params$T<-tre$T
	print(params$tr$edge)
			print(lo)
			print(s)
	mvec<-c(rep(1,nB),2,3) ;
	
	Nacc=1;
	Nrej=0 ;
	while (Nacc<N)
	{
		# make new proposal
		m<-sample(seq(1,M),1) ;
		#m<-sample(mvec,1) ;
		#m<-4;
		#m<-sample(c(2,4),1)
		m<-1
		switch(as.character(m),
		"1"={ 
				params<-cSIR_Bupdate(params=params, dat.params=dat.params, dat=dat, hyperparams=hyperparams)
			},
		"5"={
				params<-cSIR_drupdate(params, dat.params, dat)
			
		},
		"3"={
				params<-cSIR_mrupdate(params,dat.params,dat) 
			},
		"2"={
			# narrow exchange
			#res<-cSIR_narrowexchange(params$tr,SN, params$lo, params$s)
			res<-cSIR_ne(params$tr,SN, params$lo, params$s)
			
			lo<-res$lo ;
			s<-res$s ;
			tr<-res$tr ;
			ll<-pml(tr,x, rate=params$mr, model="GTR") ; 
			newll<-ll$logLik ;
			pacc=min(ll$logLik-params$ll,0) ;	
			if (log(runif(1))<=pacc)
			{
				params$lo<-lo
				params$s<-s
				params$tr<-tr
				params$ll<-ll$logLik
			}
		},
		"4"={
				params<-cSIR_Tupdate(params, dat.params, dat)
		}
		)
		#print(params)
		#print(sprintf("%8.4f",pacc)) ;
			trchain[[Nacc]]<-params$tr ;
			llchain[[Nacc]]<-params$ll ;
			drchain[[Nacc]]<-params$dr ;
			Bchain[[Nacc]]<-params$B ;
			mrchain[[Nacc]]<-params$mr ;
			#sirchain[[Nacc]]<-sir ;
			
			print(sprintf("%i loglikelihood   %8.4f   Mutation rate   %8.4f  Death rate %8.4f Move %i",Nacc,params$ll,params$mr,params$dr, m))
			#print(params$T)
			#if ((Nacc %% 10)==0)
			#{
			#		plot.phylo(params$tr, cex = 0.8) ;
			#}
			
			dvec<-c(params$dr,dvec)
			#print(params)
			Nacc<-Nacc+1 ;
		    			#print(params$tr$edge)
			#print(lo)
			#print(s)
	}
	hist(dvec,50)
	rList<-list(Nacc=Nacc,Nrej=Nrej,B=Bchain,dr=drchain,mr=mrchain,tr=trchain,ll=llchain) ;
	return(rList)
}

cSIR_drawTre<-function(params, dat.params, dat, Ntries=1000)
{
	B<-params$B 
	dr<-params$dr
	lo<-params$lo
	s<-params$s
	I0<-dat.params$I0
	NS<-dat.params$NS
	NHosts<-dat.params$NHosts
	
	
	Nacc<-0
	
	T<-params$T
	tre<-list(tr=params$tr,T=T,s=s,lo=lo)
	nt<-1
	while (Nacc==0 & nt<Ntries)
	{
		sir<-sample_cSIR_S_C(I0=I0, NS=NS, NHosts = NHosts, B=B, dr=dr, SN=dat$SN, ST=dat$ST) 
		#print(sir)
		if (cSIR_eval(sir,dat$SN)>0)
		{
			Nacc<-1 ;
			
			# sample tree (branching times)
			tre<-cSIR_phylosample(sir=sir,nS=NS,lo=lo,s=s,dat)
			tre$T<-diff(tre$T[,1])
			
		}
		nt<-nt+1
	}
	print(Nacc)
	tre$Nacc<-Nacc
	return(tre)
}

cSIR_Tupdate<-function(params, dat.params, dat)
{
	res<-cSIR_drawTre(params=params,dat.params=dat.params,dat=dat)
	T<-res$T
	tr<-res$tr
	lo<-res$lo
	s<-res$s
	llcur<-params$ll
	ll<-pml(tr,dat.params$x, rate=params$mr, model="GTR")  
	pacc=min(ll$logLik-llcur,0)			
	if (log(runif(1))<=pacc)
	{
		params$lo<-lo
		params$s<-s
		params$tr<-tr
		params$ll<-ll$logLik
		params$T<-T
	}
	return(params)
}

cSIR_Bupdate<-function(params, dat.params, dat, hyperparams, th=2)
{
	newB<-cSIR_Bproposal(params$B) 
	newB<-newB*params$Bcon
	oldB<-params$B
	params$B<-newB
	res<-cSIR_drawTre(params=params,dat.params=dat.params,dat=dat)
	params$B<-oldB
	T<-res$T ;
	d<-sum((T-params$T)^2)/(length(T)-1) 
	print(d)
	
	if (d<th & res$Nacc==1)
	{
		pold<-cSIR_Bdrprior(oldB,params$Bcon,params$dr,hyperparams$l1,hyperparams$l2,hyperparams$l3)  ;
		pnew<-cSIR_Bdrprior(newB,params$Bcon,params$dr,hyperparams$l1,hyperparams$l2,hyperparams$l3)  ;
		h<-min(1,exp(pnew - pold)) ;
		print(h)
		
		if (runif(1)<=h)
		{ 
			params$B<-newB
		}
	}
	
	return(params)			
}

cSIR_drupdate<-function(params,dat.params,dat, hyperparams, th=0.1)
{
	newdr<-cSIR_drproposal(params$dr) 
	olddr<-params$dr
	params$dr<-newdr
	res<-cSIR_drawTre(params,dat.params,dat)
	T<-res$T ;
	d<-sum((T-params$T)^2)/(length(T)-1) 
	print(d)
	
	if (d>th | res$Nacc==0)
	{
		params$dr<-olddr
	}
	
>>>>>>> 9efa3e2269dea9d8828b5ae1f850e3feb146b714
	return(params)	
}

cSIR_mrupdate<-function(params,dat.params,dat)
{
	newmr<-cSIR_mrproposal(params$mr) 
	oldmr<-params$mr
	params$mr<-newmr
	ll<-pml(params$tr, dat.params$x, rate=params$mr, model="GTR")  
	params$mr<-oldmr
	pacc=min(ll$logLik-params$ll,0)			
	if (log(runif(1))<=pacc)
	{
		params$mr<-newmr
		params$ll<-ll$logLik
	}
	
	return(params)	
}




cSIR_phylosample<-function(sir,nS,lo,s,dat)
{
	#tr<-cSIRtree_reconstruct(sir,N=nS,SN=dat$SN,ST=dat$ST) ;
	NHosts<-length(dat$SN) ;
	#xprint(sir$T)
	tr_list<-.Call("tree_reconstruct",sir, NHosts, dat) ;
	#print(tr_list)
	tr<-tr_list[[1]]
	if (length(which(tr$edge.length < 0)) >0) 
	{
		print("negative edges!") ;
	}

	T=tr_list[[2]]
	#print(T)
	#Tp<-cSIR_treeparams(tr,SN=dat$SN,ST=dat$ST) ;
	#T=Tp$T
	
	# get new branching order based on current and leaf order
	bt<-cSIR_stimes(T=T,lo=lo,s=s) ;
	
	
	
	tr2<-cSIR_phylofroms(s=bt$s,lo=bt$lo, T=T, dat=dat) ;
	#print(tr2) ;
	#print(tr2$edge) 
	#print(tr2$edge.length)
	return(list(tr=tr2,s=bt$s, lo=bt$lo, T=T)) ;
	
}

cSIR_phylofroms<-function(s,lo,T,dat)
{
		Nnode<-nrow(T) ;
		NHosts<-length(dat$SN) ;
		ntips<-Nnode+1 ;
		Inode<-2*ntips -1 ;
		lo2<-lo[1,] ;
		# vector of node times
		t<-NULL ;
		for (i in seq(1,NHosts))
		{
			t<-c(t,rep(dat$ST[i],dat$SN[i]));
			
		}
		t<-c(t,(T[,1])) ;
		
		edge<-NULL ;
		el <- NULL ;
		for (i in seq(Nnode,1))
		{
			s1<-s[Nnode+1-i] ;
			s2<-s1+1 ;
			edge<-rbind(c(Inode,lo2[s2]),edge) ;
			edge<-rbind(c(Inode,lo2[s1]),edge) ;	
			j<-which(lo2==lo2[s2]);
			k<-which(lo2==lo2[s1]);
			el<-c(t[lo2[s1]]-t[Inode], t[lo2[s2]]-t[Inode],el)
			lo2[j]<-Inode ;
			lo2[k]<-Inode ;
			Inode<-Inode-1 ;
			#print(edge)
			
		}
		tr<-list() ;
		tr$Nnode<-Nnode ;
		tr$edge <- edge ;
		tr$edge.length <- el ;
		tr$node.label <- as.character(T[,2]) ;
		tr$tip.label <- NULL ;
		for (i in seq(1,NHosts))
		{
			tr$tip.label<-c(tr$tip.label,dat$info[[i]]$key) ;
		}
		class(tr)<-"phylo" ;
		return(tr) ;
}


cSIR_plottraj<-function(sir)
{
	# plot infected trajectories
	I<-sir$I ;
	T<-sir$T ;
	nHosts<-ncol(I) ;
	mI<-max(I)+10 ;
	plot(T[,1],I[,1],ylim=c(0,mI),type="s") ;
	for (i in seq(2,nHosts))
	{
		lines(T[,1],I[,i],type="s") ;
	}

}

cSIR_treeparams<-function(tre,SN,ST)
{
	# read in tree and return topology and branching times
	ntips<-length(tre$tip.label) ;
	T<-matrix(0,nrow=ntips-1,ncol=4) ;
	# current times and node labels
	nHosts<-length(SN) ;
	tcur<-NULL ;
	nl<-NULL ;
	for (i in seq(1,nHosts))
	{
		tcur<-c(tcur,rep(ST[i],SN[i])) ;
		nl<-c(nl,rep(i,SN[i])) ;
	}
	nl<-c(nl,as.integer(tre$node.label)) ;
	tcur<-c(tcur,rep(0,ntips-1)) ;
	
	# internal nodes from most recent to oldest
	inode<-seq((ntips+(ntips-1)),(ntips+1),by=-1) ;
	
	
	for (i in inode)
	{
		
		inode_i<-which(tre$edge[,1]==i) ;
		
		# find tips that i connects
		ch<-tre$edge[inode_i,2] ;
		
		T[i-ntips,1]<-tcur[ch[1]]-tre$edge.length[inode_i[1]] ;
		tcur[i]<-T[i-ntips,1] ;
		T[i-ntips,3]<-nl[ch[1]] ;
		T[i-ntips,4]<-nl[ch[2]] ;
		T[i-ntips,2]<-nl[i] ;
	}
	
	# process T
	
	j<-union(which(T[,2]!=T[,3]),which(T[,2]!=T[,4])) ;
	migNodes<-rev(tre$edge[,1])[j*2] ;
	
	
	# leaf ordering
	#tre2<-reorder(tre) ;
	#lo<-tre2$edge[tre2$edge[,2]<=ntips,2] ;
	#lo<-rbind(lo,nl[lo]);
	
	#lo2<-lo[1,] ;
	#s<-rep(0,length(lo)-1) ;
	# branching time ordering
	#for (i in inode)
	#{
	#	inode_i<-which(tre2$edge[,1]==i) ;	
	#	
	#	# find tips that i connects
	#	ch<-tre2$edge[inode_i,2] ;
	#	
	#	lo_1<- which(lo2==ch[1]) ;
	#	lo_2<- which(lo2==ch[2]) ;
	#	mlo1<-max(lo_1) ;
	#	mlo2<-max(lo_2) ;
	#	if (mlo1 > mlo2)
	#	{
	#		s[mlo2]<-i ;
	#	}
	#	else
	#	{
	#		s[mlo1]<-i ;
	#	}
	#	lo2[lo_1]<-i ;
	#	lo2[lo_2]<-i ;		
	#}
	
	
	#return(list(T=T,mn=migNodes,s=s,lo=lo)) ;
	return(list(T=T,mn=migNodes)) ;
}

cSIR_updatebh<-function(lbc,rbc,bh,h,rt)
{
	
		# lbc - cluster of leaves on left branch
		# rbc - cluster of leaves on right branch
		# h - mrcas of each leaf (ordered)
		# bh - matrix of coalescent event ordering, and which hosts they occur in	
	
	

		# update bh from h 
		
		lbcm<-max(lbc) ;
		rbcm<-max(rbc) ;
		
		# check if l-r ordering is right
		if (rbcm>lbcm)
		{
			lbcm=min(lbc) ;
		}
		else
		{
			rbcm=lbcm ;
			lbcm=min(rbc) ;
		}
			
			
		if (rbcm<length(h))
		{
			bh_i<-which(bh[3,]==rbcm) ;
			bh[1,bh_i]<-rt ;
		}
		
		if (lbcm>1)
		{
			bh_i<-which(bh[3,]==lbcm-1) ;
			bh[2,bh_i]<-rt ;
		}
		
		return(bh) ;
}


cSIR_updatesnew<-function(lo,snew)
{
	snew<-match(snew,lo[1,]) ;
}

cSIR_stimes<-function(T,lo,s)
{
	N<-nrow(T) ;
	ntips<-ncol(lo) ;
	
	# keep track of which host each internal node is in.
	Inodes<-(ntips+1):(2*ntips-1) ;
	Inodes<-rbind(Inodes,rep(0,length(Inodes))) ;
	Inode<-2*ntips-1 ;
	lbh<-lo[2,s] ;
	rbh<-lo[2,s+1] ;
	
	h<-lo[1,] ;
	bh<-rbind(lbh,rbh,s) ;
	lo_s<-lo[,order(lo[1,])] ;
	
	snew<-NULL ;
	for (i in seq(N,1))
	{

		ha<-T[i,3] ;
		hb<-T[i,4] ;
		j<-which((bh[1,]==ha & bh[2,]==hb) |  (bh[1,]==hb & bh[2,]==ha)) ;
		if (length(j)>0)
		{
			e<-min(j) ;
			snew<-c(snew,bh[3,e]) ;
			lb<-h[bh[3,e]] ;
			rb<-h[bh[3,e]+1] ;
			lbc<-which(h==lb) ;
			rbc<-which(h==rb) ;
			h[lbc]<-Inode ;
			h[rbc]<-Inode ;
			Inode<-Inode-1 ;
			Inodes[2,i]<-T[i,2] ;
			
			bh<-cSIR_updatebh(lbc,rbc,bh,h,T[i,2]) ;
			bh<-matrix(bh[,-e],nrow=3) ;
		}

		else
		{
			# have to change the leaf ordering
			# check which nodes are available
			
			hh<-cbind(lo_s,Inodes) ;
			hu<-unique(h) ;
			ho<-hh[2,hu] ;
			
			lb<-which(ho==ha) ;
			rb<-which(ho==hb) ;
			
			if (ha!=hb)
			{
				lbi<-ifelse(length(lb)>1,sample(lb),lb) ;
				rbi<-ifelse(length(rb)>1,sample(rb),rb) ;
			}
			else
			{
				lbi<-ifelse(length(lb)>1,sample(lb),lb) ;
				lbt<-setdiff(lb,lbi) ;
				rbi<-ifelse(length(lbt)>1,sample(lbt),lbt) ;
				
			}
			lbv<-hh[1,hu[lbi]] ;
			rbv<-hh[1,hu[rbi]] ;
			
			# shift rb to be next to lb
			
			hri<- which(h==rbv) ;
			hli<- which(h==lbv) ;
		
			# test which is on the lhs
			if (min(hli)>min(hri))
			{
				tmp<-hli ;
				hli<-hri ;
				hri<-tmp ;
			}
			

			seg1<-(max(hli)+1):(min(hri)-1) ;
			lo_o1<-lo[1,hri] ;
			lo_o2<-lo[1,seg1];
			# update leaf ordering
			if (max(hri)<ntips)
			{
				seg2<-(max(hri)+1):ntips;
				h<-c(h[-c(seg1,seg2)],h[seg1],h[seg2]) ;
				lo<-cbind(lo[,-c(seg1,seg2)],lo[,seg1],lo[,seg2]) ;
			}
			else
			{
				h<-c(h[-seg1],h[seg1]) ;
				lo<-cbind(lo[,-seg1],lo[,seg1]) ;
			}

			# get breakpoints
			bp<-array() ;
			bp[1]<-max(hli) ;
			bp[2]<-bp[1]+length(hri) ;
			bp[3]<-bp[2]+length(seg1) ;
			
			#li1<-(bh[3,] %in% lo_o1 ) ;
			#li2<-(bh[3,] %in% lo_o2) ;
			li1<-(bh[3,] %in% hri ) ;
			li2<-(bh[3,] %in% seg1) ;
			bh[3,li2]<-bh[3,li2]+length(hri) ;
			bh[3,li1]<-bh[3,li1]-length(seg1) ;
		
			
			#li1<-(snew %in% lo_o1) ;
			#li2<-(snew %in% lo_o2) ;
			li1<-(snew %in% hri) ;
			li2<-(snew %in% seg1) ;
			snew[li2]<-snew[li2]+length(hri) ;
			snew[li1]<-snew[li1]-length(seg1) ;
		
			
			# update bh at the breakpoints.
			
			for (j in seq(1,3))
			{
				bpi<- which(bh[3,]==bp[j]) ;
				if (length(bpi)==0)
				{
					bpi<-ncol(bh)+1 ;
					bh<-cbind(bh,matrix(0,3,1)) ;
					bh[3,bpi]<-bp[j] ;
				}
				if (bp[j]<ntips)
				{
					bh[1,bpi]<-lo[2,bp[j]] ;
					bh[2,bpi]<-lo[2,bp[j]+1] ;
					bh[1,bpi]<-hh[2,h[bp[j]]] ;
					bh[2,bpi]<-hh[2,h[bp[j]+1]] ;
				
					
				}
				else{
					bh<-bh[,-bpi] ;
				}
			}
	
			
			j<-which((bh[1,]==ha & bh[2,]==hb) |  (bh[1,]==hb & bh[2,]==ha)) ;
			
			e<-min(j) ;
			snew<-c(snew,bh[3,e]) ;
			lb<-h[bh[3,e]] ;
			rb<-h[bh[3,e]+1] ;
			lbc<-which(h==lb) ;
			rbc<-which(h==rb) ;
			h[lbc]<-Inode ;
			h[rbc]<-Inode ;
			Inode<-Inode-1 ;
			Inodes[2,i]<-T[i,2] ;
			bh<-cSIR_updatebh(lbc,rbc,bh,h,T[i,2]) ;
			bh<-matrix(bh[,-e],nrow=3) ;
		}
	}	
	return(list(s=snew,lo=lo)) ;
}