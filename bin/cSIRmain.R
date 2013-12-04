cSIR <- function(I0=1,nS=100,br=2,dr=1,maxiters=50)
{
  
  # Generates population trajectory following stochastic
  # continuous time SIR dynamics. 
  #
  # Args:
  #  I0: Initial number of infected
  #  nS: Population size
  #  br: birth/transmission rate
  #  dr: death/removal rate
  #  maxiters: maximum number of iterations.
  #
  # Returns:
  #  List containing:
  #
  #   S: trajectory of number of susceptibles
  #   I: trajectory of number of infected
  #   T: trajectory of event times
  
	
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

cSIR_multi <- function(nHosts=2,I0=1,nS=100,br=2,br2=1,dr=1,maxiters=50)
{
	# Generates population trajectory following stochastic
  # continuous time SIR dynamics, for multiple hosts 
  #
  # Args:
  #  nHosts: Number of hosts
  #  I0: Initial number of infected
  #  nS: Population size
  #  br: birth/transmission rate within host
  #  br2: birth/transmission rate between hosts
  #  dr: death/removal rate
  #  maxiters: maximum number of iterations.
  #
  # Returns:
  #  List containing:
  #   S: trajectory of number of susceptibles
  #   I: trajectory of number of infected
  #   T: trajectory of event timesB<-matrix(br2,nrow=nHosts,ncol=nHosts) 
	
  for (i in seq(1, nHosts))
  {
	  B[i, i] <- br 
  }
	S <- matrix(nS,nrow=1, ncol=nHosts) 
  S[1, 1] <- nS-I0 
  I <- matrix(0,nrow=1, ncol=nHosts) 
  I[1, 1] <- I0 
  T <- matrix(0,nrow=1, ncol=4) 
  T[1, 1] <- 0 
  R <- matrix(0,nrow=nHosts, ncol=nHosts+1) 
  n <- 1 
  while(any(I[n, ] > 0) & n < maxiters)
  {
    R[, nHosts + 1] <- I[n, ] * dr 
    R[1:nHosts, 1:nHosts] <- matrix(I[n, ], nHosts, nHosts, byrow=TRUE) * matrix(S[n, ], nHosts,nHosts)*B / nS 
	  # draw time interval
	  T <- rbind(T, rep(0,4)) 
	  T[n+1,1] <- rexp(1, sum(R)) + T[n]	
    # draw event
    e <- sample(1:length(R), 1, prob=as.vector(R)) 
    ec <- ceiling(e / nHosts) 
    er <- e - (ec - 1) * nHosts 
    I <- rbind(I, I[n, ])
    S <- rbind(S, S[n, ]) 
	  if (ec == (nHosts + 1))
    {
	    I[n+1, er] <- I[n+1, er] - 1
		  T[n+1, 4] <- -1 
		  T[n+1, 2:3] <- er 
		}else{
		  I[n+1, er] <- I[n+1, er] + 1
			S[n+1, er] <- S[n+1, er] - 1
			T[n+1, 4] <- 1 
			T[n+1, 2] <- ec 
			T[n+1, 3] <- er 
		}
		n<-n+1 
	}
	return(list(S=S,I=I,T=T)) 
}

cSIR_eval <- function(sir,SN)
{
  # Checks if the final infected population sizes in 
  # hosts are larger than number of observed samples
  #  
  # Args:
  #  sir: list object containing information about SIR dynamics
  #  SN: vector of population sizes for each host
  #
  # Returns:
  #  p: (1: valid SIR, 0: invalid SIR)
  
  p <- 1
	Iend <- sir$Iend 
	NHosts <- length(Iend) ;
	i <- 1 
	while (p == 1 & i <= NHosts)
	{
	 	if (Iend[i] < SN[i]) 
       p <- 0 
	 	i <- i + 1 
	} 
	return(p) 
}

cSIR_reconstructedtree<-function(tre, T)
{
  # Reconstructs phylogeny from matrix of birth and death events
  #  
  # Args:
  #  tre: 
  #  T: 
  #
  # Returns:
  #  
	
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

cSIRfull_tree <- function(sir)
{
  # generate possible underlying tree (including extinctions)
  T <- sir$T
  Iend <- sir$Iend
  print(Iend)
  ntips <- -sum(T[which(T[,4] < 0),4]) + sum(sir$I[nrow(sir$I),])
  print(ntips)
  nT <- nrow(T)
  e <- NULL
  el <- vector()
  nl <- vector()
  Inode <- ntips + 1
  Leafnode <- 1
  leaf_i <- -1
  Nodes <- matrix(c(Inode,1), nrow=2, ncol=1)
  NodeTimes <- rbind(seq(1,2*ntips-1), rep(0,2*ntips-1))
  #NodeTimes[ntips-1,2] <- 0
  for (i in seq(2, nT))
  {
    if (T[i,4] > 0)
    {
      for (j in seq(1,T[i,4]))
      {
        h <- which(Nodes[2,] == T[i,2])
        parent <- ifelse(length(h) > 1, sample(Nodes[1,h],1), Nodes[1,h]) 
        l <- which(e[, 2]==parent)
        if (length(l)>0)
        {
          e[l ,2] <- Inode
        }
        nl <- c(nl,T[i,2])
        Nodes <- matrix(Nodes[,-(which(Nodes[1,]==parent))],nrow=2)
        parent <- Inode
        
        e <- rbind(e,c(parent, leaf_i))
        e <- rbind(e,c(parent, leaf_i - 1))
        NodeTimes[parent] <- T[i, 1]
        Nodes <- cbind(Nodes, matrix(c(leaf_i, T[i,2]), nrow=2, ncol=1))
        Nodes <- cbind(Nodes, matrix(c(leaf_i - 1, T[i,3]), nrow=2, ncol=1))
        leaf_i <- leaf_i - 2
        
        Inode <- Inode + 1
      }
    }
    if (T[i,4] < 0)
    {
      
      for (j in seq(1, -T[i,4]))
      {
        h <- which(Nodes[2,] == T[i,2])
        
        leaf <- ifelse(length(h) > 1, sample(Nodes[1,h],1), Nodes[1,h]) 
        l <- which(e[, 2]==leaf)
        if (length(l)>0)
        {
          e[l ,2] <- Leafnode
          Nodes <- matrix(Nodes[,-(which(Nodes[1,]==leaf))],nrow=2)
          leaf <- Leafnode
          NodeTimes[Leafnode] <- T[i, 1]
          Leafnode <- Leafnode + 1
          
          
        }
        
        NodeTimes[leaf] <- T[i,1]
        
      }
    }
    #print(e)
    #print(Nodes)
    # construct the tree
  }
  
  Irem <- sir$Iend[nrow(sir$Iend),]
    for (i in seq(1,length(Irem)))
    {
      if (Irem[i] > 0)
      {
        for (j in seq(1, Irem[i]))
        {
          h <- which(Nodes[2,] == i)
          leaf <- ifelse(length(h) > 1, sample(Nodes[1,h],1), Nodes[1,h]) 
          l <- which(e[, 2]==leaf)
          if (length(l)>0)
          {
            e[l ,2] <- Leafnode
            Nodes <- matrix(Nodes[,-(which(Nodes[1,]==leaf))],nrow=2)
            leaf <- Leafnode
            NodeTimes[Leafnode] <- T[nrow(T), 1]
            Leafnode <- Leafnode + 1
          }
        
        }
      }  
    }
    
     
  
  
  for (i in seq(1,nrow(e)))
  {
    el <- c(el, NodeTimes[e[i,2]] - NodeTimes[e[i,1]])
  }
  tr <- list()
  tr$edge <- e
  tr$edge.length <- el
  tr$Nnode <- ntips - 1
  tr$node.label <- as.character(nl)
  tr$tip.label <- seq(1,ntips)
  class(tr) <- "phylo"
  return(list(tr=tr,nt=NodeTimes))
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
	
	# edge lengths
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

sample_cSIR_S_C <- function(I0, NS, NHosts, B, dr, SN, ST, bnprob)
{
  # Samples SIR population trajectory for given parameters and data
  #
  # Args:
  #   I0: initial number of infected in host 1
  #   NS: size of intrahost populations
  #   NHosts: number of hosts
  #   B: transmission matrix
  #   dr: death rate
  #   SN: vector of sample size in hosts
  #   ST: vector of sampling times for hosts
  # Returns:
  #   sir: list of:
  #     I:  infected size trajectory
  #     S:  susceptible size trajectory
  #     T:  matrix of events: T[i,]={time, from host, to host, [1=birth,-1=death]}
  #     Iend: vector of infected size for each host at times ST
  
	sir1 <- .Call("sample_cSIR_S_R", I0, NS, NHosts, B, dr, ST, SN, bnprob) ;
	sir <- list(I=sir1[[1]], S=sir1[[2]], T=sir1[[3]], Iend=sir1[[4]]) ;
	sir$T[,2:3] <- sir$T[,2:3] + 1 ;
	sir$B <- B ;
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
			#if (ST[j]-ST[i] >= 0 & ST[j]-ST[i] <=14)
			if (abs(ST[j]-ST[i]) <=14)
			{
			Bcon[i,j]=1 ;
			}
		}
	}

	return(Bcon) ;
}

cSIR_ne<-function(tr,SN,lo,s,mode)
{
	nHosts<-length(SN)
	nl<-NULL
	for (i in seq(1,nHosts))
	{
		nl<-c(nl,rep(i,SN[i])) ;
	}
	nl<-c(nl,as.integer(tr$node.label)) ;
	switch(as.character(mode),
		"1"={
			tr<-phylo_ne_2(tr,nl)
		},
		"2"={
			tr<-phylo_ne_3(tr,nl)
		}
	)
	res<-s_lo_from_phylo(tr)
	h<-match(res$lo,lo[1,])
	res$lo<-rbind(res$lo,lo[2,h])
	return(list(lo=res$lo,s=res$s,tr=tr))
}

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
	print(paste("aunt=",aunt,"niece=", niece))
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
	print(lo)
	print(s)
	return(list(tr=tr,lo=lo,s=s))
		
}	

cSIR_modelinit <- function(x, init, dat, dat.params, mcmc.params)
{
  # Initialises model with initialisation file
  #
  # Args:
  #   init: list of initial values
  #   dat: data info
  #   x: sequence data
  #
  # Returns:
  #   params
  
  
  
  dat$ST <- dat$STin + init$t_off
  ST <- dat$ST
  print(ST)
  SN <- dat$SN
  Bcon <- cSIR_Bstruct(ST=ST) 
  B <- cSIR_Binit(init, Bcon)
  dr <- cSIR_drinit(init) 
  mr <- cSIR_mrinit(init) 
  trinit <- cSIR_trinit(init, x, dat)
  print(trinit)
  llcur <- log(0)  
  params <- list(B=B, dr=dr, mr= mr, ll=llcur, Bcon=Bcon, tr_list=list(lo=trinit$lo, s=trinit$s), bn=init$bn, t_off=init$t_off, is.acc=0)

  tr_list <- cSIR_drawtr_list(params, dat.params, dat, mcmc.params$tr)
  
  params$tr_list <- tr_list
  params$is.acc <- tr_list$is.acc
  return(params)
}

cSIR_Bprior <- function(B=NULL, nHosts, hp1, hp2, mode="d", is.log=TRUE)
{
  # Prior on B
  # Args:
  #   B: current value of B
  #   hp: information about the prior
  #   mode: 
  #     "d": evaluates density for current value B
  #     "r": draws value of B from prior
  #   is.log: indicates whether log(d) is returned
  # Returns:
  #   rv: density or random value
  
  nHosts2 <- nHosts * nHosts
  intra <- seq(1, nHosts2, by=nHosts + 1)
  switch(mode,
    "d"={
      args <- hp1$args
      args$x <- B[intra]
      args$log <- is.log
      rv1 <- do.call(paste("d", hp1$dens, sep=""), args)       
      args <- hp2$args
      args$x <- B[-intra]
      args$log <- is.log
      rv2 <- do.call(paste("d", hp2$dens, sep=""), args)       
      # assume logs for now
      rv <- rv1 + rv2
    },
    "r"={
    }  
  )
  return(rv)
}

cSIR_mrprior <- function(mr=NULL, hp, mode="r", is.log=TRUE)
{
  # Prior on mr
  # Args:
  #   mr: current value of mr
  #   hp: information about the prior
  #   mode: 
  #     "d": evaluates density for current value mr
  #     "r": draws value of mr from prior
  #   is.log: indicates whether log(d) is returned
  # Returns:
  #   rv: density or random value
  
  args <- hp$args
  switch(mode,
    "d"={
      args$x <- mr
      args$log <- is.log
      rv <- do.call(paste("d", hp$dens, sep=""), args)       
    },
    "r"={
      args$n <- 1
      rv <- do.call(paste("r", hp$dens, sep=""), args)       
    }  
  )
  return(rv)
}

cSIR_bnprior <- function(bn=NULL, hp, mode="r", is.log=TRUE)
{
  # Prior on mr
  # Args:
  #   mr: current value of mr
  #   hp: information about the prior
  #   mode: 
  #     "d": evaluates density for current value mr
  #     "r": draws value of mr from prior
  #   is.log: indicates whether log(d) is returned
  # Returns:
  #   rv: density or random value
  
  args <- hp$args
  switch(mode,
         "d"={
           args$x <- bn
           args$log <- is.log
           rv <- do.call(paste("d", hp$dens, sep=""), args)       
         },
         "r"={
           args$n <- 1
           rv <- do.call(paste("r", hp$dens, sep=""), args)       
         }  
  )
  return(rv)
}


cSIR_drprior <- function(dr=NULL, hp, mode="r", is.log=TRUE)
{
  # Prior on dr
  # Args:
  #   dr: current value of dr
  #   hp: information about the prior
  #   mode: 
  #     "d": evaluates density for current value dr
  #     "r": draws value of mr from prior
  #   is.log: indicates whether log(d) is returned
  # Returns:
  #   rv: density or random value
  
  args <- hp$args
  switch(mode,
         "d"={
           args$x <- dr
           args$log <- is.log
           rv <- do.call(paste("d", hp$dens, sep=""), args)       
         },
         "r"={
           args$n <- 1
           rv <- do.call(paste("r", hp$dens, sep=""), args)       
         }  
  )
  return(rv)
}

cSIR_t_offprior <- function(t_off=NULL, hp, mode="r", is.log=TRUE)
{
  # Prior on t_off
  # Args:
  #   dr: current value of dr
  #   hp: information about the prior
  #   mode: 
  #     "d": evaluates density for current value dr
  #     "r": draws value of mr from prior
  #   is.log: indicates whether log(d) is returned
  # Returns:
  #   rv: density or random value
  
  args <- hp$args
  switch(mode,
         "d"={
           args$x <- t_off
           args$log <- is.log
           rv <- do.call(paste("d", hp$dens, sep=""), args)       
         },
         "r"={
           args$n <- 1
           rv <- do.call(paste("r", hp$dens, sep=""), args)       
         }  
  )
  return(rv)
}

cSIR_Binit <- function(init, Bcon)
{
  # Initialises transmission parameter matrix B
  #
  # Args:
  #   init: list of initial values
  #   Bcon: constraint on B matrix
  #
  # Returns:
  #   B: B initialised
  
  NHosts <- ncol(Bcon)
  B <- matrix(init$B$br2, nrow=NHosts, ncol=NHosts)
  for (i in seq(1, NHosts))
  {
    B[i,i] <- init$B$br1
  }
  B <- B * Bcon
  return(B)
}

cSIR_drinit <- function(init)
{
  # Initialises death rate dr
  #
  # Args:
  #   init: list of initial values
  #
  # Returns:
  #   dr: dr initialised
  
  dr <- init$dr
  return(dr)
}

cSIR_mrinit <- function(init)
{
  # Initialises mutation rate mr
  #
  # Args:
  #   init: list of initial values
  #
  # Returns:
  #   mr: mr initialised
  
  mr <- init$mr
  return(mr)  
}

cSIR_trinit <- function(init, x, dat)
{
  # Initialises tree structure: leaf order and branching order
  #
  # Args:
  #   init: list of initial values
  #   x: data
  #   dat: metadata
  # Returns:
  #   ini$lo  leaf order and host labels
  #   ini$s:  branching order
  
  SN <- dat$SN
  ST <- dat$ST
  
  switch(init$tr,
    "upgma"={
      ini <- cSIR_upgma(x, SN, ST)
    },
    "random"={
      ini <- list(s=NULL,lo=NULL)
      l2 <- 0
      
      for (i in seq(1,length(SN)))
      {
        l1 <- l2 + 1
        l2 <- l2 + SN[i]
        #lo <- ifelse(SN[i]>1, sample(seq(l1,l2)), l1)
        lo <- sample(seq(l1,l2))
        lo <- rbind(lo,rep(i,SN[i]))
        ini$lo <- cbind(ini$lo, lo)
      }
      ini$s <- sample(seq(1,sum(SN)-1))
    }
         
  )
  return(ini) 
}

cSIR_Bproposal <- function(oldB, Bcon, prop.params)
{
  # Propose a new value for B 
  #
  # Args:
  #   oldB: current value of B
  #   Bcon: constraint on B
  #   prop.params: parameters of the proposal:
  #     density
  #     type: 1=update all elements
  #           2=1 value for within host
  #           3=1 value for within host, and 1 for between hosts
  # Returns:
  #   B
  
  B <- oldB
  nHosts <- nrow(oldB) 
  nHosts2 <- nHosts * nHosts
  args <- prop.params$args
  switch(as.character(prop.params$type),
    "1"={
      args$n <- nHosts2
      m <- do.call(paste("r", prop.params$dens, sep=""), args)       
      B <- exp(log(oldB) + m)   
    },
    "2"={
      intra <- seq(1, nHosts2, by=nHosts + 1)
      Bii <- oldB[intra]
      Bii <- Bii[Bii>0][1]
      
      args$n <- nHosts2 - nHosts
      m <- do.call(paste("r", prop.params$dens, sep=""), args)       
      B[-intra] <- exp(log(oldB[-intra]) + m)
      args$n <- 1
      m <- do.call(paste("r", prop.params$dens, sep=""), args)       
      B[intra] <- exp(log(Bii) + m)
    },
    "3"={
      intra <- seq(1, nHosts2, by=nHosts + 1)
      args$n <- 2
      Bii <- oldB[intra]
      Bii <- Bii[Bii>0][1]
      Bij <- oldB[-intra]
      Bij <- Bij[Bij>0][1]
      intra <- seq(1, nHosts2, by=nHosts + 1)
      m <- do.call(paste("r", prop.params$dens, sep=""), args)
      B[-intra] <- exp(log(Bij) + m[1])
      B[intra] <- exp(log(Bii) + m[2])
    }
  )
  B <- B * Bcon
  return(B) ;
}

cSIR_drproposal <- function(olddr, prop.params)
{
  # Propose a new value for dr 
  #
  # Args:
  #   olddr: current value of dr
  #   prop.params: parameters of the proposal
  # Returns:
  #   dr
  
  args <- prop.params$args
  args$n <- 1
  m <- do.call(paste("r", prop.params$dens, sep=""), args)       
  dr <- exp(log(olddr) + m)
  return(dr) ;
}

cSIR_mrproposal <- function(oldmr, prop.params)
{
  # Propose a new value for mr 
  #
  # Args:
  #   oldmr: current value of mr
  #   prop.params: parameters of the proposal
  # Returns:
  #   mr
  
  args <- prop.params$args
  args$n <- 1
  m <- do.call(paste("r", prop.params$dens, sep=""), args)       
  mr <- exp(log(oldmr) + m)
  return(mr) ;
}

cSIR_t_offproposal <- function(oldt_off, prop.params)
{
  # Propose a new value for mr 
  #
  # Args:
  #   oldmr: current value of mr
  #   prop.params: parameters of the proposal
  # Returns:
  #   mr
  
  args <- prop.params$args
  args$n <- 1
  m <- do.call(paste("r", prop.params$dens, sep=""), args)       
  t_off <- exp(log(oldt_off) + m)
  return(t_off) ;
}

cSIR_bnproposal <- function(oldbn, prop.params)
{
  # Propose a new value for mr 
  #
  # Args:
  #   oldmr: current value of mr
  #   prop.params: parameters of the proposal
  # Returns:
  #   mr
  
  args <- prop.params$args
  args$n <- 1
  m <- do.call(paste("r", prop.params$dens, sep=""), args)      
  logitbn <- log(oldbn/(1-oldbn))
  bn <- exp(logitbn + m)
  bn <- bn / (bn+1)
  return(bn) ;
}

cSIR_choosemove <- function(moves, prob=NULL)
{
  # Picks a move at random
  #
  # Args:
  #   moves: vector of possible moves
  # Returns;
  #   m: move
  
  m <- ifelse(length(moves)>1,sample(moves, 1, prob=prob), moves)
  return(m)
}
  
cSIR_movesinit <- function(movevars, vars)
{
  # Creates a vector of the possible proposals
  # Args:
  #   movevars: which moves are allowed
  # Returns:
  #   moves: vector of move keys.
  
  moves <- NULL
  
  
  for (i in seq(1:length(vars)))
  {
    if (vars[i] %in% movevars)
    {
      moves <- c(moves, i)
    }
  }
  return(moves)
}

cSIR_chainsinit <- function(chainvars)
{
  # Initialises the structure to save samples
  # Args:
  #   chainvars: names of samples to save
  # Returns:
  #   ch: list of sample chains
  
  ch <- list()
  for (i in seq(1,length(chainvars)))
  {
      ch[[chainvars[i]]] <- list()    
  }
  
  return(ch)
}

cSIR_chainsupdate <- function(ch, chainvars, params)
{
  n <- length(ch[[1]])
  for (i in seq(1,length(chainvars)))
  {
      if (chainvars[i]=="tr")
      {
        ch$tr[[n+1]] <- params$tr_list$tr 
      }else{
        ch[[chainvars[i]]][[n+1]] <- params[[chainvars[i]]]
      }
  }
  return(ch)
}

cSIR_runmcmc <- function(x, dat, opt, init, mcmc.params, hp.params)
{
	# infer posterior parameters
  # initialise the parameters
  
  SN<-dat$SN 
  ST<-dat$ST 
  dat.params<-list(I0=init$I0, NS=init$NS, NHosts=length(SN))
  
  
  
  dat$STin = dat$ST - min(dat$ST)
    
  
  vars <- c("mr","dr", "T", "B", "T2", "B_T", "bn", "t_off", "T3")
  
  params <- cSIR_modelinit(x, init, dat, dat.params, mcmc.params)
  moves <- cSIR_movesinit(opt$movevars, vars) 
  prob <- rep(0.5/8, 8)
  prob[3] <- 0.5 
	ch <- cSIR_chainsinit(opt$chainvars)
  ch <- cSIR_chainsupdate(ch, opt$chainvars, params)
  
  # convert data 
  x <- phyDat(x) 
  if (mcmc.params$acc.rate==1)
  {
    # calculate acceptance rates
    ch$acc.rate <- list()
    oldacc.rate <- list()
    acc.rate <- list()
    for (i in seq(1, length(opt$movevars)))
    {
      ch$acc.rate[[opt$movevars[i]]] <- NULL
      oldacc.rate[[opt$movevars[i]]] <- NULL
      acc.rate[[opt$movevars[i]]] <- 1
    }
    print(acc.rate)
  }  
    
  
  t <- 1
	while (t <= mcp$Niters)
	{
		# make new proposal
		m <- cSIR_choosemove(moves, prob)

   
		switch(as.character(m),
		  "1"={
		    params <- cSIR_mrupdate(x, params, dat.params, dat, hp.params, mcmc.params$mr)
		  },
		  "2"={
		    params <- cSIR_drupdate(params, dat.params, dat, hp.params, mcmc.params)
		  },
		  "3"={
		    params <- cSIR_Tupdate(x, params, dat.params, dat, hp.params, mcmc.params$tr)
		  },
		  "4"={ 
		    params <- cSIR_Bupdate(params, dat.params, dat, hp.params, mcmc.params)
    	},
		  "5"={
		    params <- cSIR_T2update(x, params, dat.params, dat, hp.params, mcmc.params$tr, 1)
        
		  },
		  "6"={ 
		    params <- cSIR_B_T_update(x, params, dat.params, dat, hp.params, mcmc.params)
		  },
		  "7"={ 
		    params <- cSIR_bnupdate(params, dat.params, dat, hp.params, mcmc.params)
		  },
		  "8"={ 
		    params <- cSIR_t_offupdate(params, dat.params, dat, hp.params, mcmc.params)
		  },
		  "9"={
		    params <- cSIR_T2update(x, params, dat.params, dat, hp.params, mcmc.params$tr, 2)
        
		  }
		)
		ch <- cSIR_chainsupdate(ch, opt$chainvars, params)
    if (mcmc.params$acc.rate==1)
    {
      v <- vars[m]
      if (length(acc.rate[[v]])<100)
      {
        acc.rate[[v]] <- c(acc.rate[[v]], params$is.acc)
      }else{
        acc.rate[[v]] <- c(acc.rate[[v]][2:100], params$is.acc)
      }  
      
      for (i in seq(1, length(opt$movevars)))
      {
       
        ch$acc.rate[[opt$movevars[i]]] <- c(ch$acc.rate[[opt$movevars[i]]], oldacc.rate[[opt$movevars[i]]])
        #oldacc.rate[[opt$movevars[i]]] <- ch$acc.rate[[opt$movevars[i]]][length()] 
      }
      oldacc.rate[[v]] <- mean(acc.rate[[v]])
      ch$acc.rate[[v]] <- c(ch$acc.rate[[v]], mean(acc.rate[[v]]))

    }
    
    print(sprintf("%i loglikelihood   %8.4f   Mutation rate   %8.4f  Death rate %8.4f BN %8.4f t_off %8.4f Move %i", t, params$ll, params$mr, params$dr,params$bn, params$t_off, m))
	
    
		if (opt$saveevery > 0 & ((t %% opt$saveevery) == 0))
		{
		  vars2 <- opt$chainvars
		  for (i in seq(1,length(vars2)))
		  {
		    assign(vars2[i],ch[[vars2[i]]])
		    save(list=vars2[i],file=paste(opt$outdir,"/",vars2[i], "_", t,".mcmc",sep="")) ;
		  }
		  if (mcmc.params$acc.rate==1)
		  {
		    acc_rate <- ch$acc.rate
		    save(acc_rate, file=paste(outdir,"/acc.rate", "_", t,".mcmc",sep=""))
		  }  
      
      for (i in seq(1,length(ch)))
      {
        ch[[i]] <- list()
      }
		}
		t <- t + 1
	}
	
	#rList<-list(Nacc=Nacc,Nrej=Nrej,B=Bchain,dr=drchain,mr=mrchain,tr=trchain,ll=llchain) ;
	return(ch)
}

cSIR_demo_1<-function()
{
	B<-matrix(0,nrow=3,ncol=3)
	for (i in seq(1,3))
	{
		B[i,i]<-1
	}
	B[2,3]<-0.01
	B[1,2]<-0.01
	B[2,1]<-0
	I0<-1
	NS<-100
	NHosts<-3
	dr<-0.6
	SN<-c(10,10,10)
	ST<-c(7,9,11)
	Nacc<-0
	while(Nacc==0)
	{
		sir<-sample_cSIR_S_C(I0=I0, NS=NS, NHosts = NHosts, B=B, dr=dr, SN=SN, ST=ST, 0.1) 
		
		if (cSIR_eval(sir,SN)>0)
		{Nacc<-1}
	}
	return(sir)
}

cSIR_drawtr_list <- function(params, dat.params, dat, mcmc.params)
{
  # Draws ancestral tree of infected intrahost populations, based on SIR dynamics
  # Args:
  #   params: current parameter estimates of the model
  #   dat.params: parameter estimates of the data
  #   dat:  information about the data
  #   mcmc.params: parameters of the optimisation of tr
  #
  # Returns:
  #     tr_list:  list of phylo parameters:
  #       tr: ancestral phylo object 
  #       T:  matrix of branching times T[i,] = {time, host of parent, host of child 1, host of child 2}
  #       lo: ordering of the tips of new tr
  #       s: coalescent order of new tr
  #       is.acc: 1: new phylo has been generated, 0: new phylo hasn't been generated
  
  tr_list <- params$tr_list 
  dat$ST <- dat$STin  + params$t_off
  sir <- cSIR_trdraw1(params, dat.params, dat, mcmc.params)
  if (sir$is.acc==1)
  {
    tr_list <- cSIR_trdraw2(sir, dat, dat.params)

    tr_list <- cSIR_trdraw3(tr_list$T, dat, params$tr_list)

  }

  tr_list$is.acc <- sir$is.acc
  return(tr_list)
}

cSIR_trdraw1 <- function(params, dat.params, dat, mcmc.params)
{
  # draw valid sir process with non-zero probability of generating data
  # Args:
  #   params: current parameter estimates of the model
  #   dat.params: parameter estimates of the data
  #   dat:  information about the data
  #   mcmc.params: parameters of the optimisation for tr
  # Returns:
  #   sir: list of:
  #     I:  infected size trajectory
  #     S:  susceptible size trajectory
  #     T:  matrix of events: T[i,]={time, from host, to host, [1=birth,-1=death]}
  #     Iend: vector of infected size for each host at times ST  
  #     is.acc: 1=new sir sample accepted, 0= not accepted
  
  I0 <- dat.params$I0
  NS <- dat.params$NS
  NHosts <- dat.params$NHosts
  B <- params$B
  dr <- params$dr
  bnprob <- params$bn
  SN <- dat$SN
  st <- dat$ST
  t <- 0
  is.acc <- 0
  
  while (t < mcmc.params$Ntries_sir & is.acc==0 )
  {
    # draw SIR sample
    sir <- sample_cSIR_S_C(I0=I0, NS=NS, NHosts = NHosts, B=B, dr=dr, SN=dat$SN, ST=dat$ST, bnprob)
    if (cSIR_eval(sir,dat$SN)>0)
    {
      # if accepted
      is.acc <- 1 
    }
    t <- t + 1
    
  }
  sir$is.acc <- is.acc
  return(sir)
}

cSIR_trdraw2 <- function(sir, dat, dat.params)
{
  # get branching times for reconstructed process based on birth-death process event times
  # Args:
  #   sir: list of:
  #     I:  infected size trajectory
  #     S:  susceptible size trajectory
  #     T:  matrix of events: T[i,]={time, from host, to host, [1=birth,-1=death]}
  #     Iend: vector of infected size for each host at times ST  
  #     is.acc: 1=new sir sample accepted, 0= not accepted
  # Returns:
  #   tr_list: list of:  
  #     tr: phylo object
  #     T: matrix of branching times T[i,] = {time, host of parent, host of child 1, host of child 2}
    
  tr_list <- .Call("tree_reconstruct", sir, dat.params$NHosts, dat)    
  names(tr_list) <- c("tr", "T")
  return(tr_list)
}

cSIR_trdraw3 <- function(T, dat, tr.params)
{
  # get reconstructed tree based on input topology and branching times
  # Args:
  #     T:  matrix of events: T[i,]={time, from host, to host, [1=birth,-1=death]}
  #     dat: data parameters
  #     tr.params: list of phylo parameters:
  #       lo: ordering of the tips
  #       s: coalescent order
  #
  # Returns:
  #     tr_list:  list of phylo parameters:
  #       tr: phylo object 
  #       T:  matrix of branching times T[i,] = {time, host of parent, host of child 1, host of child 2}
  #       lo: ordering of the tips of new tr
  #       s: coalescent order of new tr
  
  lo <- tr.params$lo
  s <- tr.params$s
  bt <- cSIR_stimes(T=T,lo=lo,s=s) 
  tr <- cSIR_phylofroms(s=bt$s, lo=bt$lo, T=T, dat=dat) 
  tr$edge.length[tr$edge.length<0] <- 1e-10
  tr_list <- list(tr=tr, T=T, lo=bt$lo, s=bt$s)

  return(tr_list)
}

cSIR_trdraw4 <- function()
{
  # get statistics of tree
  #tre$TS<-cSIR_getTstats(tre$T, NHosts)
  #tre$T<-diff(tre$T[,1])
}

cSIR_Tupdate<-function(x, params, dat.params, dat, hp.params, mcmc.params)
{
  # Updates ancestral tree of infected intrahost populations, based on SIR dynamics
  # Args:
  #   params: current parameter estimates of the model
  #   dat.params: parameter estimates of the data
  #   dat:  information about the data
  #   mcmc.params: parameters of the optimisation 
  #
  # Returns:
  #      
  #
  
	tr_list <- cSIR_drawtr_list(params=params, dat.params=dat.params, dat=dat, mcmc.params=mcmc.params)
	params$is.acc <- 0
  if (tr_list$is.acc==1)
	{
	  llcur<-params$ll
    print(min(tr_list$tr$edge.length))
    #print(tabulate(tr_list$tr$edge[,1]))
    #print(tabulate(tr_list$tr$edge[,2]))
    #print(tr_list$tr$edge)
	  ll <- pml(tr_list$tr, x, rate=params$mr, model=hp.params$mut$model)  
	  print(ll)
	  pacc <- min(ll$logLik - llcur, 0)			
	  if (log(runif(1)) <= pacc)
	  {
	    params$tr_list <- tr_list
	    params$ll <- ll$logLik
	    params$is.acc <- 1
	  }  
	}
  
	return(params)
}

cSIR_T2update<-function(x, params, dat.params, dat, hp.params, mcmc.params, mode=1)
{
  # Updates ancestral tree of infected intrahost populations, based on narrow exchange
  # Args:
  #   params: current parameter estimates of the model
  #   dat.params: parameter estimates of the data
  #   dat:  information about the data
  #   mcmc.params: parameters of the optimisation 
  #
  # Returns:
  #   params
  # narrow exchange

  res <- cSIR_ne(params$tr_list$tr, dat$SN, params$tr_list$lo, params$tr_list$s, mode=mode)
  lo <- res$lo 
  s <- res$s 
  tr <- res$tr 
  ll <- pml(res$tr, x, rate=params$mr, model=hp.params$mut$model)  
  newll <- ll$logLik 
  pacc <- min(ll$logLik - params$ll, 0) 
  params$is.acc <- 0
  rv <- log(runif(1))
  #print(sprintf("pacc=%8.8f, rv=%8.8f, newll=%8.8f, oldll=%8.8f",pacc, rv, newll, params$ll))
  
  if (rv < pacc)
  {
    params$tr_list$lo <- lo
    params$tr_list$s <- s
    
    params$tr_list$tr <- tr
    
    params$ll <- ll$logLik
    params$is.acc <- 1
  }
  
  return(params)
}

cSIR_B_T_update<-function(x, params, dat.params, dat, hp.params, mcmc.params)
{
	newparams<-cSIR_Bupdate(params, dat.params, dat, hp.params, mcmc.params)
	newB<-newparams$B
	newparams2 <- cSIR_Tupdate(x, params, dat.params, dat, hp.params, mcmc.params$tr)
	if (newparams2$ll > params$ll)
	{
		newparams2$B<-newparams$B
    newparams2$is.acc <- 1
		print("T,B accept")
	}
	else
	{
		newparams2$B<-params$B
		newparams2$is.acc <- 0
	}
	params<-newparams2
	return(params)
}

cSIR_Bupdate<-function(params, dat.params, dat, hp.params, mcmc.params)
{
	newB <- cSIR_Bproposal(params$B, params$Bcon, mcmc.params$B)  
	oldB<-params$B
	params$B<-newB
	tr_list <- cSIR_drawtr_list(params=params, dat.params=dat.params, dat=dat, mcmc.params=mcmc.params$tr)
	params$B<-oldB
  params$is.acc <- 0
  
  if (tr_list$is.acc == 1)
  {
    tstats <- cSIR_getTstats(tr_list, params$tr_list, mcmc.params)
    print(tstats)
    th <- mcmc.params$abc$th
    if (tstats$d1<th & tstats$d2<th)
    {
      pold <- cSIR_Bprior(B=oldB, dat.params$NHosts, hp.params$B.br1, hp.params$B.br2, mode="d", is.log=TRUE)
      pnew <- cSIR_Bprior(B=newB, dat.params$NHosts, hp.params$B.br1, hp.params$B.br2, mode="d", is.log=TRUE)
      h<-min(1,exp(pnew - pold)) ;
      if (runif(1)<=h)
      {
        params$B<-newB
        params$is.acc <- 1
      }
    }
    
  }
	return(params)			
}

cSIR_bnupdate<-function(params, dat.params, dat, hp.params, mcmc.params)
{
  newbn <- cSIR_bnproposal(params$bn, mcmc.params$bn)  
  oldbn<-params$bn
  params$bn<-newbn
  tr_list <- cSIR_drawtr_list(params=params, dat.params=dat.params, dat=dat, mcmc.params=mcmc.params$tr)
  params$bn<-oldbn
  params$is.acc <- 0
  
  if (tr_list$is.acc == 1)
  {
    tstats <- cSIR_getTstats(tr_list, params$tr_list, mcmc.params)
    th <- mcmc.params$abc$th
    if (tstats$d1<th & tstats$d2<th)
    {
      pold <- cSIR_bnprior(bn=oldbn, hp.params$bn, mode="d", is.log=TRUE) 
      pnew <- cSIR_bnprior(bn=newbn, hp.params$bn, mode="d", is.log=TRUE) 
      h<-min(1,exp(pnew - pold)) 
      if (runif(1)<=h)
      {
        params$bn<-newbn
        params$is.acc <- 1
      }
    }
    
  }
  return(params)			
}

cSIR_drupdate <- function(params, dat.params, dat, hp.params, mcmc.params)
{
  newdr <- cSIR_drproposal(params$dr, mcmc.params$dr)  
  olddr<-params$dr
  params$dr<-newdr
  tr_list <- cSIR_drawtr_list(params=params, dat.params=dat.params, dat=dat, mcmc.params=mcmc.params$tr)
  params$dr<-olddr
  params$is.acc <- 0
  
  if (tr_list$is.acc == 1)
  {
    tstats <- cSIR_getTstats(tr_list, params$tr_list, mcmc.params)
    th <- mcmc.params$abc$th
    if (tstats$d1<th & tstats$d2<th)
    {
      pold <- cSIR_drprior(dr=olddr, hp.params$dr, mode="d", is.log=TRUE)
      pnew <- cSIR_drprior(dr=newdr, hp.params$dr, mode="d", is.log=TRUE)
      h<-min(1,exp(pnew - pold)) ;
      if (runif(1)<=h)
      {
        params$dr<-newdr
        params$is.acc <- 1
      }
    }
    
  }
  return(params)			
}

cSIR_mrupdate <- function(x, params, dat.params, dat, hp.params, mcmc.params)
{
  # Updates mutation rate mr
  #
  # Args:
  #   x: sequence data
  #   params: current model parameters
  #   dat.params: other parameters
  #   dat: data info
  #   prop.params:  proposal parameters for mr
  #   hp.params: parameters for the priors
  # Returns:
  #   params
  #

  oldmr <- params$mr
  oldp <- cSIR_mrprior(mr=oldmr, hp=hp.params$mr, mode="d", is.log=TRUE)
  newmr <- cSIR_mrproposal(params$mr, mcmc.params) 
  newp <- cSIR_mrprior(mr=newmr, hp=hp.params$mr, mode="d", is.log=TRUE)
  params$mr <- newmr
  tr <- params$tr_list$tr
	ll<-pml(tr, x, rate=params$mr, model=hp.params$mut$model)  
	newll <- ll$logLik
  params$mr <- oldmr
  pacc=min(newll - params$ll + oldp - newp,0)	
  params$is.acc <- 0
	if (log(runif(1)) <= pacc)
	{
		params$mr <- newmr
		params$ll <- newll
    params$is.acc <- 1
	}
	
	return(params)	
}

cSIR_t_offupdate <- function(params, dat.params, dat, hp.params, mcmc.params)
{
  newt_off <- cSIR_t_offproposal(params$t_off, mcmc.params$t_off)  
  oldt_off<-params$t_off
  params$t_off<-newt_off
  tr_list <- cSIR_drawtr_list(params=params, dat.params=dat.params, dat=dat, mcmc.params=mcmc.params$tr)
  params$t_off<-oldt_off
  params$is.acc <- 0
  
  if (tr_list$is.acc == 1)
  {
    tstats <- cSIR_getTstats(tr_list, params$tr_list, mcmc.params)
    th <- mcmc.params$abc$th
    if (tstats$d1<th & tstats$d2<th)
    {
      pold <- cSIR_t_offprior(t_off=oldt_off, hp.params$t_off, mode="d", is.log=TRUE)
      pnew <- cSIR_t_offprior(t_off=newt_off, hp.params$t_off, mode="d", is.log=TRUE)
      h<-min(1,exp(pnew - pold)) ;
      if (runif(1)<=h)
      {
        params$t_off<-newt_off
        params$is.acc <- 1
      }
    }
    
  }
  return(params)  		
}
cSIR_phylosample<-function(sir,nS,lo,s,dat)
{
  #     
  #  
  # 
  #   
  #
  
	tr<-cSIRtree_reconstruct(sir,N=nS,SN=dat$SN,ST=dat$ST) ;
	NHosts<-length(dat$SN) ;
	#xprint(sir$T)
	
	
	#tr_list<-.Call("tree_reconstruct",sir, NHosts, dat) ;
	#tr<-tr_list[[1]]


	#T=tr_list[[2]]
	#print(T)
	Tp<-cSIR_treeparams(tr,SN=dat$SN,ST=dat$ST) ;
	T=Tp$T
	
	# get new branching order based on current and leaf order
	bt<-cSIR_stimes(T=T,lo=lo,s=s) ;
	
	
	print(T)
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

cSIR_plottraj<-function(sir, cols, fname)
{
	# plot infected trajectories
	I<-sir$I ;
	T<-sir$T ;
  
	nHosts<-ncol(I) ;
	mI<-max(I)+10 ;
  xlim=c(0,max(T[,1]))
	#pdf(fname)
  postscript(file=fname, width=10, height=5, pointsize=12 , horizontal = FALSE, onefile = FALSE, paper = "special")
  plot(T[I[,1]>0,1],I[I[,1]>0,1],xlim=xlim, ylim=c(0,mI),type="s", col=cols[1],xlab=expression(t), ylab=expression(I(t)), lwd = 3, cex=2, axes=FALSE) ;
  axis(1)
  axis(2)
  
	for (i in seq(2,nHosts))
	{
		lines(T[I[,i]>0,1],I[I[,i]>0,i],type="s", col=cols[i], lwd = 3) ;
	}
  dev.off()
 
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

cSIR_getTstats<-function(tr_list1, tr_list2, mcmc.params)
{
	s1_1 <- cSIR_trstat1(tr_list1)
    s1_2 <- cSIR_trstat1(tr_list2)
	s2_1 <- diff(tr_list1$T[,1])
	s2_2 <- diff(tr_list2$T[,1])
	
  d1<-sum((s1_1 - s1_2)^2) 
	d2<-sum((s2_1 - s2_2)^2)/(length(s2_1)-1) 
	return(list(d1=d1, d2=d2))
}

cSIR_trstat1 <- function(tr_list)
{
  NHosts <- max(tr_list$lo[2,])
  T <- tr_list$T
  N_T <- nrow(T)
  TS <- matrix(0, nrow=NHosts, ncol=NHosts)
  for (i in seq(1, N_T))
  {
    TS[T[i, 2], T[i, 3]] <- TS[T[i, 2], T[i, 3]] + 1
    TS[T[i, 2], T[i, 4]] <- TS[T[i, 2], T[i, 4]] + 1
    
  }
  return(TS / sum(TS))
}

cSIR_savesmps <- function(opt, smp)
{
  vars<-c("mr","dr","B","ll","tr")
  for (i in seq(1,length(vars)))
  {
    assign(vars[i],smp[[vars[i]]])
    save(list=vars[i],file=paste(outdir,"/",vars[i],".mcmc",sep="")) ;
  }  
}