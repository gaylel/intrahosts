cSIR<-function(I0=1,nS=100,br=2,dr=1,maxiters=50)
{
	S<-nS-I0 ;
	I<-I0 ;
	T<-0 ;

	n<-1 ;
	while(I[n]>0 & n<maxiters) {
		r<-c((br/nS)*I[n]*S[n],dr*I[n]) ;
		
		# draw time interval
		T[n+1]<-rexp(1,sum(r)) + T[n];
		
		# draw event
		e<- sample(c(1,2),1,prob=r) ;
		
		I[n+1]<-I[n] ;
		S[n+1]<-S[n] ;
		if (e==1)
		{
			I[n+1]<-I[n+1] + 1;
			S[n+1]<-S[n+1] - 1;
		}
		else
		{
			I[n+1]<-I[n+1] - 1;
		}
		
		
		n<-n+1 ;
		
	}
	return(list(S=S,I=I,T=T)) ;
}

cSIR_multi<-function(nHosts=2,I0=1,nS=100,br=2,br2=1,dr=1,maxiters=50)
{
	# between host transmission rates
	B<-matrix(br2,nrow=nHosts,ncol=nHosts) ;
	for (i in seq(1,nHosts))
	{
		B[i,i]<-br ;
	}
	
	S<-matrix(nS,nrow=1,ncol=nHosts) ;
	S[1,1]<-nS-I0 ;
	
	I<-matrix(0,nrow=1,ncol=nHosts) ;
	I[1,1]<-I0 ;
	
	T<-matrix(0,nrow=1,ncol=4) ;
	T[1,1]<-0 ;

	R<-matrix(0,nrow=nHosts,ncol=nHosts+1) ;
	n<-1 ;
	
	while(any(I[n,]>0) & n<maxiters) {
		
		R[,nHosts+1]<-I[n,]*dr ;
		R[1:nHosts,1:nHosts]<-matrix(I[n,],nHosts,nHosts,byrow=TRUE)*matrix(S[n,],nHosts,nHosts)*B/nS ;
		#r<-c((br/nS)*I[n]*S[n],dr*I[n]) ;
		
		# draw time interval
		T<-rbind(T,rep(0,4)) ;
		T[n+1,1]<-rexp(1,sum(R)) + T[n];
		
		# draw event
		e<- sample(1:length(R),1,prob=as.vector(R)) ;
		ec<-ceiling(e/nHosts) ;
		er<-e-(ec-1)*nHosts ;
		
		I<-rbind(I,I[n,]);
		#I[n+1,]<-I[n,] ;
		
		S<-rbind(S,S[n,]) ;
		#S[n+1,]<-S[n,] ;
		
		if (ec==(nHosts+1))
		{
			I[n+1,er]<-I[n+1,er] - 1;
			T[n+1,4]<- -1 ;
			T[n+1,2:3]<-er ;
		}
		else
		{
		
			I[n+1,er]<-I[n+1,er] + 1;
			S[n+1,er]<-S[n+1,er] - 1;
			T[n+1,4]<-1 ;
			T[n+1,2]<-ec ;
			T[n+1,3]<-er ;
		}
		
		
		
		
		n<-n+1 ;
		
	}
	return(list(S=S,I=I,T=T)) ;
}

cSIR_tree_sample<-function(sir,N=100)
{
	I<-sir[[1]];
	S<-sir[[2]] ;
	T<-sir[[3]] ;
	nT<-nrow(T) ;
	nHosts<-ncol(I) ;
	g<-NULL ;
	el<-NULL ;

	nl<-matrix(1,nrow=1,ncol=nHosts) ;
	Nodes<-list() ;
	
	for (i in seq(1,nHosts))
	{
		Nodes[[i]]<-integer(0) ;
	}
	i<-which(I[1,]>0) ;
	nl[i]<-nl[i]+1;
	Nodes[[i]]<-paste(i,1,sep=":") ;
	
	for (t in seq(2,nT))
	{
		#if (T[t,4]==1)
		#{
		#	hb<-T[t,3] ;
		#	hind<-(hb-1)*N ;
		#	node_p<-sample(Nodes[[hb]],1) ;
		#	node_notp<-setdiff(Nodes[[hb]],node_p) ;
		#	m<-match(node_notp,g[2,]) ;
		#	el[m]<-el[m]+T[t,1]-T[t-1,1] ;
		#	g<-cbind(g,hind+c(node_p,nl[hb]));
		#	el<-cbind(el,T[t,1]-T[t-1,1]) ;
		#	Nodes[[hb]]<-c(Nodes[[hb]],nl[hb]) ;
		#	nl[hb]<-nl[hb]+1 ;
		#	g<-cbind(g,hind+c(node_p,nl[hb]));
		#	el<-cbind(el,T[t,1]-T[t-1,1]) ;
		#	Nodes[[hb]]<-c(Nodes[[hb]],nl[hb]) ;
		#	nl[hb]<-nl[hb]+1 ;
		#	node_notp<-setdiff(Nodes[[hb]],node_p) ;
		#	Nodes[[hb]]<-node_notp ;
		#}
		
		if (T[t,4]==1)
		{
			ha<-T[t,2] ;
			hb<-T[t,3] ;
			hbind<-(hb-1)*N ;
			haind<-(ha-1)*N ;
			node_p<-sample(Nodes[[ha]],1) ;
			node_notp<-setdiff(Nodes[[ha]],node_p) ;
			#m<-match(haind+node_notp,g[2,]) ;
			m<-match(node_notp,g[2,]) ;
			el[m]<-el[m]+T[t,1]-T[t-1,1] ;
			
			
			#g<-cbind(g,c(haind + node_p,hbind+nl[hb]));
			g<-cbind(g,c(node_p,paste(hb,nl[hb],sep=":")));
			el<-cbind(el,T[t,1]-T[t-1,1]) ;
			#Nodes[[hb]]<-c(Nodes[[hb]],nl[hb]) ;
			Nodes[[hb]]<-c(Nodes[[hb]],paste(hb,nl[hb],sep=":")) ;
			nl[hb]<-nl[hb]+1 ;
			
			#g<-cbind(g,haind+c(node_p,nl[ha]));
			g<-cbind(g,c(node_p,paste(ha,nl[ha],sep=":")));
			el<-cbind(el,T[t,1]-T[t-1,1]) ;
			Nodes[[ha]]<-c(Nodes[[ha]], paste(ha,nl[ha],sep=":"));
			#Nodes[[ha]]<-c(Nodes[[ha]],nl[ha]) ;
			nl[ha]<-nl[ha]+1 ;
			node_notp<-setdiff(Nodes[[ha]],node_p) ;
			Nodes[[ha]]<-node_notp ;
		}
		
		if (T[t,4]==-1)
		{
			hb<-T[t,3] ;
			hbind<-(hb-1)*N ;
			
			node_p<-sample(Nodes[[hb]],1) ;
			node_notp<-setdiff(Nodes[[hb]],node_p) ;
			#m<-match(hbind+node_notp,g[2,]) ;
			m<-match(node_p,g[2,]) ;
			el[m]<-el[m]+T[t,1]-T[t-1,1] ;
			m<-match(node_notp,g[2,]) ;
			el[m]<-el[m]+T[t,1]-T[t-1,1] ;
			Nodes[[hb]]<-setdiff(Nodes[[hb]],node_p) ;
		}
		
		
	}
	
	
	
	tre1<-list() ;
	tre1$edge<-t(g) ;
    tre1$edge.length<-el;
    Nnodes<-length(intersect(unlist(g[1,]),unlist(g[2,]))) ;
	tre1$Nnode<-Nnodes + 1;
	tips<-setdiff(unlist(g[2,]),unlist(g[1,])) ;
	tips<-sort(tips) ;
	tre1$tip.label<-as.character(tips) ;
	
	e<-matrix(0,nrow(tre1$edge),2) ;
	m<-match(tips,tre1$edge[,2]) ;
	e[m,2]<-seq(1,length(tips)) ;
	
	inode<-setdiff(c(unlist(g[2,]),unlist(g[1,])),tips) ;
	inode<-sort(inode) ;
	nlabel<-length(tips)+1 ;
	for (j in 1:length(inode))
	{
		m<-which(tre1$edge==inode[j]) ;
		e[m]<-nlabel ;
		nlabel<-nlabel+1 ;
	}
	tre1$edge<-e ;
	
	tlab<-tre1$tip.label ;
	nlab<-inode ;
	
	tlab2<-matrix(unlist(strsplit(tlab,":")),ncol=2,byrow=TRUE) ;
	nlab2<-matrix(unlist(strsplit(nlab,":")),ncol=2,byrow=TRUE) ;
	
	# sort out labels
	for (i in seq(1, nHosts))
	{
		j<-which(tlab2[,1]==as.character(i)) ;
		hname<-paste("H",i,sep="") ;
		tlab[j]<-paste(hname,paste("S",1:length(j),sep=""),sep=":")
		
		j<-which(nlab2[,1]==as.character(i)) ;
		nlab[j]<-as.character(i) ;
	}
	tre1$tip.label<-tlab ;
	tre1$node.label<-nlab ;
	
	
	class(tre1)<-"phylo" ;
	return(tre1);
	
}


cSIR_prunetree<-function(t,T)
{
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
	
	
	
	
	edges<-t$edge ;
	pr<-unique(c(which(d[edges[,2]]<T),sp)) ;
	edges<-edges[pr,] ;
	el<-t$edge.length[pr] ;
	nodes<-setdiff(edges[,2],tips2) ;
	nodes<-c(nodes,t$edge[sp,1]) ;
	nodes<-unique(nodes) ;
	no_ch<-setdiff(nodes,edges[,1]) ;
	while (length(no_ch)>0)
	{
	
		# nodes which have no children
	
	
		# get parents
		m<-match(nodes,t$edge[,2],nomatch=0);
		m<-m[m>0] ;
		nodes<-t$edge[m,1] ;
		nodes<-unique(nodes) ;
		
		# delete childless nodes
		m<-match(no_ch,edges[,2],nomatch=0) ;
		m<-m[m>0] ;
		edges<-edges[-m,] ;
		el<-el[-m] ;
		edges
		nodes
		no_ch<-setdiff(nodes,edges[,1]) ;
	
	}
	
	# collapse single nodes
	nc<-tabulate(edges[,1]) ;
	
	# single nodes
	
	s<-which(nc==1) ;
	if (length(s)>0)
	{
		for (i in seq(1,length(s)))
		{
			# get child of single node
			ch_i<-match(s[i],edges[,1],nomatch=0) ;
			if (ch_i>0)
			{
			ch<-edges[ch_i,2] ;
			
			# get parent of single node
			p_i<-match(s[i],edges[,2],nomatch=0) ;
			if (p_i>0)
			{
				pa<-edges[p_i,1] ;
				
				# join two edges together
				edges[p_i,2]<-ch ;
				el[p_i]<-el[p_i]+el[ch_i] ;
				edges<-edges[-ch_i,] ;
				el<-el[-ch_i] ;
			}
			else
			{
				# get grandchildren
				gc_i<-which(edges[,1]==ch) ;
				gc<-edges[gc_i,2] ;
				edges[gc_i,1]<-s[i] ;
				el[gc_i]<-el[gc_i]+el[ch_i] ;
				edges<-edges[-ch_i,] ;
				el<-el[-ch_i] ;
			}
			}
		}
	}
	# relabel
	ntips2<-length(tips2) ;
	m<-match(tips2,edges[,2]) ;
	
	
	nodes2<-setdiff(c(edges[,1],edges[,2]),tips2) ;
	edges2<-edges ;
	for (i in seq(1,length(nodes2)))
	{
		m2<-which(edges==nodes2[i]) ;
		edges2[m2]<- i+ntips2;
	}
	
	
	edges2[m,2]<-seq(1,ntips2) ;
	#nodelabel2<-t$node.label[nodes2-ntips] ;
	
	tiplabel2<-t$tip.label[tips2] ;
	Nnode2<-length(nodes2) ;
	t2<-list() ;
	t2$edge <- edges2 ;
	t2$Nnode<-Nnode2 ;
	#t2$node.label<-nodelabel2 ;
	t2$edge.length<-el ;
	t2$tip.label<-tiplabel2 ;
	class(t2)<-"phylo" ;
	return(t2) ;
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

cSIR_tree_sample2<-function(sir,N=100)
{
	I<-sir[[1]];
	S<-sir[[2]] ;
	T<-sir[[3]];
	nT<-nrow(T) ;
	nHosts<-ncol(I) ;
	g<-NULL ;
	el<-NULL ;

	nl<-matrix(1,nrow=1,ncol=nHosts) ;
	Nodes<-list() ;
	
	for (i in seq(1,nHosts))
	{
		Nodes[[i]]<-integer(0) ;
	}
	
	i<-which(I[1,]>0) ;
	node_p<-paste(i,nl[i],sep=":") ;
	nl[i]<-nl[i]+1;
	
	g<-rbind(node_p,paste(i,nl[i],sep=":")) ;
	
	for (j in i)
	{
	Nodes[[i[j]]]<-paste(i[j],nl[i[j]],sep=":") ;
	}
	nl[i]<-nl[i]+1;
	
	el<-rep(0,length(i)) ;
	
	for (t in seq(2,nT))
	{
		
		
		for (n in seq(1,nHosts))
		{
			m<-match(Nodes[[n]],g[2,]) ;
			el[m]<-el[m]+T[t,1]-T[t-1,1] ;
		}
		
		if (T[t,4]==1)
		{
			ha<-T[t,2] ;
			hb<-T[t,3] ;
			node_p<-sample(Nodes[[ha]],1) ;
			#node_notp<-setdiff(Nodes[[ha]],node_p) ;
			#m<-match(node_notp,g[2,]) ;
			
			#m<-match(Nodes[[ha]],g[2,]) ;
			#el[m]<-el[m]+T[t,1]-T[t-1,1] ;
			
			
			#g<-cbind(g,c(haind + node_p,hbind+nl[hb]));
			g<-cbind(g,c(node_p,paste(hb,nl[hb],sep=":")));
			el<-cbind(el,0) ;
			#Nodes[[hb]]<-c(Nodes[[hb]],nl[hb]) ;
			Nodes[[hb]]<-c(Nodes[[hb]],paste(hb,nl[hb],sep=":")) ;
			nl[hb]<-nl[hb]+1 ;
			
			#g<-cbind(g,haind+c(node_p,nl[ha]));
			g<-cbind(g,c(node_p,paste(ha,nl[ha],sep=":")));
			el<-cbind(el,0) ;
			Nodes[[ha]]<-c(Nodes[[ha]], paste(ha,nl[ha],sep=":"));
			#Nodes[[ha]]<-c(Nodes[[ha]],nl[ha]) ;
			nl[ha]<-nl[ha]+1 ;
			node_notp<-setdiff(Nodes[[ha]],node_p) ;
			Nodes[[ha]]<-node_notp ;
		}
		
		if (T[t,4]==-1)
		{
			hb<-T[t,3] ;
		
			
			node_p<-sample(Nodes[[hb]],1) ;
			#node_notp<-setdiff(Nodes[[hb]],node_p) ;
			#m<-match(hbind+node_notp,g[2,]) ;
			
			#m<-match(node_p,g[2,]) ;
			#el[m]<-el[m]+T[t,1]-T[t-1,1] ;
			#m<-match(node_notp,g[2,]) ;
			
			#m<-match(Nodes[[hb]],g[2,]) ;
			#el[m]<-el[m]+T[t,1]-T[t-1,1] ;
			Nodes[[hb]]<-setdiff(Nodes[[hb]],node_p) ;
		}
		
		
	}
	
	
	
	tre1<-list() ;
	tre1$edge<-t(g) ;
    tre1$edge.length<-el;
    Nnodes<-length(intersect(unlist(g[1,]),unlist(g[2,]))) ;
	tre1$Nnode<-Nnodes + 1;
	tips<-setdiff(unlist(g[2,]),unlist(g[1,])) ;
	tips<-sort(tips) ;
	tre1$tip.label<-as.character(tips) ;
	
	e<-matrix(0,nrow(tre1$edge),2) ;
	m<-match(tips,tre1$edge[,2]) ;
	e[m,2]<-seq(1,length(tips)) ;
	
	inode<-setdiff(c(unlist(g[2,]),unlist(g[1,])),tips) ;
	inode<-sort(inode) ;
	nlabel<-length(tips)+1 ;
	for (j in 1:length(inode))
	{
		m<-which(tre1$edge==inode[j]) ;
		e[m]<-nlabel ;
		nlabel<-nlabel+1 ;
	}
	tre1$edge<-e ;
	
	tlab<-tre1$tip.label ;
	nlab<-inode ;
	
	tlab2<-matrix(unlist(strsplit(tlab,":")),ncol=2,byrow=TRUE) ;
	nlab2<-matrix(unlist(strsplit(nlab,":")),ncol=2,byrow=TRUE) ;
	
	# sort out labels
	for (i in seq(1, nHosts))
	{
		j<-which(tlab2[,1]==as.character(i)) ;
		hname<-paste("H",i,sep="") ;
		tlab[j]<-paste(hname,paste("S",1:length(j),sep=""),sep=":")
		
		j<-which(nlab2[,1]==as.character(i)) ;
		nlab[j]<-as.character(i) ;
	}
	tre1$tip.label<-tlab ;
	tre1$node.label<-nlab ;
	
	
	class(tre1)<-"phylo" ;
	return(tre1);
	
}

cSIR_p_recover<-function(nS,br,dr,s)
{
	p<-dr/(dr + ((br/nS)*s)) ;
	return(p) ;
}

cSIR_tmx<-function(nHosts=1,I0=1,nS=5,br=2,br2=1,dr=1,maxiters=50)
{

	# construct the transition matrix
	#nstates<-(nS+1)*(nS+2)/2
	nstates<-(nS*(nS+1))/2 ;
	T<-matrix(0,nrow=nstates,ncol=nstates) ;
	for (i in seq(1,nS))
	{
		T[i,i]<-1;
	} 
	
	c1<-c(nS,2*nS -1) ;
	
	r1<-c(0,nS) ;
	r2<-c(2*nS-1, 2*nS-1 + nS-2);
	
	
	
	
	for (i in seq(nS,2,by=-1))
	{
		c1[1]<-c1[1]+1 ;
		r1[1]<-r1[1]+1 ;
		r2[1]<-r2[1]+1 ;
		A<-matrix(0,nrow=i,ncol=i-1) ;
		B<-matrix(0,nrow=i-2,ncol=i-1) ;
		
		if (i-2 >=1)
		{
		for (j in seq(i-2,1,by=-1))
		{
			p<-cSIR_p_recover(nS=nS,br=br,dr=dr,s=j) ;
			A[i-j,i-1-j]<-p ;
			B[i-j-1,i-j-1] <-(1-p) ;	
		}
		}
		A[i,i-1]<-cSIR_p_recover(nS=nS,br=br,dr=dr,s=0) ;
		if (r1[2] <=nstates && c1[2]<=nstates)
		{
			T[r1[1]:r1[2],c1[1]:c1[2]]<-A ;
		}
		
		if (r2[1] <=nstates )
		{
			T[r2[1]:r2[2],c1[1]:c1[2]]<-B ;
		}
		
		r1[1]<-r1[2] ;
		r2[1]<-r2[2] ;
		c1[1]<-c1[2] ;
		r1[2]<-r1[2]+i-1 ;
		r2[2]<-r2[2]+i-3 ;
		c1[2]<-c1[2]+i-2 ;
	}
	return(T) ;
}

cSIR_rmx<-function(nHosts=1,I0=1,nS=5,br=2,br2=1,dr=1,maxiters=50)
{

	# construct the transition matrix
	nstates<-(nS+1)*(nS+2)/2
	#nstates<-(nS*(nS+1))/2 ;
	Q<-matrix(0,nrow=nstates,ncol=nstates) ;
	
	c1<-c(nS,2*nS -1) ;
	
	r1<-c(0,nS) ;
	r2<-c(2*nS-1, 2*nS-1 + nS-2);
	
	
	
	
	for (i in seq(nS,2,by=-1))
	{
		c1[1]<-c1[1]+1 ;
		r1[1]<-r1[1]+1 ;
		r2[1]<-r2[1]+1 ;
		A<-matrix(0,nrow=i,ncol=i-1) ;
		B<-matrix(0,nrow=i-2,ncol=i-1) ;
		
		if (i-2 >=1)
		{
		for (j in seq(i-2,1,by=-1))
		{
			p<-cSIR_p_recover(nS=nS,br=br,dr=dr,s=j) ;
			A[i-j,i-1-j]<-p ;
			B[i-j-1,i-j-1] <-(1-p) ;	
		}
		}
		A[i,i-1]<-cSIR_p_recover(nS=nS,br=br,dr=dr,s=0) ;
		if (r1[2] <=nstates && c1[2]<=nstates)
		{
			T[r1[1]:r1[2],c1[1]:c1[2]]<-A ;
		}
		
		if (r2[1] <=nstates )
		{
			T[r2[1]:r2[2],c1[1]:c1[2]]<-B ;
		}
		
		r1[1]<-r1[2] ;
		r2[1]<-r2[2] ;
		c1[1]<-c1[2] ;
		r1[2]<-r1[2]+i-1 ;
		r2[2]<-r2[2]+i-3 ;
		c1[2]<-c1[2]+i-2 ;
	}
	return(T) ;
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
	sir<-list(I=sir1[[1]],S=sir1[[2]],T=sir1[[3]]) ;
	sir$T[,2:3]<-sir$T[,2:3]+1 ;
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
cSIR_sample<-function(nHosts=2,I0=1,NS=100,br=2,br2=1,dr=1,maxiters=50)
{
	sir<-cSIR(nHosts,I0,NS,br,br2,dr) ;
	tr<-cSIR_tree_sample2(sir,NS) ;
	#prtr<-cSIR_prunetree(tr,10) ;
	return(tr) ;
}





cSIRtest<-function(N=10)
{
	sir_all=list() ;
	for (i in seq(1,N))
	{
		sir<-cSIR(nHosts=1,nS=4,br=3) ;
		sir_all[[i]]<-cbind(sir[[1]],sir[[2]],sir[[3]][,c(1,4)]) ;
	}
	return(sir_all) ;
}

cSIRteststats<-function(sir)
{
	N<-length(sir) ;
	fs<-matrix(0,nrow=N,ncol=3) ;
	
	for (i in seq(1,N))
	{
		fs[i,]<-sir[[i]][nrow(sir[[i]]),c(1,2,3)] ;
		
	}
	
	
	return(fs) ;
}

cSIRteststats2<-function(sir,s0,i0)
{
	N<-length(sir) ;
	ts<-NULL ;
	
	for (i in seq(1,N))
	{
		a<-sir[[i]][,c(1,2,3)] ;
		inter<-intersect(which(a[,1]==i0),which(a[,2]==s0)) ;
		if (length(inter)>0)
		{ 
			
			ts<-rbind(ts,c(i,a[inter,])) ;
		}
	}
	return(ts) ;
	
	

}

cSIRteststats3<-function(sir,s,i)
{
	N<-length(sir) ;
	ts<-NULL ;
	
}

progenytest<-function(br1,br2,bd)
{
	l=length(bd) ;
	br<-c(br1,br2) ;
	for (i in seq(1,l))
	{
		# choose branch
		bi<-sample(c(1,2),1,prob=c(br[1]/sum(br),br[2]/sum(br))) ;
		
		br[bi]<-br[bi]+bd[i] ;
	}
	
	
	return(br) ;
}

progenytests<-function(br1,br2,bd,N)
{
	br<-NULL ; for (i in seq(1,N)){br<-rbind(br,progenytest(br1,br2,bd))}
    bins<-seq(0,max(br));
	pbr<-matrix(0,nrow=length(bins),2);
	for (i in bins)
	{
		pbr[i+1,1]<-length(which(br[,1]==i) );
		pbr[i+1,2]<-length(which(br[,2]==i) );	
	}	
	
	return(pbr/N);
	
}

progenydist<-function(br1,br2,bd)
{
	
	
	l<-length(bd) ;
	br<-c(br1,br2);
	lmax<-max(br1,br2)+length(which(bd>0))
	pbr<-matrix(0,nrow=1,ncol=lmax+1) ;
	p<-c(br[1]/sum(br),br[2]/sum(br)) ;
	pbr[br[1]+1]<-p[2] ;
	pbr[br[1]+bd[1]+1]<-p[1] ;
	
	if (l>1)
	{
	for (i in seq(2,l)) 
	{
		Ni<-br1+br2+sum(bd[1:(i-1)]) ;
		pbr2<-pbr ;
		for (c in seq(0,lmax)) 
		{
			pbr2[c+1]<-ifelse(c-bd[i] > 0 && c-bd[i]<=lmax,((c-bd[i])/Ni)*pbr[c+1-bd[i]],0);
			pbr2[c+1]<-pbr2[c+1]+(pbr[c+1]*(Ni-c)/Ni) ;
		}
		pbr<-pbr2 ; 
	}
	}
	return(pbr) ;
}

progenymxdist<-function(br1,br2,bd)
{
	l<-length(bd) ;
	br<-c(br1,br2) ;
	lmax<-max(br1,br2)+length(which(bd>0))
	pbr<-matrix(0,ncol=1,nrow=lmax+1) ;
	p<-c(br[1]/sum(br),br[2]/sum(br)) ;
	pbr[br[1]+1]<-p[2] ;
	pbr[br[1]+bd[1]+1]<-p[1] ;
	print(pbr) ;
	pmx2<-matrix(0,nrow=lmax+1,ncol=lmax+1) ;
	for (i in seq(1,lmax+1))
	{
		pmx2[i,i]=1 ;
	}
	
	if (l>1)
	{
	for (i in seq(2,l)) 
	{
		Ni<-br1+br2+sum(bd[1:(i-1)]) ;
		pmx<-matrix(0,nrow=lmax+1,ncol=lmax+1) ;
		for (m in seq(1, lmax+1))
		{
			pmx[m,m]<-(Ni-(m-1))/Ni ;
			if (m-1-bd[i] > 0 && m-1-bd[i]<=lmax)
			{
				pmx[m,m-bd[i]]<-(m-1-bd[i])/Ni ;
			}
		}
		pmx2<-pmx%*%pmx2;
		print(pmx) ;
		print(pmx2) ;
		pbr2<-pmx %*% pbr ;
		pbr<-pbr2 ; 
	}
	}
	return(pbr) ;
}

survivalsim<-function(lin_in,bd)
{
	lins<-NULL ;
	for (i in seq(1,length(lin_in)))
	{
		lins<-c(lins,rep(i,lin_in[i])) ;
	} 
	
	l<-length(bd) ;
	for (i in 1:l)
	{
		s<-sample(1:length(lins),1) ;
		if (bd[i]==1)
		{
			lins<-c(lins,lins[s]) ;
		}
		else
		{
			
			lins<-lins[-s] ;
		}
	}
	return(lins);
}
