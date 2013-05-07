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
	I<-sir$I ;
	S<-sir$S ;
	T<-sir$T ;
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
}

cSIR_tree_sample2<-function(sir,N=100)
{
	I<-sir$I ;
	S<-sir$S ;
	T<-sir$T ;
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
	nstates<-(nS+1)*(nS+2)/2
	#nstates<-(nS*(nS+1))/2 ;
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