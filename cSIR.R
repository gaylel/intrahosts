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
cSIR_eval<-function(sir,SN)
{

	p<-1;
	
	Iend<-sir$Iend ;
	NHosts<-length(Iend) ;
	i=1 ;
	while (p==1 & i<=NHosts)
	{
	 	#if (Iend[i]<SN[i] || Iend[i]>(10+SN[i])) p<-0 ;
	 	if (Iend[i]<SN[i] ) p<-0 ;
	 	i<-i+1 ;
	} 
	return(p) ;
}

cSIR_Keval<-function(sir,SN,h)
{

	p<-1;
	dis<-0 ;
	Iend<-sir$Iend ;
	NHosts<-length(Iend) ;
	 
	for (i in seq(1,NHosts))
	{
			dis=dis+(SN[i]-Iend[i])^2 ;
	}
	
	p<-exp(-dis/h) ;
	 
	return(p) ;
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

cSIRtree_reconstruct<-function(sir,N=100,SN,ST)
{
	I<-sir[[1]];
	S<-sir[[2]] ;
	T<-sir[[3]];
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
	T<-cSIR_treeparams()
	
	return(tr) ;
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
		if (T[t,4]< -1)
		{
			hb<-T[t,3] ;
		
			
			node_p<-sample(Nodes[[hb]],-T[t,4]) ;
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
	B<-exp(log(oldB)+rnorm(nHosts*nHosts,0,0.1)) ;
	return(B) ;
}

cSIR_drproposal<-function(olddr)
{
	dr<-exp(log(olddr)+rnorm(1,0,1)) ;
	return(dr) ;
}

cSIR_Bdrprior<-function(B,Bcon,dr,l1,l2,l3) 
{
	p<-0 ;
	nHosts<-nrow(B) ;
	# use exponential distributions
	for (i in seq(1,nHosts))
	{
		for (j in seq(1,nHosts))
		{
			if (Bcon[i,j]==1)
			{
			if (i==j)
			{
				p=p+(l1*exp(-l1*B[i,j]))			
			}
			else
			{
				p=p+(l2*exp(-l2*B[i,j])) ;
			}
			}
		}
	}
	p=p+(l3*exp(-l3*dr))
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

cSIR_SB_metrop<-function(nHosts=2,I0=1,nS=100,dr=1,ST=c(1,1),SN=c(30,30),N=1000)
{
	
	# initialise the parameters
	dr=0.1 ;
	B<-matrix(0.001*runif(nHosts*nHosts),nrow=nHosts,ncol=nHosts) ;
	for (i in seq(1,nHosts))
	{
		B[i,i]=1 ;
	}
	Bcon<-cSIR_Bstruct(ST=ST) ;
	
	B<-B*Bcon ;
	l1=1 ;
	l2=0.1;
	l3=0.1;
	Nacc=0 ;
	Nrej=0 ;
	iend<-NULL ;
	Bchain<-list() ;
	drchain<-list() ;
	sirchain<-list() ;
	trchain<-list() ;
	pold<-cSIR_Bdrprior(B,Bcon,dr,l1,l2,l3)  ;
	while (Nacc<N)
	{
		# make new proposal
		newdr<-cSIR_drproposal(dr) ;
		newdr<-0.1 ;
		newB<-cSIR_Bproposal(B) ;
	
		sir<-sample_cSIR_S_C(I0=I0, NS=nS, NHosts = nHosts, B=newB, dr=newdr, SN=SN, ST=ST) ;
		if (cSIR_eval(sir,SN)>0)
		{ 
			pnew<-cSIR_Bdrprior(newB,Bcon,newdr,l1,l2,l3) ;
			#print(pnew/pold)
			h<-min(c(1,pnew/pold)) ; 
			
			if (runif(1)<=h)
			{
			Nacc<-Nacc+1 ;
			iend<-rbind(iend,sir$Iend) ;
			drchain[[Nacc]]<-newdr ;
			Bchain[[Nacc]]<-newB ;
			sirchain[[Nacc]]<-sir ;
			dr<-newdr ;
			B<-newB ;
			pold<-cSIR_Bdrprior(B,Bcon,dr,l1,l2,l3)  ;
			print(sir$Iend) ;
		
			print("accept")
			tr<-cSIRtree_reconstruct(sir,N=nS,SN=SN,ST=ST) ;
			trchain[[Nacc]]<-tr ;
			}
		}
		else
		{
			Nrej<-Nrej+1 ;
			#print("reject") ;
		}
	}
	rList<-list(Nacc=Nacc,Nrej=Nrej,Iend=iend,B=Bchain,dr=drchain,sir=sirchain,tr=trchain) ;
	
}

cSIR_phylosample(sir,N,SN,ST)
{
	tr<-cSIRtree_reconstruct(sir,N=nS,SN=SN,ST=ST) ;
	cSIR_treeparams(tr
}

cSIR_SB_abcmetrop<-function(nHosts=2,I0=1,nS=100,dr=1,ST=c(1,1),SN=c(30,30),N=1000)
{
	
	# initialise the parameters
	dr=0.02 ;
	B<-matrix(0.001*runif(nHosts*nHosts),nrow=nHosts,ncol=nHosts) ;
	for (i in seq(1,nHosts))
	{
		B[i,i]=1 ;
	}
	l1=1 ;
	l2=0.1;
	l3=0.1;
	Nacc=0 ;
	Nrej=0 ;
	iend<-NULL ;
	Bchain<-list() ;
	drchain<-list() ;
	wchain<-list() ;
	pold<-cSIR_Bdrprior(B,dr,l1,l2,l3)  ;
	while (Nacc<N)
	{
		# make new proposal
		newdr<-cSIR_drproposal(dr) ;
		#newdr<-0.01 ;
		newB<-cSIR_Bproposal(B) ;
	
		sir<-sample_cSIR_S_C(I0=I0, NS=nS, NHosts = nHosts, B=newB, dr=newdr, SN=SN, ST=ST) ;
		print(sir$Iend) ;
		pacc<-cSIR_Keval(sir,SN,10)
		if (pacc>runif(1))
		{ 
			pnew<-cSIR_Bdrprior(newB,newdr,l1,l2,l3) ;
			Nacc<-Nacc+1 ;
			iend<-rbind(iend,sir$Iend) ;
			drchain[[Nacc]]<-newdr ;
			Bchain[[Nacc]]<-newB ;
			dr<-newdr ;
			B<-newB ;
			pold<-cSIR_Bdrprior(B,dr,l1,l2,l3)  ;
			print("accept")
			
			
		}
		else
		{
			
			Nrej<-Nrej+1 ;
			print("reject") ;
		}
	}
	rList<-list(Nacc=Nacc,Nrej=Nrej,Iend=iend,B=Bchain,dr=drchain) ;
	
}



cSIR_SB_sample<-function(nHosts=2,I0=1,nS=100,dr=1,ST=c(1,1),SN=c(30,30),N=1000)
{
	
	Nacc=0 ;
	Nrej=0 ;
	iend<-NULL ; 
	
	while (Nacc<N)
	{
		
		sir<-cSIR_SB(nHosts=nHosts,I0=I0,nS=nS,dr=dr,ST=ST,SN=SN) ;
		if (cSIR_eval(sir,SN)>0)
		{ 
			Nacc<-Nacc+1 ;
			iend<-rbind(iend,sir$Iend) ;
		}
		else
		{
			Nrej<-Nrej+1 ;
		}
	}
	rList<-list(Nacc=Nacc,Nrej=Nrej,Iend=iend,sir=sir) ;
	
}
cSIR_S_sample<-function(nHosts=2,I0=1,nS=100,br=2,br2=1,dr=1,ST=c(1,1),SN=c(30,30),N=1000)
{
	
	Nacc=0 ;
	Nrej=0 ;
	iend<-NULL ; 
	while (Nacc<N)
	{
		sir<-cSIR_S(nHosts=nHosts,I0=I0,nS=nS,br=br,br2=br2,dr=dr,ST=ST,SN=SN) ;
		if (cSIR_eval(sir,SN)>0)
		{ 
			Nacc<-Nacc+1 ;
			iend<-rbind(iend,sir$Iend) ;
		}
		else
		{
			Nrej<-Nrej+1 ;
		}
	}
	rList<-list(Nacc=Nacc,Nrej=Nrej,Iend=iend,sir=sir) ;
	
}

cSIR_ll<-function(nHosts=2,I0=1,nS=100,br=2,br2=1,dr=1,ST=c(1,1),SN=c(30,30))
{
	ll<-log(0) ;
	sir<-cSIR_S(nHosts=nHosts,I0=I0,nS=nS,br=br,br2=br2,dr=dr,ST=ST,SN=SN) ;
	if (cSIR_eval(sir,SN)>0)
	{ 
		ll<-log(1) ;
	}
	return(ll) ;
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
			
	}
	
	
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

cSIR_updatebh2<-function()
{
	# update bh after leaf shuffle
	
	
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
			
			
			lbi<-ifelse(length(lb)>1,sample(lb),lb) ;
			rbi<-ifelse(length(rb)>1,sample(rb),rb) ;
			
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
			li1<-(bh[3,] %in% lo_o1[1:(length(lo_o1)-1)] ) ;
			li2<-(bh[3,] %in% lo_o2[1:(length(lo_o2)-1)]) ;
			
			
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
			
			li1<-(bh[3,] %in% lo_o1 ) ;
			li2<-(bh[3,] %in% lo_o2) ;
			bh[3,li2]<-bh[3,li2]+length(hri) ;
			bh[3,li1]<-bh[3,li1]-length(seg1) ;
			
			li1<-(snew %in% lo_o1) ;
			li2<-(snew %in% lo_o2) ;
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
	return(list(bh=bh,s=snew,lo=lo)) ;
}