reed_frost<-function(I0=1,N=100, q=0.5,niters=100)
{
  I<-I0 ;
  S<-N-I ;
  R<-0 ;
  
  for (n in seq(1,niters))
  {
  p<-1-q^I ;
  I<-rbinom(1,S,p) ; 
  S<-S-I ;
  print(I) ;
  }
}



reed_frost_multi<-function(I0=1, N=100, NHosts=2, q0=0.99, q1=0.99)
{
  # evolution of viral population sizes in multiple hosts
  niters<-50;
  I<-matrix(0,nrow=50,ncol=NHosts) ;
  I[1,1]<-I0 ;
  It<-matrix(0,nrow=NHosts, ncol=NHosts+1) ;
  
  S <- matrix(N,nrow=1,ncol=NHosts) ;
  S[1,1] <- N - I[1,1] ;
  p <- matrix(0,nrow=1,ncol=NHosts+1) ;
 
  for (n in seq(2,niters))
  {
       for (j in seq(1,NHosts))
       {
	 p[1,NHosts+1]=1 ;
	 for ( i in seq(1,NHosts))
	  {
	      # probability of infection between i and j
	      q<-ifelse(i==j,q0,q1) ; 
	      p[1,i]<-1-q^I[n-1,i] ;
	      p[1,NHosts+1] = p[1,NHosts+1]*(q^I[n-1,i]) ;
	  }
	  p<-p/sum(p) ;
	  #print(p) ;
	  It[j,]<-rmultinom(1,S[1,j],p) ; 
	  
	  I[n,j]<- sum(It[j,1:NHosts]);
       }
       

       
       for ( i in seq(1,NHosts))
       {
	  S[1,i]<-It[i,NHosts+1];
       }
       
       #print(S[1,]);
       #print(I[n,]) ;
       #print(It);
  }
  return(list(I=I));
}

reed_frost_multi.2<-function(I0=1, N, NHosts, q0=0.999, q1=0.9999)
{
  # evolution of viral population sizes in multiple hosts
  niters<-50;
  I<-matrix(0,nrow=niters,ncol=NHosts) ;
  I[1,1]<-I0 ;
  Q <- matrix(q1, nrow=NHosts,ncol = NHosts) ;
  for (i in seq(1,NHosts)){
  		Q[i,i] <-q0 ;
  	
  } 
  S <- matrix(N,nrow=1,ncol=NHosts) ;
  S[1,1] <- N - I[1,1] ;
  p <- matrix(0,nrow=1,ncol=NHosts+1) ;
  rf<-.C("sample_SIR",as.integer(NHosts),as.integer(niters),as.double(as.vector(Q)),S<-as.integer(as.vector(S)),I<-as.integer(as.vector(I))) ;
  I<-matrix(rf[[5]],rf[[2]],rf[[1]]) ;
  
  return(list(I=I));
}

reed_frost_multi.3<-function(I0=1, N, NHosts, q0=0.999, q1=0.9999)
{
  # evolution of viral population sizes in multiple hosts

  Q <- matrix(q1, nrow=NHosts,ncol = NHosts) ;
  for (i in seq(1,NHosts)){
  		Q[i,i] <-q0 ;
  } 
  
  p <- matrix(0,nrow=1,ncol=NHosts+1) ;
  I<-.Call("sample_SIR2_R", I0, N, NHosts, Q) ;
  #return(list(I=I,I2=rf[[6]],I2count=rf[[7]]));
  return(I) ;
}

reed_frost_multi_demo<-function(N=1000, NHosts=10)
{
	dyn.load("SIR.so") ;
	rf<-reed_frost_multi.2(N=N,NHosts=NHosts,q0=0.999,q1=0.9999) ;
	return(list(I=rf$I)) ;
}


plot_reed_frost<-function(rf,psfile)
{
  postscript(psfile,horizontal=FALSE)
  iters<- nrow(rf$I) ;
  nHosts<-ncol(rf$I) ;
  par(mfrow=c(6,1)) ;
  for (i in seq(1,nHosts))
  {
    if (!((i-1)%%6))
    {
    	par(mfrow=c(6,1)) ;
    }
    plot(1:iters,rf$I[,i],type="b",col="black",xlab="t",ylab=i) ;
    
  }
  dev.off() 
}

tree_sample<-function(rf,N=1000)
{
	I<-rf[[1]] ;
	I2 <- rf[[2]] ;
	# base 0 to base 1
	I2[,2:3] = I2[,2:3]+1 ;
	I2[,1]=I2[,1]+1 ;
	
	NHosts<-ncol(I) ;
	Hi<-matrix(1,1,NHosts) ;
	
	# most recent time point where there was an infection.
	tnow<-max(which(rowSums(I)>0)) ;
	
	g<-NULL ;
	Nodes<-list() ;
	for (i in seq(1,NHosts))
	{
		H1<-c(1,1+cumsum(rev(I[1:tnow,i]))) ;
		H2<-c(cumsum(rev(I[1:tnow,i]))) ;
		#mx<-which(H2==N)[1] ;
		#mn<-which(H2>0)[1] ;
		mx=tnow ;
		mn = 1;
		Nodes[[i]]<-rbind(rev(H1[mn:mx]),rev(H2[mn:mx])) ;
	}
	
	T=tnow;
	for (tnow in seq(T,2,by=-1))
	{
	H<-which(I[tnow,]>0) ;
	i_now<-which(I2[,1]==(tnow)) ;
	
	I2now<-matrix(I2[i_now,],length(i_now),4) ;
	for (i in seq(1,length(H)))
	{
		# look at previous generation
		It1<-I[tnow-1,H[i]] ;
		
		# get current generation
		mxtmp<-Nodes[[H[i]]] ;
		nc<-seq(mxtmp[1,tnow],mxtmp[2,tnow]);
		
		# these are split depending on which host they come from
		indh<-which(I2now[,3]==H[i]) ;
		Iw<-I[tnow,H[i]] ;
		Ib<-0 ;
		b1<-(H[i]-1)*N ;
		if (length(indh)>0)
		{
			# came from other host
			Ib<-sum(I2now[indh,4]) ;
			Iw<-Iw - Ib;
		}
		
		if (Iw>0)
		{
		    mxtmp<-Nodes[[H[i]]] ;
		    x<-seq(mxtmp[1,tnow-1],mxtmp[2,tnow-1]) ;
		    if (length(x)>1)
		    {
			n1<-sample(x,Iw,replace=TRUE) ;
			}
			else {
				n1<-rep(x,Iw) ;
			}
			#g<-cbind(g,rbind(b1+nc[1:Iw],b1+n1));
			g<-cbind(g,rbind(b1+n1,b1+nc[1:Iw]));
			nc<-setdiff(nc,nc[1:Iw]) ;
		}
		
		if (Ib>0)
		{
			for (j in seq(1,length(indh)))
			{
				mxtmp<-Nodes[[I2now[indh[j],2]]] ;
				x<-seq(mxtmp[1,tnow-1],mxtmp[2,tnow-1]) ;
				if (length(x)>1)
				{
					n2<-sample(x,I2now[indh[j],4],replace=TRUE) ;
				}
				else{
					n2<-rep(x,I2now[indh[j],4]) ;
				}
				b<-(I2now[indh[j],2]-1)*N ;
				#g<-cbind(g,rbind(b1 + nc[1:I2now[indh[j],4]],b+n2));
				g<-cbind(g,rbind(b+n2,b1 + nc[1:I2now[indh[j],4]]));
				nc<-setdiff(nc,nc[1:I2now[indh[j],4]]) ;
			}
		}
		
		
	}
	}
	
	
	tre1<-list() ;
	tre1$edge<-t(g) ;
    tre1$edge.length<-rep(1,ncol(g)) ;
    Nnodes<-length(intersect(unlist(g[1,]),unlist(g[2,]))) ;
	tre1$Nnode<-Nnodes ;
	tips<-setdiff(unlist(g[2,]),unlist(g[1,]))
	tre1$tip.label<-as.character(tips)
	return(tre1);
	
}

getchildren<-function(node,edges)
{
	nbrs<-which(edges[,1]==node) ;
	nbrs_i<-edges[nbrs,2] ;
	return(nbrs_i) ;
}

getparent<-function(node,edges)
{
	nbrs<-which(edges[,2]==node) ;
	p<-edges[nbrs,1] ;
	return(p) ;
}

dfs_tree<-function(node,edges,nl,tl,v)
{
	# if node is discovered and children explored : go up tree
	if (v[2,node]==0)
	{
		node<-getparent(node,edges) ;
	}
	else 
	{
		if (v[2,node]==2)
		{
			# node is undiscovered
			ch<-getchildren(node,edges) ;
			
			# if has more than one child mark as unexplored
			if (length(ch) > 1)
			{		
				v[2,node]<-1 ;
				v[3,node]<-nl ;
				nl<-nl+1 ;
				node<-ch[1];
			}
			else 
			{
				# if no children then a leaf node
				if (length(ch)==0)
				{
					v[2,node]<-0 ;
					v[3,node]<-tl ;
					tl<-tl+1 ;
				}
				else
				{
					v[2,node]<-0 ;
					#v[3,node]<-nl ;
					#nl<-nl+1 ;
					node<-ch[1] ;
				}
			}
		}
		else
		{
			# node is discovered but children not explored
			ch<-getchildren(node,edges) ;
			vch<-v[2,ch] ;
			if (sum(vch) == 0)
			{
				# all children are discovered ...
				v[2,node]<-0 ;
				
			}
			else
			{
				node<-ch[which(vch>0)[1]] ;
			}
			
		}	
	}	
	return(list(node=node,nl=nl,tl=tl,v=v)) ;
}

dfs_bintree<-function(node,edges,nl,tl,v)
{
	# if node is discovered and children explored : go up tree
	if (v[2,node]==0)
	{
		node<-getparent(node,edges) ;
	}
	else 
	{
		if (v[2,node]==2)
		{
			# node is undiscovered
			ch<-getchildren(node,edges) ;
			
			# if has more than one child mark as unexplored
			if (length(ch) > 1)
			{		
				v[2,node]<-1 ;
				v[3,node]<-nl ;
				nl<-nl+1 ;
				node<-ch[1];
			}
			else 
			{
				# if no children then a leaf node
				if (length(ch)==0)
				{
					v[2,node]<-0 ;
					v[3,node]<-tl ;
					tl<-tl+1 ;
				}
				else
				{
					v[2,node]<-0 ;
					#v[3,node]<-nl ;
					#nl<-nl+1 ;
					node<-ch[1] ;
				}
			}
		}
		else
		{
			# node is discovered but children not explored
			ch<-getchildren(node,edges) ;
			vch<-v[2,ch] ;
			if (sum(vch) == 0)
			{
				# all children are discovered ...
				v[2,node]<-0 ;
				
			}
			else
			{
				node<-ch[which(vch>0)[1]] ;
			}
			
		}	
	}	
	return(list(node=node,nl=nl,tl=tl,v=v)) ;
}




tree_convert_to_phylo<-function(tre,N)
{
 
 tips<-as.integer(tre$tip.label) ;
 ntips<-length(tips) ;
 nnode<-tre$Nnode ;
 tl<-1 ;

 # get root
 rt<-setdiff(unlist(tre$edge[,1]),unlist(tre$edge[,2])) ;
 
 nl<-ntips+1 ;
 
 # nodes
 v<-seq(min(tre$edge),max(tre$edge)) ;
 # flags
 v<-rbind(v,rep(2,length(v))) ;
 # new node names
 v<-rbind(v,rep(0,ncol(v))) ;
 
 dtr<-list(node=rt,nl=nl,tl=tl,v=v) ;
 
 while (sum(dtr$v[2,]>0))
 {
 	dtr<-dfs_tree(dtr$node,tre$edge,dtr$nl,dtr$tl,dtr$v) ;
 }
 
 # relabel nodes
 te<-matrix(0,nrow(tre$edge),2) ;
 for (i in seq(min(tre$edge),max(tre$edge)))
 {
 	te[which(tre$edge==i)]<-dtr$v[3,i] ;
 }
 #tre$edge<-te ;
 elen<-tre$edge.length ;
 te2<-tre$edge ;
 
 #m1<-intersect(which(te[,1]>0),which(te[,2]==0)) ;
 
 while (length(which(te==0)>0))
 {
 	m1<-intersect(which(te[,1]>0),which(te[,2]==0)) ;
  	e1<-te2[m1,2] ;
 	m2<-match(e1,te2) ; 
 	te[m1,2]<-te[m2,2] ;
 	te2[m1,2]<-te2[m2,2] ;
 	elen[m1]<-elen[m1]+1 ;
 	te<-te[-m2,] ;
 	elen<-elen[-m2] ;
 	te2<-te2[-m2,] ;
 	nnode<-nnode-length(m2) ;
 	
 }
 nnode<-nnode+1 ;
 # label internal nodes
 m1<-match(seq(ntips+1,ntips+nnode),dtr$v[3,]) ;
 
 tre$node.label <- as.character(ceiling(m1/N))
 tre$edge.length <- elen ;
 tre$edge <- te;
 tre$Nnode <-nnode ; 
  
  
 class(tre)<-"phylo" ;
 return(tre) ;
 }



phylo_convert_to_binphylo<-function(tre)
{
	tree<-reorder(tre) ; 
	edges<-tree$edge ;
	el<-tree$edge.length ;
	nn<-tree$node.label ;
	Nnode<-tree$Nnode ;
	ntips<-length(tree$tip.label) ;
	
	m<-tabulate(tree$edge[,1])
	nl<-max(tree$edge[,1]) ;
	
	
	target<-which(m>2) ;
	for (i in 1:length(target))
	{
		j<-which(edges[,1]==target[i]) ;
		nnew<-length(j)-2 ;
		newnodes<-seq(nl+1,nl+nnew) ;
		
		e<-matrix(0,2*(length(j)-1),2)  ;
		e[c(seq(1,nrow(e),by=2),nrow(e)),2]<-edges[j,2] ;
		e1<-c(edges[j[1],1],newnodes) ;
		e[,1]<-e1[rep(seq(1:(nnew+1)),each=2)] ;
		e[2*(1:nnew),2]<-newnodes ;
		
		newel<-rep(0,nrow(e)) ;
		newel[c(seq(1,nrow(e),by=2),nrow(e))]<-el[j] ;
		
		nn<-c(nn,(rep(nn[target[i]-ntips],nnew))) ;

		edges<-edges[-j,] ;
		el<-el[-j] ;
		Nnode <-Nnode+nnew ;
		edges<-rbind(edges,e) ;
		el<-c(el,newel) ;
		nl<-nl+nnew ;
	}
	
	tree$edge<-edges ;
	tree$edge.length<-el ;
	tree$Nnode <-Nnode ;
	tree$node.label<-nn ;
	
	
	return(tree) ;
}

phylo_convertlabels<-function(phy,N)
{
	tlab<-phy$tip.label ;
	lab<-as.numeric(tlab) ;
	for (i in seq(1, length(lab)))
	{
		j<-which(ceiling(lab/N)==i) ;
		hname<-paste("H",i,sep="") ;
		tlab[j]<-paste(hname,paste("S",1:length(j),sep=""),sep=":")
	}
	phy$tip.label<-tlab ;
	return(phy) ;

}