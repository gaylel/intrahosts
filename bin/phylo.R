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
      if (length(k1)==0|length(k2)==0)
      {
        print(t$edge[k,])
        print(t$edge.length[k])
      }
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
  ntips <- length(tr$tip.label)
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

phylo_ne_2<-function(tr,h=-1)
{

	# swap any 2 subtrees
	
	# list of nodes
	ntips <- length(tr$tip.label)
	Nodes <- seq(1 , (2*ntips - 1))

	# choose i
	i <- sample(Nodes, 1)
	ip <- tr$edge[tr$edge[,2]==i, 1]
	
	# check that i is valid
	rootno <- ntips + 1
	rootch <- tr$edge[tr$edge[,1]==rootno, 2]
	if  (!(i %in% c(rootno,rootch)))
	{
		# pick j at uniform over tree
		j <- sample(Nodes[-i], 1)
		
		d <- dist.nodes(tr)
		if ((d[ip,rootno] <= d[j,rootno]) & (d[j,rootno] <= d[i, rootno]))
		{
			if (h[i] == h[j])
			{
				tr <- ne_swapnodes_2(tr,j,i)
			}
		}
		
	}
	
	return(tr)
}

phylo_ne_3<-function(tr,h=-1)
{

	# swap any 2 subtrees that have aunt niece relationship
	
	# list of nodes
	ntips <- length(tr$tip.label)
	Nodes <- seq(1 , (2*ntips - 1))

	# choose i
	i <- sample(Nodes, 1)
	ip <- tr$edge[tr$edge[,2]==i, 1]
	
	# check that i is valid
	rootno <- ntips + 1
	rootch <- tr$edge[tr$edge[,1]==rootno, 2]
	if  (!(i %in% c(rootno,rootch)))
	{
		# choose j as aunt
		jp <- tr$edge[tr$edge[, 2]==ip, 1]
		j <- tr$edge[tr$edge[,1]==jp, 2]
		j <- setdiff(j, ip)
		d <- dist.nodes(tr)
		if ((d[ip,rootno] <= d[j,rootno]) & (d[j,rootno] <= d[i, rootno]))
		{
			if (h[i] == h[j])
			{
				tr <- ne_swapnodes_2(tr,j,i)
			}
		}
		
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
  ntips <- length(tr$tip.label)
	ia<-which(tr$edge[,2]==a)
	ib<-which(tr$edge[,2]==b)
	ic<-which(tr$edge[,2]==tr$edge[ib,1])
	if ((b > ntips & tr$edge[ia,1] < b) | b <= ntips)
	{
    if ((a > ntips & tr$edge[ib,1] < a) | a <= ntips)
    {
	    tr$edge[ia,2]<-b
	    tr$edge[ib,2]<-a
	    ea<-tr$edge.length[ia]
	    eb<-tr$edge.length[ib]
	    ec<-tr$edge.length[ic]
	    tr$edge.length[ib]<-ea-ec
	    tr$edge.length[ia]<-eb+ec
    }
    else
    {
      print("skipping swap")
    }
	}
  return(tr)
	
}



ne_swapnodes_2<-function(tr,j,i)
{
	# swap nodes j ('aunt') and i ('niece')
  	ntips <- length(tr$tip.label)
	ij<-which(tr$edge[,2]==j)
	ii<-which(tr$edge[,2]==i)
	
	d <- dist.nodes(tr)
	rootnode <- length(tr$tip.label) + 1
	
	if (i >= rootnode & tr$edge[ij,1] < i | i<rootnode)
	{
		if (j >= rootnode & tr$edge[ii,1] < j | j<rootnode)
		{
		tr$edge[ij,2]<-i
		tr$edge[ii,2]<-j
		ei <- d[rootnode,j] - d[rootnode,tr$edge[ii, 1]]
		ej <-  d[rootnode, i] - d[rootnode, tr$edge[ij, 1]]
		tr$edge.length[ij]<-ej
		tr$edge.length[ii]<-ei
    	}
    }
  return(tr)
	
}
ne_filterauntniecebylength<-function(tr,fam)
{
	# only retain valid swaps i.e. parent-aunt has to be longer than parent-sister
	# parent aunt
	e1<-tr$edge.length[match(fam[4,],tr$edge[,2])]
	# parent child
	e2<-tr$edge.length[match(fam[2,],tr$edge[,2])]
	fam<-fam[,e1>=e2]
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

pml.gl <- function (tree, data, bf = NULL, Q = NULL, inv = 0, k = 1, shape = 1, 
                    rate = 1, model = NULL, ...) 
{
  call <- match.call()
  extras <- match.call(expand.dots = FALSE)$...
  pmla <- c("wMix", "llMix")
  existing <- match(pmla, names(extras))
  wMix <- ifelse(is.na(existing[1]), 0, eval(extras[[existing[1]]], 
                                             parent.frame()))
  llMix <- ifelse(is.na(existing[2]), 0, eval(extras[[existing[2]]], 
                                              parent.frame()))
  if (class(tree) != "phylo") 
    stop("tree must be of class phylo")
  if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
    tree <- reorder(tree, "postorder")
  if (any(tree$edge.length < 0)) {
    tree$edge.length[tree$edge.length < 0] <- 1e-08
    warning("negative edges length changed to 0!")
  }
  if (class(data)[1] != "phyDat") 
    stop("data must be of class phyDat")
  if (is.null(tree$edge.length)) 
    stop("tree must have edge weights")
  if (any(is.na(match(tree$tip, attr(data, "names"))))) 
    stop("tip labels are not in data")
  data <- subset(data, tree$tip.label)
  levels <- attr(data, "levels")
  weight <- attr(data, "weight")
  nr <- attr(data, "nr")
  type <- attr(data, "type")
  if (type == "AA" & !is.null(model)) {
    model <- match.arg(model, .aamodels)
    getModelAA(model, bf = is.null(bf), Q = is.null(Q))
  }
  if (type == "CODON") 
    Q <- as.numeric(.syn > 0)
  if (is.null(bf)) 
    bf <- rep(1/length(levels), length(levels))
  if (is.null(Q)) 
    Q <- rep(1, length(levels) * (length(levels) - 1)/2)
  m <- 1
  eig <- edQt(bf = bf, Q = Q)
  w <- rep(1/k, k)
  if (inv > 0) 
    w <- (1 - inv) * w
  if (wMix > 0) 
    w <- wMix * w
  g <- discrete.gamma(shape, k)
  if (inv > 0) 
    g <- g/(1 - inv)
  g <- rate * g
  INV <- Matrix(lli(data, tree), sparse = TRUE)
  ll.0 <- as.matrix(INV %*% (bf * inv))
  if (wMix > 0) 
    ll.0 <- ll.0 + llMix
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  nTips <- as.integer(length(tree$tip.label))
  on.exit(.C("ll_free"))
  .C("ll_init", nr, nTips, nc, as.integer(k))
  tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q, 
                 levels = attr(data, "levels"), inv = inv, rate = rate, 
                 g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, 
                 wMix = wMix, site = TRUE)
  if (type == "CODON") {
    df <- length(tree$edge.length) + (k > 1) + (inv > 0) + 
      length(unique(bf)) - 1
  }
  else df = length(tree$edge.length) + (k > 1) + (inv > 0) + 
    length(unique(bf)) - 1 + length(unique(Q)) - 1
  result = list(logLik = tmp$loglik, inv = inv, k = k, shape = shape, 
                Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik, weight = weight, 
                g = g, w = w, eig = eig, data = data, model = model, 
                INV = INV, ll.0 = ll.0, tree = tree, lv = tmp$resll, 
                call = call, df = df, wMix = wMix, llMix = llMix)
  if (type == "CODON") {
    result$dnds <- 1
    result$tstv <- 1
  }
  class(result) = "pml"
  result
}