
cSIR_runmodel <- function(outdir, datadir, paramfile)
{
	source(paramfile)
	set.seed(opt$seed)

	# load in data
	load(file=paste(datadir,"/dat.RData",sep="")) ;
	load(file=paste(datadir,"/seqs.RData",sep="")) ;
	if (!is.null(opt$firstN))
	{
  		dat$SN<-dat$SN[opt$firstN] ;
  		dat$ST<-dat$ST[opt$firstN] ;
  		dat$info<-dat$info[opt$firstN] ;
	}

	opt$outdir <- outdir
	
	print(mcp)
	smp <- cSIR_runmcmc(seqs, dat, opt, init, mcp, hp)
	
	#vars<-c("mr","dr","B","ll","tr")
	#for (i in seq(1,length(vars)))
	#{
	#	assign(vars[i],smp[[vars[i]]])
	#	save(list=vars[i],file=paste(outdir,"/",vars[i],".mcmc",sep="")) ;
	#}
	#if (mcp$acc.rate==1)
	#{
  	#	acc_rate <- smp$acc.rate
  	#	save(acc_rate, file=paste(outdir,"/acc_rate.mcmc",sep=""))
	#}  
	warnings()
}

cSIR_runmcmc <- function(x, dat, opt, init, mcmc.params, hp.params)
{
  # Runs mcmc chain targeting posterior distribution
  # Args:
  #   x: sequences (in DNAbin format)
  #	  dat:	structure holding information about the data:
  #		SN: number of sequences for each host
  #		ST: time of sampling
  #		opt:	options for run
  #		init:	initial conditions
  #		mcmc.params: 
  #		hp: hyperparameters 
  #
  # Returns:
  #   ch: list of sample chains
  
  SN<-dat$SN 
  ST<-dat$ST 
  dat.params<-list(I0=init$I0, NS=init$NS, NHosts=length(SN))
  dat$STin = dat$ST - min(dat$ST)
  
  vars <- c("mr",		# 1: Mutation rate + tr
  			"dr", 		# 2: Recovery rate + tr
  			"T", 		# 3: tr
  			"B", 		# 4: Transmission matrix + tr
  			"bn", 		# 5: Bottleneck parameter + tr
  			"t_off")	# 6: time offset + tr
  
  # initialise model
  params <- cSIR_modelinit(x, init, dat, dat.params, mcmc.params)
   # initialise priors
  priors <- list()
  priors$mr <- cSIR_mrprior(mr=params$mr, hp.params$mr, mode="d", is.log=TRUE)
  priors$dr <- cSIR_drprior(dr=params$dr, hp.params$dr, mode="d", is.log=TRUE)
  priors$B <- cSIR_Bprior(B=params$B, dat.params$NHosts, hp.params$B.br1, hp.params$B.br2, mode="d", is.log=TRUE)
  priors$bn <- cSIR_bnprior(bn=params$bn, hp.params$bn, mode="d", is.log=TRUE)        
  priors$t_off <- cSIR_t_offprior(t_off=params$t_off, hp.params$t_off, mode="d", is.log=TRUE)
  params$priors <- priors
  
  # initialise possible moves
  moves <- cSIR_movesinit(opt$movevars, vars) 
  
  prob <- rep(1/length(moves), length(moves))

  # initialise structure for saving chain		
  ch <- cSIR_chainsinit(opt$chainvars)
  ch <- cSIR_chainsupdate(ch, opt$chainvars, params)
  
  # convert data 
  x <- phyDat(x) 
  
  Ninstances = 1 
  #Ninstances = mcmc.params$Np 
  rv <-  .Call("smc_init_R", Ninstances, x, attr(x, "weight"), sum(dat$SN), attr(x, "nr"))
 
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
	while (t <= mcmc.params$Niters)
	{
		# make new proposal
		m <- cSIR_choosemove(moves, prob)

   
		switch(as.character(m),
		  "1"={
		    params <- cSIR_mrupdate(x, params, dat.params, dat, hp.params, mcmc.params)
		  },
		  "2"={
		    params <- cSIR_drupdate(x, params, dat.params, dat, hp.params, mcmc.params)
		  },
		  "3"={
		    params <- cSIR_Tupdate(x, params, dat.params, dat, hp.params, mcmc.params)
		  },
		  "4"={ 
		    params <- cSIR_Bupdate(x, params, dat.params, dat, hp.params, mcmc.params)
    	  },
		  "5"={ 
		    params <- cSIR_bnupdate(x, params, dat.params, dat, hp.params, mcmc.params)
		  },
		  "6"={ 
		    params <- cSIR_t_offupdate(x, params, dat.params, dat, hp.params, mcmc.params)
		  }
		  )
		ch <- cSIR_chainsupdate(ch, opt$chainvars, params)
    	if (mcmc.params$acc.rate==1)
    	{
      		# running estimate of acceptance rate
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
    
    	# print status
    	if (t %% 1 == 0)
    	{
    		print(sprintf("%s %i loglikelihood   %8.4f   Mutation rate   %8.8f  Death rate %8.4f BN %8.4f t_off %8.4f %i Move %i", format(Sys.time(), "%X"), t, params$ll, params$mr, params$dr,params$bn, params$t_off, params$K, m))
			#print(params$B)
			#print(t(params$tr_list$Itraj))
		}
	
    
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
		
	.Call("smc_free_instances_R", Ninstances) ;
	return(ch)
}

# Initialisation functions
cSIR_modelinit <- function(x, init, dat, dat.params, mcmc.params)
{
  # Initialises model with initialisation file
  #
  # Args:
  #   init: list of initial values
  #   dat: data info
  #   x: sequence data
  #	  dat.params
  #   mcmc.params
  #
  # Returns:
  #   params
  
  
  dat$ST <- dat$STin + init$t_off
  ST <- dat$ST
  
  SN <- dat$SN
  Bcon <- cSIR_Bstruct(ST=ST) 
  Bcon <- init$Bcon
  if (!is.matrix(init$B))
  {
  	B <- cSIR_Binit(init, Bcon)
  }
  else
  {
  	B <- init$B
  }
  dr <- cSIR_drinit(init) 
  mr <- cSIR_mrinit(init) 
  trinit <- cSIR_trinit(init, x, dat)
  
  llcur <- log(0)  
  params <- list(B=B, dr=dr, mr= mr, ll=llcur, Bcon=Bcon, tr_list=list(lo=trinit$lo, s=trinit$s), bn=init$bn, t_off=init$t_off, is.acc=0, K=init$K)
  
  a <- .Call("smc_draw_prior_R", 1, dat.params$I0, dat.params$NS, dat.params$NHosts, params$B, params$dr, dat$ST, dat$SN, params$bn, params$K) 
  #tr_list <- cSIR_drawtr_list(params, dat.params, dat, mcmc.params$tr)
  tr_list <- list(tr = a[[1]])
  sir <- list(I=a[[2]], S=a[[3]], T=a[[4]], Iend=a[[5]])
  sir$T[,2:3] <- sir$T[,2:3] + 1 ;
  tr_list$Itraj <- cSIR_Itraj(sir, Tmax=10) 
  params$tr_list <- tr_list
  params$is.acc <- 0
  params$ll <- llcur
  
  
  return(params)
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
 # B <- matrix(0, nrow=NHosts, ncol=NHosts)
  for (i in seq(1, NHosts))
  {
    B[i,i] <- init$B$br1
  }
  #B[1,2] <- init$B$br1
  #B[2,3] <- init$B$br2
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

# Move functions

cSIR_chainsupdate <- function(ch, chainvars, params)
{
  n <- length(ch[[1]])
  for (i in seq(1,length(chainvars)))
  {
  	  switch(chainvars[i],
  	  	"tr"={
      	  	ch$tr[[n+1]] <- params$tr_list$tr 
      	},
      	"Itraj"={
      		ch$Itraj[[n+1]] <- params$tr_list$Itraj
      	},
      	{
        	ch[[chainvars[i]]][[n+1]] <- params[[chainvars[i]]]
      	}
      	)
  }
  return(ch)
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
  
# updating variables  
  
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
  #      params
  #
  
	#tr_list <- cSIR_drawtr_list_smc2(x, params, dat.params, dat, mcmc.params)
	tr_list <- cSIR_trproposal(x, params, dat.params, dat, hp.params, mcmc.params)
	params$is.acc <- 0
	priorsum <- sum(unlist(params$priors))
  	if (tr_list$is.acc==1)
	{
	  	llcur<-params$ll
	  	print(tr_list$ll)
      	pacc <- min(tr_list$ll + priorsum - llcur, 0)	
      	
	  	if (log(runif(1)) <= pacc)
	  	{
	    	params$tr_list <- tr_list
	    	params$ll <- tr_list$ll + priorsum
	    	params$is.acc <- 1
	  	}  
	}
  
	return(params)
}

cSIR_Bupdate <- function(x, params, dat.params, dat, hp.params, mcmc.params)
{
  # Updates matrix of tranmission rates B
  #
  # Args:
  #   x: sequence data
  #   params: current model parameters
  #   dat.params: other parameters
  #   dat: data info
  #   hp.params: parameters for the priors
  #   mcmc.params:  proposal parameters for B
  # Returns:
  #   params
  #

 
	oldB <- params$B
	newB <- cSIR_Bproposal(params$B, params$Bcon, mcmc.params$B)  
	params$B <- newB
	tr_list <- cSIR_trproposal(x, params, dat.params, dat, hp.params, mcmc.params)
	params$is.acc <- 0
	params$B <- oldB
  	if (tr_list$is.acc==1)
	{
	  	llcur<-params$ll
	  	print(tr_list$ll)
	  	oldprior <- params$priors$B
	  	newprior <- cSIR_Bprior(newB, dat.params$NHosts, hp.params$B.br1, hp.params$B.br2, mode="d", is.log=TRUE)
	  	newpriors <- newprior - params$priors$B + sum(unlist(params$priors))
	  	
      	pacc <- min(tr_list$ll + newpriors - llcur, 0)	
      	
	  	if (log(runif(1)) <= pacc)
	  	{
	    	params$tr_list <- tr_list
	    	params$priors$B <- newprior 
	    	params$ll <- tr_list$ll + sum(unlist(params$priors))
	    	params$is.acc <- 1
	    	params$B <- newB
	  	}  
	}
	
	
	return(params)	
}

cSIR_bnupdate <- function(x, params, dat.params, dat, hp.params, mcmc.params)
{
  # Updates bottleneck parameter bn
  #
  # Args:
  #   x: sequence data
  #   params: current model parameters
  #   dat.params: other parameters
  #   dat: data info
  #   hp.params: parameters for the priors
  #   mcmc.params:  proposal parameters for bn
  # Returns:
  #   params
  #

 	oldbn <- params$bn
  	oldK <- params$K
  	newbn <- cSIR_bnproposal(params$bn, mcmc.params$bn)  
  	newK <- rbinom(1, dat.params$NS, params$bn)
  	params$bn <- newbn
  	params$K <- newK 
  	tr_list <- cSIR_trproposal(x, params, dat.params, dat, hp.params, mcmc.params)
	params$bn <- oldbn
  	params$K <- oldK
 	if (tr_list$is.acc==1)
	{
	  	llcur<-params$ll
	  	print(tr_list$ll)
	  	oldprior <- params$priors$bn
	  	newprior <- cSIR_bnprior(newbn, hp.params$bn, mode="d", is.log=TRUE)
	  	newpriors <- newprior - params$priors$bn + sum(unlist(params$priors))
	  	
      	pacc <- min(tr_list$ll + newpriors - llcur, 0)	
      	
	  	if (log(runif(1)) <= pacc)
	  	{
	    	params$tr_list <- tr_list
	    	params$priors$bn <- newprior
	    	params$ll <- tr_list$ll + sum(unlist(params$priors))
	    	params$is.acc <- 1
	    	params$bn <- newbn
	    	params$K <- newK
	  	}  
	}
	
	
	return(params)	
}

cSIR_drupdate <- function(x, params, dat.params, dat, hp.params, mcmc.params)
{
  # Updates removal rate dr
  #
  # Args:
  #   x: sequence data
  #   params: current model parameters
  #   dat.params: other parameters
  #   dat: data info
  #   hp.params: parameters for the priors
  #   mcmc.params:  proposal parameters for dr
  # Returns:
  #   params
  #

 
	olddr <- params$dr
	newdr <- cSIR_drproposal(params$dr, mcmc.params$dr) 
    params$dr <- newdr
	tr_list <- cSIR_trproposal(x, params, dat.params, dat, hp.params, mcmc.params)
	params$is.acc <- 0
	params$dr <- olddr
  	if (tr_list$is.acc==1)
	{
	  	llcur<-params$ll
	  	print(tr_list$ll)
	  	oldprior <- params$priors$dr
	  	newprior <- cSIR_drprior(newdr, hp.params$dr, mode="d", is.log=TRUE)
	  	newpriors <- newprior - params$priors$dr + sum(unlist(params$priors))
	  	
      	pacc <- min(tr_list$ll + newpriors - llcur, 0)	
      	
	  	if (log(runif(1)) <= pacc)
	  	{
	    	params$tr_list <- tr_list
	    	params$priors$dr <- newprior
	    	params$ll <- tr_list$ll + sum(unlist(params$priors))
	    	params$is.acc <- 1
	    	params$dr <- newdr
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
	newmr <- cSIR_mrproposal(params$mr, mcmc.params$mr) 
    params$mr <- newmr
	tr_list <- cSIR_trproposal(x, params, dat.params, dat, hp.params, mcmc.params)
	params$is.acc <- 0
	params$mr <- oldmr
  	if (tr_list$is.acc==1)
	{
	  	llcur<-params$ll
	  	print(tr_list$ll)
	  	
	  	oldprior <- params$priors$mr
	  	newprior <- cSIR_mrprior(newmr, hp.params$mr, mode="d", is.log=TRUE)
	  	newpriors <- newprior - params$priors$mr + sum(unlist(params$priors))
	  	pacc <- min(tr_list$ll + newpriors - llcur , 0)	
      	
	  	if (log(runif(1)) <= pacc)
	  	{
	    	params$tr_list <- tr_list
	    	params$priors$mr <- newprior
	    	params$ll <- tr_list$ll + sum(unlist(params$priors))
	    	params$is.acc <- 1
	    	params$mr <- newmr
	  	}  
	}
	
	
	return(params)	
}

cSIR_t_offupdate <- function(x, params, dat.params, dat, hp.params, mcmc.params)
{
  # Updates offset t_off
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

 
	
	oldt_off <- params$t_off
	newt_off <- cSIR_t_offproposal(params$t_off, mcmc.params$t_off) 
    params$t_off <- newt_off 
	tr_list <- cSIR_trproposal(x, params, dat.params, dat, hp.params, mcmc.params)
	params$t_off <- oldt_off
  	if (tr_list$is.acc==1)
	{
	  	llcur<-params$ll
	  	print(tr_list$ll)
	  	oldprior <- params$priors$mr
	  	newprior <- cSIR_t_offprior(newt_off, hp.params$t_off, mode="d", is.log=TRUE)
	  	newpriors <- newprior - params$priors$t_off + sum(unlist(params$priors))
	  	pacc <- min(tr_list$ll + newpriors - llcur , 0)	
      
	  	if (log(runif(1)) <= pacc)
	  	{
	    	params$tr_list <- tr_list
	    	params$priors$t_off <- newprior
	    	params$ll <- tr_list$ll + sum(unlist(params$priors))
	    	params$is.acc <- 1
	    	params$t_off <- newt_off
	  	}  
	}
	
	
	return(params)	
}

# Proposals
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

cSIR_trproposal <- function(x, params, dat.params, dat, hp.params, mcmc.params)
{
	dat$ST <- dat$STin  + params$t_off
	Np <- mcmc.params$Np		# number of particles
	Tmax = mcmc.params$Tmax
	#Np <- 100
	
	w <- attr(x, "weight")		# weight for each site
	a <- .Call("smc_draw_R", Np, dat.params$I0, dat.params$NS, dat.params$NHosts, params$B, params$dr, dat$ST, dat$SN, params$bn, params$K, x, sum(dat$SN), attr(x, "nr"), w, params$mr)
	tr_list <- list(tr = a[[1]], ll=a[[6]])
  	sir <- list(I=a[[2]], S=a[[3]], T=a[[4]], Iend=a[[5]])
  	sir$T[,2:3] <- sir$T[,2:3] + 1 
  	tr_list$Itraj <- cSIR_Itraj(sir, Tmax=Tmax)
	tr_list$is.acc <- 1		# ?
	return(tr_list)
}

# Priors

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


# other

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

sample_cSIR_C <-function(I0, NS, NHosts, B, dr)
{
	sir1<-.Call("sample_cSIR_R",I0, NS, NHosts, B, dr) ;
	sir<-list(I=sir1[[1]],S=sir1[[2]],T=sir1[[3]]) ;
	sir$T[,2:3]<-sir$T[,2:3]+1 ;
	return(sir) ;
}

sample_cSIR_S_C <- function(I0, NS, NHosts, B, dr, SN, ST, bnprob, K)
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
  
	sir1 <- .Call("sample_cSIR_S_R", I0, NS, NHosts, B, dr, ST, SN, bnprob, K) ;
	sir <- list(I=sir1[[1]], S=sir1[[2]], T=sir1[[3]], Iend=sir1[[4]]) ;
	sir$T[,2:3] <- sir$T[,2:3] + 1 ;
	sir$B <- B ;
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
  	tr_list <- .Call("tree_reconstruct_with_partialll2", sir, dat.params$NHosts, dat, params$B, dat.params$NS, params$bn, rep(0,dat.params$NHosts), dat$SN, 0) 
	names(tr_list) <- c("tr", "bt", "ll", "mig")	
	tr_list$tr$tip.label <- NULL
	for (i in seq(1,dat.params$NHosts))
	{
		tr_list$tr$tip.label<-c(tr_list$tr$tip.label,dat$info[[i]]$key) ;		
	}		
    #tr_list <- cSIR_trdraw2(sir, dat, dat.params, params)
	#ll <- tr_list$ll
    #tr_list <- cSIR_trdraw3(tr_list$T, dat, params$tr_list)
	#tr_list$Itraj <- cSIR_Itraj(sir)
	#tr_list$ll <- ll
  }

  tr_list$is.acc <- sir$is.acc
  return(tr_list)
}

cSIR_obs <- function(SN, tot)
{
	L <- length(SN)
	o <- vector(mode="numeric", length=L)
	i <- 1 
	while ((i <= L) & (tot > 0))
	{
		o[i] <- ifelse(tot > SN[i], SN[i], tot)
		tot <- tot - SN[i]
		i <- i + 1
	}
	return(o)
}

cSIR_obs_labs <- function(SN, tot, nodeinfo)
{
	L <- length(SN)
	labs <- NULL
	i <- 1 
	while ((i <= L) & (tot > 0))
	{
		if (tot > SN[i])
		{
			labs <- c(labs, nodeinfo[[i]]$key)
		}
		else
		{
			labs <- c(labs, nodeinfo[[i]]$key[1:tot])
		}
		tot <- tot - SN[i]
		i <- i + 1
	}
	return(labs)
}

cSIR_Itraj <- function(sir, N=100, Tmax=15)
{
  # Find average infected population sizes over time for each host
  # Args:
  #   sir: list of:
  #     I:  infected size trajectory
  #     S:  susceptible size trajectory
  #     T:  matrix of events: T[i,]={time, from host, to host, [1=birth,-1=death]}
  #     Iend: vector of infected size for each host at times ST  
  #		N
  #		Tmax
  #     is.acc: 1=new sir sample accepted, 0= not accepted
  # Returns:
  #     
  # get most recent time point
  
  T_end <- sir$T[nrow(sir$T),1]
  
  # transform time points to be relative to present
  T <- sir$T[,1] - T_end
  
  # ti is the intervals
  ti <- seq(-Tmax, 0, length.out=N+1)
  Itraj <- matrix(0,nrow=length(ti),ncol=ncol(sir$I))
  for (i in seq(1,ncol(sir$I)))
  {
  	ap <- approx(T, sir$I[,i], ti, method="constant")
  	Itraj[,i] <- ap$y
  }
  return(Itraj)
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

cSIR_savesmps <- function(opt, smp)
{
  vars<-c("mr","dr","B","ll","tr")
  for (i in seq(1,length(vars)))
  {
    assign(vars[i],smp[[vars[i]]])
    save(list=vars[i],file=paste(outdir,"/",vars[i],".mcmc",sep="")) ;
  }  
}

sumofexp <- function(loga, sort=0)
{
	## returns log(sum(exp(loga)))
	if (sort==1)
		loga<-sort(loga)
	
	L <- loga[1]	
	if (length(loga) > 1)
	{
	for (i in seq(2, length(loga)))
	{
		Lv <- c(L,loga[i])
		mll <- max(Lv)
		L <- log(sum(exp(Lv - mll))) + mll 
	}
	}
	#L <- L- log(length(loga))
	return(L)

}






