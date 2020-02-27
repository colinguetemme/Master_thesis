ncsim = function(n, T_total = FALSE){
  
  if (n<=0 || length(n)!=1 || round(n)!=n){
    stop('n should be a positive integer')
  }
  
  if (!is.logical(T_total)){
    stop('T_total should be a logical value')
  }
  
  # Initialisation
  times = NULL; T_tot = NULL
  
  k = n  # number of distinct lineages
  count = 1  # because k will not goes just from 1 to size 
  Matstates = matrix(0, nrow=1, ncol=n)
  Matstates[1, ] = c(n, rep(0, n-1))
  Vecnext=Matstates[1,]
  Vecmut = rep(0, n-1)  # the number of mutation (in SFS)
  
  # This loop will continue until reaching the MRCA
  # (until just one lineage)
  while (k > 1){
    
    # plot(Vecrate)
    # Each coalescent time follow a exp dsitribution
    # with a parameter lambda given by the rate matrix
    
    lambda = k*(k-1)/2
    times[count] = rexp(n=1, rate=lambda)
    
    if(T_total == TRUE){T_tot[count] = times[count]*k} 
    
    # Need this condition because the last coalescent leads to MRCA
    # and is not represented in the state matrix
    # Vector of probability to go to each other states
    
    inext = NULL
    # sample two times instead of one (with size 2)
    # to remove 2 of vecnext if the sample choose 
    # two times the same lineage size
    for (i in c(1,2)){
      inext[i] = sample(x = 1:n, size = 1, prob = Vecnext, replace = T) 
      Vecnext[inext[i]]=Vecnext[inext[i]]-1
    }
    Vecnext[sum(inext)] = Vecnext[sum(inext)]+1
    
    count = count + 1
    Matstates = rbind(Matstates, Vecnext, deparse.level = 0)
    k = k - 1
  }
  # give a vector of mutation rate (corresponding to the SFS)
  
  
  
  if(T_total == T){times = T_tot}
  return(list(times, Matstates))
}

mutsim = function(simtree, theta = 10){
  
  mutrate = matrix(simtree[[1]]*(simtree[[2]][-nrow(simtree[[2]]),-ncol(simtree[[2]])])*(theta/2), ncol=ncol(simtree[[2]])-1)
  mutrate = apply(mutrate, 2, sum)

  Vecmut = rpois(n = ncol(simtree[[2]])-1, lambda = mutrate)
  return (Vecmut)
}
  
LambdaB = function(m, k, alpha=1.9999){
  lmk=beta(k-alpha, m-k+alpha)/beta(alpha, 2-alpha)
  return(lmk)
}




bcsim = function(n, alpha = 1.9999, T_total = FALSE){
  
  if (n<=0 || length(n)!=1 || round
      
      (n)!=n){
    stop('n should be a positive integer')
  }
  
  if (!is.logical(T_total)){
    stop('T_total should be a logical value')
  }
  
  if (alpha<1 || alpha >= 2 || length(alpha)!= 1){
    stop('alpha should be a positive numeric between [1, 2)')
  }
  
  times = NULL; T_tot = NULL; k = n; count = 1 
  Matstates = matrix(0, nrow=1, ncol=n)
  Matstates[1, ] = Vecnext = c(n, rep(0, n-1))
  Vecmut = rep(0, n-1)  
  
  while (k > 1){
    Vecrate = rep(0,n)
    for (i in 2:k){Vecrate[i]=(LambdaB(k, i, alpha = alpha))*choose(k, i)}
    lambda=sum(Vecrate)
    times[count] = rexp(n=1, rate=lambda)
    kmerger=sample(1:n, 1, prob = Vecrate/lambda)
    
    
    if(T_total == TRUE){T_tot[count] = times[count]*k} 
    
    inext = NULL
    for (i in 1:kmerger){
      inext[i] = sample(x = 1:n, size = 1, prob = Vecnext, replace = T) 
      Vecnext[inext[i]]=Vecnext[inext[i]]-1
    }
    Vecnext[sum(inext)] = Vecnext[sum(inext)]+1
    count = count + 1
    Matstates = rbind(Matstates, Vecnext, deparse.level = 0)
    k = k - (kmerger-1)
  }
  
  if(T_total == T){times = T_tot}
  return(list(times, Matstates))
}

nvbclass = function(alpha, n=100, replicate=1000){
  nsimulationSFS = matrix(0,ncol=n-1, nrow=replicate)
  bsimulationSFS = matrix(0,ncol=n-1, nrow=replicate)
  
  for (i in 1:replicate){
    simul = ncsim(n = n)
    nsimulationSFS[i,] = mutsim(simul)
  }
  for (i in 1:replicate){
    simul = bcsim(n = n, alpha = alpha)
    bsimulationSFS[i,] = mutsim(simul)
  }
  
  xn = nsimulationSFS[,1]/(apply(nsimulationSFS,1,sum)+0.00001)
  #yn = apply(nsimulationSFS[,100:(n-1)],1,sum)/apply(nsimulationSFS,1,sum)
  yn = apply(nsimulationSFS,1,sum)
  
  xb = bsimulationSFS[,1]/(apply(bsimulationSFS,1,sum)+0.00001)
  #yb = apply(bsimulationSFS[,100:(n-1)],1,sum)/apply(bsimulationSFS,1,sum)
  yb = apply(bsimulationSFS,1,sum)
  
  #plot(xn,yn, xlim=c(0,1),ylim=c(0,200))
  #points(xb,yb, col='red')
  
  #hist(xn, breaks = 30)
  #points(x = seq(0,0.6,length.out = 100),dgamma(seq(0,10,length.out = 100), 3.5, 1)*350)
  return(c(xn,xb,yn,yb))
}
