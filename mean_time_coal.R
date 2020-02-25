
ncsim = function(n, theta = 10, model = 'n', T_total = FALSE){
  
  if (n<=0 || length(n)!=1 || round(n)!=n){
    stop('n should be a positive integer')
  }
  
  if (theta<0 || length(theta)!= 1){
    stop('theta should be a positive numeric')
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
    times[count] = 1/lambda
    
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
  
  mutrate = matrix(times*(Matstates[-count,-n])*(theta/2), ncol=n-1)
  mutrate = apply(mutrate, 2, sum)
  
  
  Vecmut = rpois(n = n-1, lambda = mutrate)
  
  if(T_total == T){times = T_tot}
  return(list(times, Matstates, Vecmut))
}




LambdaB = function(m, k, alpha=1.9999){
  lmk=beta(k-alpha, m-k+alpha)/beta(alpha, 2-alpha)
  return(lmk)
}


bcsim = function(n, alpha = 1.9999, theta = 10, T_total = FALSE){
  
  if (n<=0 || length(n)!=1 || round(n)!=n){
    stop('n should be a positive integer')
  }
  
  if (theta<0 || length(theta)!= 1){
    stop('theta should be a positive numeric')
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
    times[count] = 1/lambda
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
  
  mutrate = matrix(times*(Matstates[-count,-n])*(theta/2), ncol=n-1)
  mutrate = apply(mutrate, 2, sum)
  
  Vecmut = rpois(n = n-1, lambda = mutrate)
  
  if(T_total == T){times = T_tot}
  return(list(times, Matstates, Vecmut))
}
