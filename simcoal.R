
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
  times = NULL; T_tot = NULL ; k = n ; count = 1
  Matstates = matrix(0, nrow=1, ncol=n)
  Matstates[1, ] = c(n, rep(0, n-1))
  Vecnext=Matstates[1,]
  Vecmut = rep(0, n-1)  
  
  while (k > 1){
    
    lambda = k*(k-1)/2
    times[count] = rexp(n=1, rate=lambda)
    if(T_total == TRUE){T_tot[count] = times[count]*k} 
    inext = NULL
    
    for (i in c(1,2)){
      inext[i] = sample(x = 1:n, size = 1, prob = Vecnext, replace = T) 
      Vecnext[inext[i]]=Vecnext[inext[i]]-1
    }
    
    Vecnext[sum(inext)] = Vecnext[sum(inext)]+1
    count = count + 1
    Matstates = rbind(Matstates, Vecnext, deparse.level = 0)
    k = k - 1
  }
  
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
  
  mutrate = matrix(times*(Matstates[-count,-n])*(theta/2), ncol=n-1)
  mutrate = apply(mutrate, 2, sum)
  Vecmut = rpois(n = n-1, lambda = mutrate)
  if(T_total == T){times = T_tot}
  
  return(list(times, Matstates, Vecmut))
}