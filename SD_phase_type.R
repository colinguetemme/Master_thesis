SD = function(n, i, lambda = 1){
  res <- BlockCountProcess(n)
  
  ## The rate matrix
  Tmat <- res$Rate_Mat
  ## and the corresponding inital distribution
  pi.vec <- c(1,replicate(nrow(Tmat)-1,0))
  
  ## We define an object of type contphasetype
  obj <- contphasetype(pi.vec, Tmat)
  
  ## In order to find the distribution for the site
  ## frequency xi_i, we need a reward vector that
  ## correpsonds to xi_i. Hence
  
  k <- length(which(res$StateSpace_Mat[,1]>=i))
  r.vec <- c(rep(1,k),rep(0,nrow(Tmat)-k))
  
  
  ## In this case, some of the enties in the reward vector
  ## are zero. Therefore, we have to use the function
  ## RewTransDistribution in order to get the transformed
  ## distribution
  newobj <- RewTransDistribution(obj, r.vec)
  
  ## Now we can compute the distribution of the
  ## tail statistic by using the descretization:
  newobj <- discretization(newobj, a=NULL, lambda=lambda)
  return(newobj)
}

n = 5
a = rphasetype(SD(n,1,2), 10000)-2
maxou = max(a)

pl = matrix(0, nrow = n, ncol = maxou)
pl[1,] = hist(a, breaks = 0:maxou)$counts
for (i in 2:n){
  a = rphasetype(SD(n,i,2), 10000)-2
  pl[i,] = hist(a, breaks = 0:maxou)$counts
}

require(plot3D)
hist3D(x = 1:n, y = 1:maxou, z = pl, theta = 240, phi = 45 ,scale = FALSE, expand = 0.01) 

meanpl = t(t(pl) * (1:maxou))
meanpl = apply(meanpl, 1, mean) - mean(pl)
plot(meanpl)
