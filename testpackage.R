library(PhaseTypeGenetics)

initDist <- c(1,0,0)

T_Mat <- matrix(c(-6,6,0,
                  0,-3,3,
                  0,0,-1), nrow = 3, ncol = 3, byrow = TRUE)

obj <- contphasetype(initDist, T_Mat)


initDist <- c(1,0,0,0)

P_Mat <- matrix(c(0.4, 0.3, 4/30, 2/30,
                  0, 0.5, 2/9, 1/9,
                  0, 0, 2/3, 0,
                  0, 0, 0, 2/3), nrow = 4, ncol = 4, byrow = TRUE)

S_Total <- discphasetype(initDist, P_Mat)

summary(S_Total)


## Defining the vector of quantiles
t_vec <- seq(0,4, by=0.1)


a = dphasetype(obj, t_vec)

BlockCountProcess(n=10)

initdist = function(n){
  return(c(1,rep(0,n-1)))
}

T_MRCA = function(n){
  mat_MRCA = matrix(0,ncol=n-1,nrow=n-1)
  for (i in 1:(n-2)){
    print(mat_MRCA)
    mat_MRCA[i,i]=-choose(n+1-i,2)
    mat_MRCA[i,i+1]=choose(n+1-i,2)
  }
  mat_MRCA[n-1,n-1]=-1
  return(mat_MRCA)
}
PhaseTypeGenetics:::dphasetype.contphasetype
