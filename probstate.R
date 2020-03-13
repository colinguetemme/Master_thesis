probstate = function(state){
  nstate = nrow(state)
  statevec = c(1, rep(0,(nstate-1)))
  for (k in 2:nstate){
    for (i in 1:(k-1)){
      statevec[k] = statevec[k] + (statevec[i]*(state[i,k]/(-state[i,i])))
    }
  }
  return(statevec)
}
