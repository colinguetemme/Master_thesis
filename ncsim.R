# ncsim --------------------------------------------
# Aim:
#   This function leads to a simulation of a gene tree
#   using the PH dsitribution

# Inputs:
#   'n' the number of individuals 
#   'theta' the population mutation rate

# Outputs:
#   't' a vector with the time spent in each state
#   'Vecstates' a vector with the transient state for this simulation
#   'Vecmut' the simulated SFS

ncsim = function(n, theta = 10, recombrate = 0, T_total = FALSE){
  
  # Initialisation
  times = NULL; T_tot = NULL
  
  k = n # number of distinct lineages
  count = 1 # because k will not goes just from 1 to size 
  Vecstates = matrix(0, nrow=n, ncol=n)
  Vecstates[1, ] = c(n, rep(0, n-1))
  Vecmut = rep(0, n-1) # the number of mutation (in SFS)

  # This loop will continue until reaching the MRCA
  # (until just one lineage)
  while (k > 1){
    lambda = k*(k-1)/2
    
    # Each coalescent time follow a exp dsitribution
    # with a parameter lambda given by the rate matrix
    times[count] = rexp(n=1, rate=lambda)
    
    if(T_total == TRUE){T_tot[count] = times[count]*k} 
           
    # Need this condition because the last coalescent leads to MRCA
    # and is not represented in the state matrix
    # Vector of probability to go to each other states
    
    Vecnext = Vecstates[count,]
    inext = NULL
    # sample two times instead of one (with size 2)
    # to remove 2 of vecnext if the sample choose 
    # two times the same lineage size
    
    # to de the k merger needs to do it as a loop
    
    inext[1] = sample(x = 1:n, size = 1, prob = Vecnext, replace = T) 
    Vecnext[inext[1]]=Vecnext[inext[1]]-1
    inext[2] = sample(x = 1:n, size = 1, prob = Vecnext)
    Vecnext[inext[2]] = Vecnext[inext[2]]-1
    Vecnext[sum(inext)] = Vecnext[sum(inext)]+1
    
    count = count + 1
    Vecstates[count, ] = Vecnext 
    k = k - 1
  }
  # give a vector of mutation rate (corresponding to the SFS)
  mutrate = times*apply(Vecstates[,-n], 2, 'sum')*theta/2
  Vecmut = Vecmut + rpois(n = n-1, lambda = mutrate)
  
  if(T_total == T){times = T_tot}
  return(list(times, Vecstates, Vecmut))
}
