# Fri Jan 24 14:30:11 2020 ------------------------------


## Simulation of CTMC for population genetics ##
# Need to be implemented in this package 
# also using the block counting process
library(PhaseTypeGenetics)

# simul --------------------------------------------
# Aim:
#   This function leads to a simulation of a gene tree
#   using the PH dsitribution

# Inputs:
#   'n' the number of individuals 
#   'theta' the population mutation rate

# Outputs:
#   't' a vector with the time spent in each state
#   'Vecstates' a vector with the transient state for this simulation

simul = function(n, theta = 10){
  
  # Initialisation
  states = BlockCountProcess(n)$StateSpace_Mat
  Nb_states =  nrow(states)
  rate = BlockCountProcess(n)$Rate_Mat
  times = NULL
  i = 1
  count = 1 # because i will not goes just from 1 to size 
  Vecstates = matrix(0, nrow=n-1, ncol=n-1)
  Vecstates[1, ] = states[1, ]
  Vecmut = rep(0, n-1) # the number of mutation (in SFS)
  
  # This loop will continue until reaching the MRCA
  # because it should be as much coalescent event 
  # as the number of lineages - 1
  while (count < n){
    lambda = rate[i, i]
    
    # Each coalescent time follow a exp dsitribution
    # with a parameter lambda given by the rate matrix
    times[count] = rexp(n=1, rate = -lambda)
    
    # give a vector of mutation rate (corresponding to the SFS)
    mutrate = t[count]*Vecstates[count, ]*theta/2
    Vecmut = Vecmut + rpois(n = n-1, lambda = mutrate)
    
    # Need this condition because the last coalescent leads to MRCA
    # and is not represented in the state matrix
    if (count < n-1){
      # Vector of probability to go to each other states
      VecNext = rate[i, ]/(-rate[i, i])
      VecNext[i] = 0
      inext = sample(x = 1:Nb_states, size = 1, prob = VecNext)
      Vecstates[count+1, ] = states[inext, ] 
      i = inext
    }
    
    count = count + 1
  }
  return(list(times, Vecstates, Vecmut))
}

# Exemples for n=6 of 'simul'

simul(n = 4)
simul(n = 6, theta = 4)

# MakeReward -----------------------------------------
# Aim:
#   Creates the reward vector for the total branch length
#
# Inputs:
#   A matrix 'block' of all the possible states
#   A vector 'vecwhich', only if we want to consider specific
#     | branch states 
#
# Ouput:
#   A vector 'reward' with the reward for each state

MakeReward = function(block, Vecwhich = 1:ncol(block)){
  reward = matrix(block[,Vecwhich], nrow = nrow(block))
  reward = apply(reward, 1, sum)
  return(reward)
}

# Exemple of 'MakeReward'

mat = BlockCountProcess(5)$StateSpace_Mat

MakeReward(block = mat)
MakeReward(block = mat, Vecwhich = 1)  # We do not consider branch with 3 lineage

# theta_est -----------------------------------------
# Aim:
#   Estimates Theta (mutation rate) with the data
#
# Inputs:
#   A vector 'xton', with the value of site frequency we want to use
#   A vector 'vecxton' with the corresponding site frequency spectrum
#
# Ouput:
#   The estimates of Theta

theta_est = function(Vecxton, xton = 1:length(Vecxton)){
  print(Vecxton*xton)
  theta = mean(Vecxton*xton)
  return(theta)
}

# exemple of theta_est:

Vecxton=c(12, 8, 2, 3, 1, 0, 0, 1, 0)

theta_est(Vecxton)
theta_est(Vecxton[3:8], 3:8)


# MatPH_MRCA -----------------------------
# Aim:
#   Creates the contphasetype object corresponding to the MRCA for 
#   a sample of size n
# Input:
#   'n' the sample size
# Ouput:
#   The contphasetype object for the MRCA

MatPH_MRCA = function(n){
  bin = (n:2)*((n:2)-1)/2
  mat = diag(x = -bin, nrow = n-1, ncol = n-1)
  diag(mat[1:(n-2), 2:(n-1)]) = bin[1:(n-2)]
  return(contphasetype(initDist = c(1, rep(0,n-2)), T_Mat = mat))
}

# PH_mutation

PH_mutation = function(n, lambda = 1){
  Vecmut = NULL
  for (i in 1:(n-1)){
    Vecmut[i] = rphasetype(SiteFrequencies(n, 1, i), 1)
  }
  return(Vecmut)
}
