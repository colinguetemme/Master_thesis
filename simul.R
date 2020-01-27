# Fri Jan 24 14:30:11 2020 ------------------------------


## Simulation of CMC for population genetics ##
# Need to be implemented in this package 
# also using the block counting process
library(PhaseTypeGenetics)

# simul --------------------------------------------
# aim:
#   This function leads to a simulation of a gene tree
#   using the PH dsitribution

# Input :
#   Q the transition matrix 
#   a the initial distribution vector

# Output:
#   t a vector with the time spent in each state
#   Vecstates a vector with the transient state for this simulation

simul = function(Q,a){
  # Initialisation
  
  t = NULL
  size = ncol(Q)
  states = BlockCountProcess(size)$StateSpace_Mat
  i = 1
  count = 1 # because i will not goes just from 1 to size 
  Vecstates = matrix(0,nrow=size,ncol=size)
  print(Vecstates)
  print(states)
  Vecstates[1,] = states$StateSpace_Mat[1,]
  
  # 
  while (i<ncol(Q)){
    lambda = -Q[i,i]
    t[count] = log(1-runif(1))/(-lambda) 

    VecNext = Q[i,]/(sum(Q[i,])-Q[i,i])
    VecNext[i] = 0
    
    inext = sample(x = 1:size, size = 1, prob = VecNext)
    Vecstates[count+1,] = states$StateSpace_Mat[inext,]
    i = inext
    count = count + 1
  }
  return(list(t,Vecstates))
}

# exemple for n=4 of 'simul'

Q = matrix(c(-6,6,0,0,0,
             0,-3,2,1,0,
             0,0,-1,0,1,
             0,0,0,-1,1), nrow = 4, byrow = T)
a = c(1,0,0,0)
simul(Q,a)



# MakeReward -----------------------------------------
# Aim:
#   Creates the reward vector for the total branch length
#
# Input:
#   A matrix 'block' of all the possible states
#   A vector 'vecwhich', only if we want to consider specific
#     | branch states 
#
# Ouput:
#   A vector 'reward' with the reward for each state

MakeReward = function(block,vecwhich = 1:ncol(block)){
  reward = block[,vecwhich]
  reward = apply(block,1,sum)
  return(reward)
}

# Exemple of 'MakeReward'

block = BlockCountProcess(nrow(Q))
block = block$StateSpace_Mat
vecwhich = c(1,2)

MakeReward(block)
MakeReward(block, vecwhich)  # We do not consider branch with 3 lineages



# theta_est ___ (to modify) -----------------------------------------
# Aim:
#   Estimates Theta (mutation rate) with the data
#
# Input:
#   A vector 'xton', with the value of site frequency we want to use
#   A vector 'vecxton' with the corresponding site frequency spectrum
#
# Ouput:
#   The estimates of Theta

theta_est = function(vecxton, xton = 1:length(vecxton)){
  theta = sum(xton)/(sum(1/vecxton[which(vecxton!=0)]))
  return(theta)
}

# exemple of theta_est:

vecxton=c(12,8,2,3,1,0,0,1,0)

theta_est(vecxton)
theta_est(vecxton[3:8], 3:8)


# AddMutation -------------------------------

mutfromtree = function(t, vecStates){
  
}