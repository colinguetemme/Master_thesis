require(tidyr)

subint_mat = matrix(c(1,2,3,4,5,6,7,8,9), ncol=3)
init_probs =c(0.6,0.2,0.2)
R = c(2,0,5)
n = length(init_probs)

Reward_zerone <- function(DPH = NULL, init_probs = NULL, subint_mat = NULL,
                          reward_vec = NULL){
  if (is.null(DPH)){
    if (is.vector(init_probs) && is.matrix(subint_mat)){
      if (length(init_probs) == ncol(subint_mat) && 
          ncol(subint_mat) == nrow(subint_mat)){
        
        setofmat = rep(list(as.list(1:n)),n)
        matsize = rep(list(as.list(1:n)),n)
        size = R
        size[which(size==0)] = 1
        for (i in 1:n){
          for (j in 1:n){
            matij = matrix(0, nrow = size[i], ncol = size[j])
            matij[size[i], 1] = subint_mat[i, j]
            
            if(i == j){matij[-size[i],-1] <- diag(1,size[i]-1)}

            print(c(i,j))
            print(matij)
            setofmat[[i]][[j]] = matij
            matsize[[i]][[j]] = c(R[i], R[j])
          }
        }
        matsize = t(as.data.frame(matsize))
        T_tilde = matrix(0, ncol=length(size), nrow = length(size))
        # for T++
        matsize[which(matsize[,1]>0 & matsize[,2]>0),]
        T_tilde_pp = as.matrix(expand.grid(which(R>0),which(R>0)))
        T_tilde_pz = as.matrix(expand.grid(which(R>0),which(R==0)))
        T_tilde_zp = as.matrix(expand.grid(which(R==0),which(R>0)))
        T_tilde_zz = as.matrix(expand.grid(which(R==0),which(R==0)))
        positionmat = as.matrix(expand.grid(1:2,1:2))
        
        position_pp = matrix(size[T_tilde_pp]+1-min(size[T_tilde_pp]),ncol=2)
        
        
        newmat = setofmat[[rewardvec[1,1]]][[rewardvec[1,2]]]
        for (i in 2:nrow(reward_vec)){
          
        }
        
        
      }else{print('Wrong dimensions')}
    }else{print()}
  }
}


