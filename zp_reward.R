require(tidyr)

Reward_zerone <- function(DPH = NULL, init_probs = NULL, subint_mat = NULL,
                         reward_vec = NULL){
  if (is.null(DPH)){
    if (is.vector(init_probs) && is.matrix(subint_mat)){
      if (length(init_probs) == ncol(subint_mat) && 
          ncol(subint_mat) == nrow(subint_mat)){
        p = which(R==1)
        z = which(R==0)
        
        ptop <- as.matrix(subint_mat[p,p])
        ztoz <- as.matrix(subint_mat[z,z])
        ztop <- as.matrix(subint_mat[z,p])
        ptoz <- as.matrix(subint_mat[p,z])
        if(length(p)<length(z)){ptoz = t(ptoz)}
        if(length(z)<length(p)){ztop = t(ztop)}
        
        new_subint_mat = ptop + (ptoz %*% solve(diag(1,ncol(ztoz))-ztoz) %*% ztop)
        
        init_probs_p = init_probs[p]
        init_probs_z = init_probs[z]
        
        new_init_probs = init_probs_p + (init_probs_z %*% solve(diag(1,ncol(ztoz))-ztoz)%*% ztop)
        
        
      }else{print('Wrong dimensions')}
    }else{print()}
  }
}


