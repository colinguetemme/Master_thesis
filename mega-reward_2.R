require(tidyr)

subint_mat = matrix(c(-3,0,0,3,-2,0,0,1,-1), ncol=3)
init_probs =c(0.6,0.2,0.2,0)
R = c(2,0,5,1)
n = length(init_probs)

Reward_zerone <- function(DPH = NULL, init_probs = NULL, subint_mat = NULL,
                          reward_vec = NULL){
  if (is.null(DPH)){
    if (is.vector(init_probs) && is.matrix(subint_mat)){
      if (length(init_probs) == ncol(subint_mat) && 
          ncol(subint_mat) == nrow(subint_mat)){
        
        setofmat = rep(list(as.list(1:n)),n)

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

          }
        }
        
        p = which(R>0)
        z = which(R==0)
        
        vp = c(0,sumvec(R))
        vz = c(0,sumvec(size * (as.numeric(R==0))))
        WESH = NULL
        
        T_tilde = list(list(p,p), list(p,z), list(z,p), list(z,z))
        count = 1
        
        for (i in T_tilde){
          combn = as.matrix(expand.grid(i[[1]],i[[2]]))
          subT = matrix(0, ncol = sum(size[i[[1]]]),
                              nrow = sum(size[i[[2]]]))

          ifelse(i[[1]]==p, abs_pos1 <- vp, abs_pos1 <- vz)
          ifelse(i[[2]]==p, abs_pos2 <- vp, abs_pos2 <- vz)
          
          for (j in 1:nrow(combn)){
            
            selec_combn = as.vector(combn[j,])
            numrow = NULL; numcol = NULL
            numcol = (abs_pos1[selec_combn[1]]+1):
              (abs_pos1[selec_combn[1]+1])
            numrow = (abs_pos2[selec_combn[2]]+1):
              (abs_pos2[selec_combn[2]+1])


            subT[numrow, numcol] = setofmat[[selec_combn[2]]][[selec_combn[1]]]
            }
          WESH[[count]] = subT
          count = count+1
        }
        new_subint_mat = WESH[[1]] + (WESH[[3]] %*% solve(diag(1,ncol(WESH[[4]]))-WESH[[4]]) %*% WESH[[2]])
      }else{print('Wrong dimensions')}
    }else{print()}
  }
}

sumvec = function(vec){
  for (i in 2:length(vec)){
    vec[i] = vec[i-1]+vec[i]
  }
  return(vec)
}

