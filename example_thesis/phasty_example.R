
## -- R Code example of the report -- ##

# Make sure to have phasty in your package list
library(phasty)

# -- Section 2.5 --

## For discrete phase-type
subint_mat <- matrix(c(0, 0.2, 0.8,
                       0.5, 0.5, 0,
                       0, 0, 0.4), ncol = 3,
                     byrow = T)
init_probs <- c(0.7, 0.3, 0)
dph <- phase_type(subint_mat, init_probs)

print(dph)

## For continuous phase-type
subint_mat <- matrix(c(-3, 1.5, 1.5,
                       0, -2, 1,
                       1, 0, -2), ncol = 3,
                     byrow = T)
init_probs <- c(0.9, 0.1, 0)
ph <- phase_type(subint_mat, init_probs)
class(ph)

## First moments

mean(dph)
mean(ph)
var(dph)
var(ph)

## Generic functions

round(dphtype(c(0, 1, 2, 3), dph), 4)
round(pphtype(c(0, 1, 2, 3), dph), 4)
qphtype(c(0.05, 0.50, 0.95), dph)

set.seed(0)
rphtype(10, dph)
set.seed(0)
round(rphtype(10, ph), 4)

### One plot example

x <- seq(0, qphtype(0.95, ph), 0.1)
pdf <- dphtype(x, ph)
cdf <- pphtype(x, ph)

plot(x, cdf, col = 'blue')
points(x, pdf, col = 'orange')

# -- Section 2.8 --

rew_vec <- c(0, 2, 1)
rph <- reward_phase_type(ph, rew_vec)
print(rph)

rew_mat <- matrix(c(0.25, 0.75, 0, 0,
                    0, 0.5, 0.5, 0,
                    0, 0.25, 0, 0.75), ncol = 4,
                  byrow = TRUE)
rdph <- reward_phase_type(dph, rew_mat)
print(rdph)

# -- Section 2.9.1 --

reward_mat <- matrix(c(1, 0, 1,
                       2, 4, 1,
                       3, 2, 1), ncol = 3,
                     byrow = T)

mdph <- phase_type(dph$subint_mat, dph$init_probs,
                   reward_mat)

mean(mdph)
var(mdph)

# -- Section 2.9.2 --

reward_mat <- matrix(c(0.2, 0, 5.2,
                       0.7, 4, 1,
                       0.9, 2, 1.3), ncol = 3,
                     byrow = T)

mph <- phase_type(ph$subint_mat, ph$init_probs,
                   reward_mat)

mean(mph)
var(mph)

# -- Section 4.1 --

# for any integer n > 2
T_mrca_n <- function(n) {
    T_mat <- diag(-choose(n:2, 2))
    T_mat[-(n-1), -1] <- T_mat[-(n-1), -1] -
        T_mat[-(n-1), -(n-1)]
    return(phase_type(T_mat))
}

T_mrca_4 <- T_mrca_n(4)
mean(T_mrca_4)
var(T_mrca_4)

# for any integer n > 2
T_total_n <- function(n) {
    return(reward_phase_type(T_mrca_n(n), 2:n))
}

T_total_4 <- T_total_n(4)
mean(T_total_4)
var(T_total_4)

# -- Section 4.2 --

T_mrca_b <- function(n, a){
    P_mat = matrix(0, n-1, n)
    for (i in 1:(n-1)){
        for (j in (i+1):(n)){
            P_mat[i, j] <- (beta(j - i + 1 - a, n - j + a) /
                                beta(a, 2 - a)) *
                choose(n + 1 - i, 2)
        }
    }
    T_mat <- P_mat[,-n]
    diag(T_mat) <- -rowSums(P_mat)
    return(phase_type(T_mat))
}

T_mrca_b_4 <- T_mrca_b(4, 1.5)
mean(T_mrca_b_4)
var(T_mrca_b_4)

# -- Section 4.4 --

## BLOCK COUNTING PROCESS FUNCTION, NOT INCLUDE IN THE PAPER ##

block_counting_process <- function(n){
    ##----------------------------------------------------
    ## Possible states
    ##----------------------------------------------------
    ## Size of the state space (number of states)
    nSt <- partitions::P(n)
    ## Definition of the state space
    StSpM <- matrix(ncol = n, nrow = nSt)
    ## Set of partitions of [n]
    x <- partitions::parts(n)
    ## Rewriting the partitions as (a1,...,an)
    for (i in 1:nSt) {
        st <- x[, i]
        StSpM[i,] <- tabulate(x[, i], nbins = n)
    }
    ## Reordering
    StSpM <- StSpM[order(rowSums(StSpM), decreasing = TRUE),]
    ## Because of this ordering we can't 'go back', i.e.
    ## below the diagonal the entries are always zero
    ##----------------------------------------------------
    ## Intensity matrix
    ##----------------------------------------------------
    RateM <- matrix(0, ncol = nSt, nrow = nSt)
    ## Algorithm for finding rates between states
    for (i in 1:(nSt - 1)){
        for (j in (i + 1):nSt){
            # cat(i," state i",StSpM[i,])
            # cat(" ",j," state j",StSpM[j,])
            cvec <- StSpM[i, ] - StSpM[j,]
            # cat(" cvec",cvec)
            ## Two branches are merged, i.e. removed from state i
            check1 <- sum(cvec[cvec > 0]) == 2
            # cat(" check1",check1)
            ## One new branch is created, i.e. added in state from j
            check2 <- sum(cvec[cvec < 0]) == -1
            # cat(" check2",check2)
            if (check1 & check2){
                ## Size(s) of the block(s) and the corresponding rates
                tmp <- StSpM[i, which(cvec > 0)]
                RateM[i, j] <- ifelse(length(tmp) == 1, tmp * (tmp - 1) / 2, prod(tmp))
            }
        }
    }
    ## Diagonal part of the rate matrix
    for (i in 1:nSt){
        RateM[i, i] <- -sum(RateM[i,])
    }
    return(list(rate_mat = RateM[-nrow(RateM), -ncol(RateM)],
                state_mat = StSpM[-nrow(StSpM),-ncol(StSpM)]))
}

###########################################

block_counting_process(4)

# -- Section 4.4.1 --

discretize <- function(subint_mat, theta){
    S_mat <-solve(diag(ncol(subint_mat)) - 2/theta * subint_mat)
    return(S_mat)
}

S_sites <- function(n, theta){
    T_total <- T_total_n(n)
    S_mat <- discretize(T_total$subint_mat, theta)
    return(phase_type(S_mat))
}


S_4 <- S_sites(4, 2)
print(S_4)
mean(S_4)
var(S_4)

# -- Section 4.4.2 --

i_ton <- function(n, i_ton, theta = 10){
    bcp <- block_counting_process(n)
    Ti <- phase_type(discretize(bcp$rate_mat, theta))
    PHi_ton <- reward_phase_type(Ti, bcp$state_mat[, i_ton])
    return(PHi_ton)
}

SFS <- function(n, theta = 10){
    bcp <- block_counting_process(n)
    Ti <- discretize(bcp$rate_mat, theta)
    return(phase_type(Ti, reward_mat = bcp$state_mat))
}

SD <- function(n, theta = 10){
    bcp <- block_counting_process(n)
    Ti <- discretize(bcp$rate_mat, theta)
    reward_mat <- 0 * bcp$rate_mat
    for (i in 1:nrow(bcp$state_mat)){
        reward_mat[i, 1:bcp$state_mat[i, 1]] <- 1 
    }
    return(phase_type(Ti, reward_mat = reward_mat))
}

bcp_4 <- block_counting_process(4)
S_sites(4, 10)

SFS_4 <- SFS(4)

mean(SFS_4)
var(SFS_4)

bcp_10 <- block_counting_process(10)
S_sites(n = 10, theta = 2)
SFS_10 <- phase_type(bcp_10$rate_mat,
                    reward_mat = bcp_10$space_mat)
plot(mean(SFS_10))
