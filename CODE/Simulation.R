# Objective : Replicate simulations presented in Carrasco & Rossi (2016)
# Created by: Luca Poll
# Created on: 4/4/2021



## =============================================================================
## ==  Set environment      ====================================================
## =============================================================================

# setwd("..")
datapath <- paste0(getwd(), "/DATA/")
simoutpath <- paste0(getwd(), "OUTPUT/SIMULATION/")



## =============================================================================
## ==  Simulation setup     ====================================================
## =============================================================================

# F_t is (r x 1) vector of iid N(0,I)
# e_t is (N x 1) vector of iid N(0,I) random variables
# A is (N x r) matrix of iid N(0,I)
# x_t = A*F_t + e_t

# in matrix notation: X = F * A' + e
# where X = (T x N), F = (T x r), A' = (r x N), e = (T x N)

# theta is (r x 1) vector
# v_t is iid N(0, sig_v^2) and sig_v^2 = 1
# y_t = theta'F_t + v_t

# in matrix notation: y = F*theta + v
# where y = (T x 1), F = (T x r), theta = (r x 1), v = (T x 1)

# - - - - - - - - - - - - - - - - - - - -
# consider two samples:
# - large sample case: N=200, T=500
# - small sample case: N=100, T=50

# - - - - - - - - - - - - - - - - - - - -
# consider different DGPs:
# - DGP 1 (few factors structure): r = 4 in theta; r_{max} = r + 10
# - DGP 2 (many factors structure): r = 50 in theta; r_{max} = min(N, T/2)
# - DGP 3 (five factors, one relevant): r = 5; r_{max} = min(r + 10, min(N, T/2));
#          theta = (1, 0_{1x4})' (y depends on only one factor, x on five);
#          F = [F1', F2']'; F1 relevant, F2 (4 factors) not, factors are uncorrelated
#          F1 has unit variance, F2 diagonal cov matrix with coefs (2,3,3,4)
# - DGP 4 (x_t factor structure but unrelated to y_t) theta vector of zeros; r=5,
#          r_{max} = r+10; factors uncorrelated diagonal cov matrix: (1,2,3,3,4)
# - DGP 5 (Eigenvalues declining slowly): r = N; r_{max}=min(N, T/2), theta (N x 1)
#          vector of ones; A = M o e (see paper)
# - DGP 6 (Near factor model): theta = 1, r = 1, r_{max} = r+10, A' = n^{-1/2} 1_{rxN}


## =============================================================================
## ==  Generate data        ====================================================
## =============================================================================

# create sampler function
sampler <- function (small_sample = F){

  # define parameters
  if (small_sample){
    N <- 100; T <- 50
  } else {
    N <- 200; T <- 500
  }
  rs <- c(4,50,5,5,N,1)
  rmaxs <- c(rs[1]+10, min(N,T/2), min(rs[3], min(N, T/2)), rs[4]+10, min(N, T/2), rs[6]+10)
  names <- c("DGP1", "DGP2", "DGP3", "DGP4", "DGP5", "DGP6")
  samp <- list()

  # get variables
  for (i in 1:6){
    r <- rs[i]
    rmax <- rmaxs[i]

    # simulate F
    F <- matrix(rnorm(r*T), nrow = T, ncol = r)
    if (i %in% c(3, 4)){
      F <- matrix(c(rnorm(T), rnorm(T, sd=sqrt(2)),
                    rnorm(T, sd=sqrt(3)), rnorm(T, sd=sqrt(3)),
                    rnorm(T, sd=sqrt(4))),
                  nrow = T, ncol = 5)
    }

    # Simulate A'
    Aprime <- matrix(rnorm(r*N), nrow = r, ncol = N)
    if (i == 5){
      M <- matrix(NA, nrow = r, ncol = N) # in DGP5 r=N
      for (n in 1:N){
        M[n,] <- 1/n}
      Xi <- matrix(rnorm(r*N), nrow = r, ncol = N)
      A <- M*Xi
      Aprime <- t(A)
    }
    if (i == 6){
      Aprime <- matrix(rep(1/sqrt(N), r*N), nrow = r, ncol = N)
    }

    # Simulate X
    e <- matrix(rnorm(T*N), nrow = T, ncol = N)
    X <- F %*% Aprime + e

    # define theta vector
    if (i < 3 | i > 4){
      theta <- matrix(rep(1, r), nrow = r, ncol = 1)
    } else if (i == 3){
      theta <- matrix(c(1, rep(0, r-1)), nrow = r, ncol = 1)
    } else if (i == 4){
      theta <- matrix(rep(0, r), nrow = r, ncol = 1)
    }

    # simulate Y
    v <- rnorm(T)
    if (i == 3){
      v <- rnorm(T, 0, 0.1)
    }
    Y <- F %*% theta + v

    # put in list
    dgp <- list(dgp=i, r, rmax, N, T, Y, F, theta, v, X, Aprime, e)
    names(dgp) <- c("DGP", "r", "rmax", "N", "T", "Y", "F", "theta", "v", "X", "Aprime", "e")
    samp[i] <- list(dgp)

    names(samp)[i] <- names[i]
  }
  return(samp)
}





