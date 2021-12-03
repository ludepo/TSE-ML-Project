# Objective : Hand-code estimators needed
# Created by: Nikita Marini & Luca Poll
# Created on: 4/10/2021


library(pracma)

# get the alpha range for the respective data generating process and estimator
get_alpha <- function(data, estimator){

    DGP <- data$DGP

    # get alpha range for ridge estimator
    if (estimator == "ridge"){
      if (DGP %in% c(1,6)){
        alpha <- seq(0.001, 0.1, 0.001)
      } else if (DGP == 2){
        alpha <- seq(0.0001, 0.01, 0.0001)
      } else if (DGP == 3){
        alpha <- seq(0.0005, 0.1, 0.0005)
      } else if (DGP == 4){
		alpha <- seq(0.01, 1, 0.01)
	  } else if (DGP == 5){
        alpha <- seq(0.001, 0.15, 0.001)
      }
	# account for Î±/N
    alpha <- alpha * ncol(data$X)
    }

    # get alpha range for Landweber Fridman estimator
    if (estimator == "LF"){
      if (DGP == 1){
        alpha <- seq(0.000001, 0.0003, 0.00002)
      } else if (DGP == 2){
        alpha <- seq(0.00001, 0.0002, 0.00001)
      } else if (DGP == 3){
        alpha <- seq(0.000001, 0.00005, 0.0000025)
      } else if (DGP == 4){
        alpha <- seq(0.000001, 0.0004, 0.00001)
      } else if (DGP == 5){
        alpha <- seq(0.000001, 0.0004, 0.00002)
      } else if (DGP == 6){
        alpha <- seq(0.000001, 0.016, 0.001)
      }
	# account for Î±/N
    alpha <- alpha * ncol(data$X)
    }

	if (estimator == "PC"){
		alpha <- 1:data$rmax
	}

    if (estimator == "PLS"){
		if (data$T == 50){
			alpha <- seq(1, 25)
		} else if (data$T == 500){
			alpha <- seq(1,15)
		}
    }

    return(alpha)

}


# function including all estimators needed
estimate <- function(data, estimator, alpha){

	# Simplify notation
	T <- nrow(data$X)
	N <- ncol(data$X)
	X <- data$X
	Y <- data$Y
	k <- alpha

	# Compute all the necessary elements for the estimation of delta
	Sxx <- (t(X) %*% X)/T
	Sxy <- (t(X) %*% Y)/T
	temp.eigen <- eigen(Sxx)
	lambda <- temp.eigen$values
	phi <- temp.eigen$vectors
	d <- 0.018/max(lambda)


	XXt <- (X %*% t(X))/T
	temp.eigen <- eigen(XXt)
	psi <- temp.eigen$vectors*sqrt(T)

	# Compute PLS hat matrix according to equation (11)
	if (estimator == "PLS"){
		Vk <- matrix(nrow = N, ncol = k)
		Vk[,1] <- t(X) %*% Y
		if (k > 1){
			Vk[,2] <- t(X) %*% X %*% t(X) %*% Y
		}
		if (k > 2) {
			for (l in 3:k){
				Vk[,l] <- (t(X) %*% X) %*% Vk[,l-1]
			}
		}

		M <- X %*% Vk %*% pinv((t(Vk) %*% t(X) %*% X %*% Vk)) %*% t(Vk) %*% t(X)
		My <- M %*% Y
		tuning <- "k"
		tp <- k
		DoF <- NA
		delta <- NA
	}
	else {
		# Determine shrinkage parameter q_j based on estimator (excluding PLS)
		q <- c()

		if (estimator == "ridge"){
			for (j in 1:length(lambda)){
				q[j] <- lambda[j]/(lambda[j] + alpha)
			}
			tuning <- "alpha"
		}

		if (estimator == "LF"){
			for (j in 1:length(lambda)){
				q[j] <- 1-(1-d*lambda[j])^(1/alpha)
			}
			tuning <- "alpha"
		}

		if (estimator == "PC"){
			for (j in 1:length(lambda)){
				ifelse(j <= alpha, q[j] <- 1, q[j] <- 0)
			}
			tuning <- "k"
		}

		# print(q)
		# Compute hat matrix
		M <- list()

		for (j in 1:min(N, T)){
		 M[[j]] <- q[j]*matrix(psi[,j])%*%t(matrix(psi[,j]))
		}

		delta <- list()
		for (j in 1:min(N, T)){
		 delta[[j]] <- q[j]/lambda[j]*(t(Y)%*%psi[,j]/T)%*%phi[,j]
		}

		M <- Reduce('+', M)/T
		My <- M %*% Y

		delta <- matrix(Reduce('+', delta))

		DoF <- sum(diag(M))
		if (tuning == "k"){
			tp <- k
		} else if (tuning == "alpha"){
			tp <- alpha
		}
	}

	# Output results
	ests <- list(M = M, My = My, alpha = alpha, tuning = tuning, tp = tp,
	             DoF = DoF, delta = delta, eigen = temp.eigen)
	return(ests)

}

