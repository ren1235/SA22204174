#' @title A mixture model approach for cluster analysis using R
#' @description A mixture model approach to estimate the expectations and covariance matrices of each cluster using R
#' @param X the random samples generated from a Gaussian mixture model, data structure: (n,p)-dimensional matrix.
#' @param K the number of mixture components in the population, data structure: int.
#' @param times the upper bound on the number of iterations, data structure: int.
#' @param initial.tau initial value of tau, data structure: K-dimensional vector.
#' @param initial.mu initial value of mu, data structure: (p,K)-dimensional matrix, or K-dimensional vector when p=1.
#' @param initial.sigma initial value of sigma, data structure: (p,p,K)-dimensional array, or K-dimensional vector when p=1.
#' @return a list including
#'  iterations: the number of iterations, data structure: int.
#'  tau: the estimate of the proportion of each cluster, data structure: K-dimensional vector.
#'  mu: the estimate of expectation in each cluster, data structure: (p,K)-dimensional matrix.
#'  sigma: the estimate of covariance matrix in each cluster, data structure: (p,p,K)-dimensional array.
#' @examples
#' \dontrun{
#'     #parameters of the samples
#'     n <- 500; p <- 2; K <- 3
#'     mu <- matrix(c(-2,-1,0,4,3,1), nrow = p, ncol = K)
#'     sigma <- matrix(0.3, nrow = p, ncol = p) + diag(0.7, p)
#'     #generate the random samples
#'     X <- matrix(0, nrow = n, ncol = p)
#'     set.seed(1)
#'     select <- sample(1:3, size = n, replace = TRUE, prob = c(1/3,1/3,1/3))
#'     X[select==1,] <- mvtnorm::rmvnorm(sum(select==1), mean = mu[,1], sigma = sigma)
#'     X[select==2,] <- mvtnorm::rmvnorm(sum(select==2), mean = mu[,2], sigma = sigma)
#'     X[select==3,] <- mvtnorm::rmvnorm(sum(select==3), mean = mu[,3], sigma = sigma)
#'     #initialization
#'     initial.tau <- rep(1/K, K)
#'     initial.mu <- matrix(c(-2,-2,0,0,2,2), nrow = p, ncol = K)
#'     initial.sigma <- array(rep(cov(X), K), dim = c(p,p,K))
#'     #cluster analysis
#'     mixmodel_R(X, K, times = 500, initial.tau, initial.mu, initial.sigma)
#' }
#' @export
mixmodel_R <- function(X, K, times = 1e4, initial.tau, initial.mu, initial.sigma){
  #X: the random samples generated from a Gaussian mixture model, data structure: n*p-dimensional matrix
  #K: the number of mixture components in the population, data structure: int
  #times: the upper bound on the number of iterations, data structure: int
  #initial.tau: initial value of tau, data structure: K-dimensional vector
  #initial.mu: initial value of mu, data structure: p*K-dimensional matrix, K-dimensional vector(p=1)
  #initial.sigma: initial value of sigma, data structure: p*p*K-dimensional array, K-dimensional vector(p=1)
  n <- nrow(X); p <- ncol(X)#n: sample size, p: the dimension of random vectors
  tau <- matrix(0, nrow = K, ncol = times+1)#the chains from iteration
  mu <- array(0, dim = c(p, K, times+1))
  sigma <- array(0, dim = c(p, p, K, times+1))
  z <- matrix(0, nrow = n, ncol = K)
  rate <- 1e-8 #Convergence is checked when the change of elements in the sequence is small enough
  
  #initialization
  tau[,1] <- initial.tau
  mu[,,1] <- initial.mu
  sigma[,,,1] <- initial.sigma
  
  #EM iteration
  for(i in 1:times){
    for(k in 1:K){
      #      zk <- apply(X, 1, function(x) exp(-(x-mu[,k,i])%*%solve(sigma[,,k,i])%*%(x-mu[,k,i])/2)/sqrt(2*pi*det(sigma[,,k,i])) )
      if(p > 1){
        zk <- apply(X, 1, function(x) mvtnorm::dmvnorm(x, mean = mu[,k,i], sigma = sigma[,,k,i]))
        }else{
          zk <- dnorm(X, mean = mu[,k,i], sd = sqrt(sigma[,,k,i]))
        }#n-dimension
      z[,k] <- tau[k,i]*zk#n-dimension
    }
    z_column <- apply(z, 1, sum)#n-dimension
    z <- z/z_column
    nhat <- apply(z, 2, sum)#K-dimension
    tau[,i+1] <- nhat/n
    for(k in 1:K){
      mu[,k,i+1] <- apply(z[,k]*X, 2, sum)/nhat[k]#p-dimension
      sigma[,,k,i+1] <- t(X-matrix(1,n,1)%*%mu[,k,i+1])%*%diag(z[,k])%*%(X-matrix(1,n,1)%*%mu[,k,i+1])/nhat[k]#(p*n)*(n*n)*(n*p)=p*p-dimension
    }
    #cat("step=", i, "tau=", tau[,i+1], "mu=", mu[,,i+1], "sigma=", sigma[,,,i+1], '\n')
    if (max( sum(abs(tau[,i+1] - tau[,i]))/sum(abs(tau[,i])), sum(abs(mu[,,i+1]-mu[,,i]))/sum(abs(mu[,,i])), sum(abs(sigma[,,,i+1]-sigma[,,,i]))/sum(abs(sigma[,,,i])) ) < rate){
      #cat("The sequence converges after", i, "iterations.", '\n')
      break
    }
  }
  return(list(iterations = i, tau = tau[,i+1], mu = mu[,,i+1], sigma = sigma[,,,i+1]))
}

#' @title A concave pairwise fusion approach for subgroup analysis using R
#' @description A concave pairwise fusion approach to estimate the parameters of a linear statistical model in order to do subgroup analysis using R
#' @param X independent variables of linear regression model, data structure: (n,p)-dimensional matrix.
#' @param y dependent variables of linear regression model, data structure: (n,1)-dimensional matrix, or n-dimensiobal vector.
#' @param times the upper bound on the number of iterations, data structure: int.
#' @param vartheta parameter of the penalty function, data structure: numeric.
#' @param gamma parameter of the penalty function, data structure: numeric.
#' @param lambda parameter of the penalty function, data structure: numeric.
#' @return a list including
#'  iterations: the number of iterations, data structure: int.
#'  mu: the estimate of mu in linear statistical model, data structure: n-dimensional vector.
#'  sigma: the estimate of beta in linear statistical model, data structure: (p,1)-dimensional matrix.
#' @examples
#' \dontrun{
#'     #parameters of the samples
#'     n <- 50; p <- 2; K <- 3
#'     sigma <- matrix(0.3, nrow = p, ncol = p) + diag(0.7, p)
#'     mu <- c(-2,0,2)
#'     beta <- matrix(c(1,3), nrow = p, ncol = 1)
#'     #generate the samples
#'     set.seed(1)
#'     mu.sample <- sample(mu, size = n, replace = TRUE, prob = c(1/3,1/3,1/3))
#'     X <- mvtnorm::rmvnorm(n, mean = c(0,0), sigma = sigma)
#'     y <- mu.sample + X%*%beta + rnorm(n, mean = 0, sd =0.5)
#'     #parameters to be input in the function
#'     vartheta <- 1; gamma <- 3; lambda <- 0.05; times <- 1e4
#'     #using the function to implement subgroup analysis
#'     result <- concavefusion(X, y, times, vartheta, gamma, lambda)
#'     #check whether the results accord with the theoretical value
#'     result$iterations; 
#'     result$mu[mu.sample==mu[1]]; result$mu[mu.sample==mu[2]]; result$mu[mu.sample==mu[3]]; 
#'     result$beta
#' }
#' @export
concavefusion <- function(X, y, times = 1e6, vartheta, gamma, lambda){
  #X: independent variables of linear regression model, data structure: n*p-dimensional matrix
  #y: dependent variables of linear regression model, data structure: n*1-dimensional matrix, n-dimensional vector
  #times: the upper bound on the number of iterations, data structure: int
  #vartheta, gamma, lambda: parameters of the penalty function, data structure: numeric
  n <- nrow(X) #sample size
  rate <- 1e-8 #Convergence is checked when the change of elements in the sequence is small enough
  
  #the deterministic matrices that will be used in the iteration
  deltaT <- matrix(0, nrow = n, ncol = n*(n-1)/2)
  #deltaT[1,1:(n-1)] <- 1
  #for(j in 2:n){deltaT[j,j-1] <- -1}
  #for(i in 2:(n-1)){para_delta[i,(sum((n-i+1):(n-1))+1):(sum((n-i):(n-1)))] <- 1}
  #for(j in 2:n){para_delta[j,j:n] <- -1}
  for(i in 1:(n-1)){
    for(j in 1:(n-i)){
      deltaT[i,(2*n-i)*(i-1)/2+j] <- 1
      deltaT[i+j,(2*n-i)*(i-1)/2+j] <- -1
    }
  }
  deltaTdelta <- tcrossprod(deltaT)
  Qx <- X %*% solve(crossprod(X)) %*%t(X)
  para.mu <- solve(vartheta*deltaTdelta + diag(1, n) - Qx)
  para.beta <- solve(crossprod(X))
  
  #initialization
  X1 <- cbind(1, X)
  initial <- solve(crossprod(X1)) %*% (crossprod(X1, y))
  mupast <- rep(initial[1,], n)
  betapast <- initial[-1,]
  deltanew <- etanew <- etapast <- upsilonnew <- upsilonpast <- matrix(0, nrow = n, ncol = n)
  
  #iteration
  for(t in 1:times){
    #update mu
    etapast.vector <- matrix(0, 1, 1)
    upsilonpast.vector <- matrix(0, 1, 1)
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        etapast.vector <- rbind(etapast.vector, etapast[i,j])
        upsilonpast.vector <- rbind(upsilonpast.vector, upsilonpast[i,j])
      }
    }
    etapast.vector <- etapast.vector[-1,]
    upsilonpast.vector <- upsilonpast.vector[-1,]
    
    munew <- para.mu %*% ((diag(1, n)-Qx) %*% y + vartheta*deltaT %*% (etapast.vector-upsilonpast.vector/vartheta))
    
    #update beta
    betanew <- para.beta %*% crossprod(X, (y-munew))
    
    #update eta
    for(i in 1:n){
      for(j in 1:n){
        deltanew[i,j] <- munew[i] - munew[j] + upsilonpast[i,j]/vartheta
        if(abs(deltanew[i,j]) <= gamma*lambda){
          etanew[i,j] <- sign(deltanew[i,j])*max(abs(deltanew[i,j])-lambda/vartheta,0) / (1-1/gamma/vartheta)
        }else{
          etanew[i,j] <- deltanew[i,j]
        }
      }
    }
    
    #update upsilon
    for(i in 1:n){
      for(j in 1:n){
        upsilonnew[i,j] <- upsilonpast[i,j] + vartheta*(munew[i]-munew[j]-etanew[i,j])
      }
    }
    
    #Check for convergence
    if (max(sum(abs(munew - mupast))/(sum(abs(mupast))+1e-4), sum(abs(betanew-betapast))/(sum(abs(betapast))+1e-4), sum(abs(etanew-etapast))/(sum(abs(etapast))+1e-4), sum(abs(upsilonnew-upsilonpast))/(sum(abs(upsilonpast))+1e-4)) < rate){
      #cat("The sequence converges after", t, "iterations.", '\n')
      break
    }else{
      mupast <- munew
      betapast <- betanew
      etapast <- etanew
      upsilonpast <- upsilonnew
    }
  }
  return(list(iterations = t, mu = munew, beta = betanew))
}

#' @title Benchmark two approaches.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of mixture model approach (R function \code{mixmodel_R} and Rcpp function \code{mixmodel_Rcpp}) and concave pairwise fusion approach (R function \code{concavefusion}).
#' @examples
#' \dontrun{
#'     #parameters of the samples
#'     n <- 30; p <- 1; K <- 3
#'     sigma <- 1
#'     mu <- c(-2,0,3)
#'     beta <- 2
#'     #generate the samples
#'     set.seed(1)
#'     mu.sample <- sample(mu, size = n, replace = TRUE, prob = c(0.2,0.3,0.5))
#'     X <- rnorm(n, mean = 0, sd = sqrt(sigma))
#'     X <- matrix(X, nrow = n, ncol = 1)
#'     y <- mu.sample + X%*%beta + rnorm(n, mean = 0, sd =0.5)
#'     X1 <- cbind(1, X)
#'     betahat <- (solve(crossprod(X1))%*%crossprod(X1,y))[2,1]
#'     Xhat <- y-X*betahat
#'     #initialization or parameters
#'     initial_tau <- rep(1/K, K)
#'     initial_mu <- c(-1, 0, 1)
#'     initial_sigma <- rep(var(Xhat),K)
#'     vartheta <- 1; gamma <- 3; lambda <- 0.05; times <- 1e4
#'     #using the function to implement subgroup analysis
#'     mixmodel_R(Xhat, K, 1e5, initial_tau, initial_mu, initial_sigma)
#'     mixmodel_Rcpp(c(Xhat), K, 1e5, initial_tau, initial_mu, initial_sigma)
#'     concavefusion(X, y, times, vartheta, gamma, lambda)
#'     #compare the speed of these three functions
#'     compare.time <- microbenchmark::microbenchmark(
#'     mixmodel_R = mixmodel_R(Xhat, K, 1e5, initial_tau, initial_mu, initial_sigma), 
#'     mixmodel_Rcpp = mixmodel_Rcpp(c(Xhat), K, 1e5, initial_tau, initial_mu, initial_sigma), 
#'     concavefusion = concavefusion(X, y, times, vartheta, gamma, lambda))
#'     summary(compare.time)[,c(1,3,5,6)]
#' }
NULL

#' @import microbenchmark
#' @import mvtnorm
#' @import knitr
#' @import pander
#' @import bootstrap
#' @import boot
#' @import DAAG
#' @import coda
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnorm
#' @useDynLib SA22204174
NULL
