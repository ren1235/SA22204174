#include <Rcpp.h>
using namespace Rcpp;
//' @title A mixture model approach for cluster analysis using Rcpp when p=1
//' @description A mixture model approach to estimate the expectations and covariance matrices of each cluster using Rcpp when p=1
//' @param X the random samples generated from a Gaussian mixture model, data structure: n-dimensional vector.
//' @param K the number of mixture components in the population, data structure: int.
//' @param times the upper bound on the number of iterations, data structure: int.
//' @param initial_tau initial value of tau, data structure: K-dimensional vector.
//' @param initial_mu initial value of mu, data structure: K-dimensional vector.
//' @param initial_sigma initial value of sigma, data structure: K-dimensional vector.
//' @return a list including
//' iterations: the number of iterations, data structure: int.
//' tau: the estimate of the proportion of each cluster, data structure: K-dimensional vector.
//' mu: the estimate of expectation in each cluster, data structure: K-dimensional vector.
//' sigma: the estimate of covariance matrix in each cluster, data structure: K-dimensional vector.
//' @examples
//' \dontrun{
//' #parameters of the samples
//' n <- 300; p <- 1; K <- 3
//' mu <- c(-2, 0, 3)
//' sigma <- 1
//' #generate the random samples
//' X <- numeric(n)
//' set.seed(1)
//' select <- sample(1:3, size = n, replace = TRUE, prob = c(0.3,0.2,0.5))
//' X[select==1] <- rnorm(sum(select==1), mean = mu[1], sd = sqrt(sigma))
//' X[select==2] <- rnorm(sum(select==2), mean = mu[2], sd = sqrt(sigma))
//' X[select==3] <- rnorm(sum(select==3), mean = mu[3], sd = sqrt(sigma))
//' #initialization
//' initial_tau <- rep(1/K, K)
//' initial_mu <- c(-1, 0, 1)
//' initial_sigma <- rep(var(X),K)
//' #cluster analysis
//' mixmodel_Rcpp(X, K, 1e5, initial_tau, initial_mu, initial_sigma)
//' }
//' @export
// [[Rcpp::export]]
List mixmodel_Rcpp(NumericVector X, int K, int times, NumericVector initial_tau, NumericVector initial_mu, NumericVector initial_sigma){
  //X: the random samples generated from a Gaussian mixture model, data structure: matrix, nrow = sample size, ncol = the dimension of random vectors
  //K: the number of mixture components in the population, data structure: int
  //times: the upper bound on the number of iterations, data structure: int
  int n = X.size(); //n: sample size, p: the dimension of random vectors
  NumericMatrix tau(K, times+1), mu(K, times+1), sigma(K, times+1);//the chains from iteration
  NumericMatrix z(n, K); NumericVector zk(n), z_column(n), nhat(K);
  double rate = 1e-8; //Convergence is checked when the change of elements in the sequence is small enough
  
  //initialization
  tau.column(0) = initial_tau;
  mu.column(0) = initial_mu;
  sigma.column(0) = initial_sigma;
  //  print(List::create(Named("tau0") = tau(_,0), Named("mu0") = mu(_,0), Named("sigma0") = sigma(_,0)));
  
  //EM iteration
  int t;//counternumber
  for(int i=0; i<times; i++){
    for(int k=0; k<K; k++){
      //zk = Rcpp::apply(X, 1, function(x) Rcpp::dnorm(x, mean = mu(k,i), sigma = sigma(k,i)));
      zk = Rcpp::dnorm(X, mu(k,i), sqrt(sigma(k,i)));//n-dimension
      //      print(List::create(Named("times") = t, Named("zk")=zk));
      z.column(k) = tau(k,i)*zk;//n-dimension
    }
    //    print(List::create(Named("times") = t, Named("z")=z));
    //z_column = Rcpp::apply(z, 1, sum);//n-dimension
    //z = z/z_column;//n*K-dimension
    //nhat = Rcpp::apply(z, 2, sum);//K-dimension
    //taunew = nhat/n;
    for(int j=0; j<n; j++){
      z_column[j] = sum(z.row(j));
    }
    //    print(List::create(Named("times") = t, Named("taupast") = taupast, Named("mupast") = mupast, Named("sigmapast") = sigmapast));
    for(int k=0; k<K; k++){
      z.column(k) = z.column(k)/z_column;//n-dimension
      nhat[k] = sum(z.column(k));
      tau(k, i+1) = nhat[k]/n;
      mu(k, i+1) = sum(z.column(k)*X)/nhat[k];//p-dimension
      sigma(k, i+1) = sum((X-mu(k, i+1))*z.column(k)*(X-mu(k, i+1)))/nhat[k];//(p*n)*(n*n)*(n*p)=p*p-dimension
    }
    //    print(List::create(Named("times") = t, Named("taupast") = tau(_,t-1), Named("mupast") = mu(_,t-1), Named("sigmapast") = sigma(_,t-1)));
    //    print(List::create(Named("times") = t, Named("z") = z, Named("z_column") = z_column, Named("nhat") = nhat, Named("tau") = taunew, Named("mu") = munew, Named("sigma") = sigmanew));
    //    print(List::create(Named("times") = t, Named("taunew") = tau(_,t), Named("munew") = mu(_,t), Named("sigmanew") = sigma(_,i)));
    //    print(List::create(Named("rate of tau") = abs(taunew-taupast), Named("rate of mu") = abs(munew-mupast), Named("rate of sigma") = abs(sigmanew-sigmapast) ));
    if( ((sum(abs(tau(_,i+1)-tau(_,i)))/sum(abs(tau(_,i))))<rate) & ((sum(abs(mu(_,i+1)-mu(_,i)))/sum(abs(mu(_,i))))<rate) & ((sum(abs(sigma(_,i+1)-sigma(_,i)))/sum(abs(sigma(_,i))))<rate) ){
      //print( List::create(Named("The sequence converges after the number of iterations = ")=i+1));
      t = i+1;
      break;
    }
    //    print(List::create(Named("times") = t, Named("taupast") = taupast, Named("mupast") = mupast, Named("sigmapast") = sigmapast));
  }
  return(List::create(Named("iterations") = t, Named("tau") = tau(_,t), Named("mu") = mu(_,t), Named("sigma") = sigma(_,t)));
}
