#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix mcmc_Gibbs_Rcpp(double a, double b, double n, NumericVector initial, int N){
  //a, b, n: 二元分布参数
  //initial: 初始值
  //N: 需要生成的样本链长度
  
  NumericMatrix X(2, N); //样本链
  X(0, 0) = initial(0); //赋初值
  X(1, 0) = initial(1);
  
  for (int t=1; t < N; t++){
    //从Y_{t-1}得到X_{t}，(X_t|Y_{t-1})服从Binomial(n,Y_{t-1})分布
    X(0, t) = Rcpp::rbinom(1, n, X(1, t-1))[0];
    //从X_{t}得到Y_{t}，(Y_t|X_t)服从Beta(X_t+a,n-X_t+b)分布
    X(1, t) = Rcpp::rbeta(1, X(0, t) + a, n - X(0, t) + b)[0];
  }
  return(X);
}