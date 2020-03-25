#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// cpp_prangles : principal angles
// cpp_pdist    : distance measure for one set  of data
// cpp_pdist2   : distance measure for two sets of data
// cpp_pairdist : pairwise distances

// [[Rcpp::export]]
arma::vec cpp_prangles(arma::mat U, arma::mat V){
  // prepare
  int u = U.n_cols;
  int v = V.n_cols;
  int n = 0;
  if (u < v){
    n = u;
  } else {
    n = v;
  }
  
  double aco1 = std::acos(static_cast<float>(1.0));
  double aval = 0.0;
  arma::vec s = arma::svd(U.t()*V);
  arma::vec output(n,fill::zeros);
  for (int i=0;i<n;i++){
    if (s(i) > 1){
      aval = aco1;
    } else {
      aval = std::acos(static_cast<float>(s(i)));
    }
    output(i) = aval;
  }
  return(output);
}
// cpp_pairdist : pairwise distances
double cpp_pair_asimov(arma::vec thetas, int k){
  return(arma::accu(thetas(k-1)));
}
double cpp_pair_binetcauchy(arma::vec thetas, int k){
  double multiple = 1.0;
  double costheta = 0.0;
  for (int i=0;i<k;i++){
    costheta = std::cos(static_cast<float>(thetas(i)));
    multiple = multiple*(costheta*costheta);
  }
  return(std::sqrt((static_cast<float>(1.0 - multiple))));
}
double cpp_pair_chordal(arma::vec thetas, int k){
  double summed   = 0.0;
  double sintheta = 0.0;
  for (int i=0;i<k;i++){
    sintheta = std::sin(static_cast<float>(thetas(i)));
    summed   = summed + (sintheta*sintheta);
  }
  return(std::sqrt(static_cast<float>(summed)));
}
double cpp_pair_fubinistudy(arma::vec thetas, int k){
  double multiple = 1.0;
  for (int i=0;i<k;i++){
    multiple = multiple*std::cos(static_cast<float>(thetas(i)));
  }
  return(std::acos(static_cast<float>(multiple)));
}
double cpp_pair_martin(arma::vec thetas, int k){
  double summed = 0.0;
  for (int i=0;i<k;i++){
    summed = summed - 2.0*std::log(std::cos(static_cast<float>(thetas(i))));
  }
  return(std::sqrt(static_cast<float>(summed)));
}
double cpp_pair_procrustes(arma::vec thetas, int k){
  double summed = 0.0;
  double sinthe = 0.0;
  for (int i=0;i<k;i++){
    sinthe  = std::sin(static_cast<float>(thetas(i)/2.0));
    summed += (sinthe*sinthe);
  }
  return(2.0*std::sqrt(static_cast<float>(summed)));
}
double cpp_pair_projection(arma::vec thetas, int k){
  return(std::sin(static_cast<float>(thetas(k-1))));
}
double cpp_pair_spectral(arma::vec thetas, int k){
  return(2.0*std::sin(static_cast<float>(thetas(k-1)/2.0)));
}
// [[Rcpp::export]]
double cpp_pairdist(arma::mat U, arma::mat V, std::string dname){
  // compute principal angles
  arma::vec thetas = cpp_prangles(U, V);
  int k = thetas.n_elem;
  double output = 0.0;
  if (dname=="asimov"){
    output = cpp_pair_asimov(thetas, k);
  } else if (dname=="binetcauchy"){
    output = cpp_pair_binetcauchy(thetas, k);
  } else if (dname=="chordal"){
    output = cpp_pair_chordal(thetas, k);
  } else if (dname=="fubinistudy"){
    output = cpp_pair_fubinistudy(thetas, k);
  } else if (dname=="martin"){
    output = cpp_pair_martin(thetas, k);
  } else if (dname=="procrustes"){
    output = cpp_pair_procrustes(thetas, k);
  } else if (dname=="projection"){
    output = cpp_pair_projection(thetas, k);
  } else if (dname=="spectral"){
    output = cpp_pair_spectral(thetas, k);
  }
  return(output);
}
// cpp_pdist    : distance measure for one set  of data
// [[Rcpp::export]]
arma::mat cpp_pdist(arma::field<arma::mat> INPUT1, std::string dname){
  // parameter
  int N = INPUT1.n_elem;
  arma::mat output(N,N,fill::zeros);
  
  arma::mat U;
  arma::mat V;
  for (int i=0;i<(N-1);i++){
    U = INPUT1(i);
    for (int j=(i+1);j<N;j++){
      V = INPUT1(j);
      output(i,j) = cpp_pairdist(U, V, dname);
      output(j,i) = output(i,j);
      V.reset();
    }
    U.reset();
  }
  return(output);
}
// cpp_pdist2   : distance measure for two sets of data
// [[Rcpp::export]]
arma::mat cpp_pdist2(arma::field<arma::mat> INPUT1, arma::field<arma::mat> INPUT2, std::string dname){
  // parameter
  int N = INPUT1.n_elem;
  int M = INPUT2.n_elem;
  arma::mat output(N,M,fill::zeros);
  
  arma::mat U;
  arma::mat V;
  for (int n=0;n<N;n++){
    U = INPUT1(n);
    for (int m=0;m<M;m++){
      V = INPUT2(m);
      output(n,m) = cpp_pairdist(U,V,dname);
      V.reset();
    }
    U.reset();
  }
  return(output);
}
