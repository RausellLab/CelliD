#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]
NumericMatrix fastOrder(NumericMatrix Ar, NumericMatrix Br) {
    int m = Ar.nrow(), 
        n = Br.nrow(),
        k = Ar.ncol();
    arma::mat A = arma::mat(Ar.begin(), m, k, false);
    arma::mat B = arma::mat(Br.begin(), n, k, false); 
    arma::colvec An =  sum(square(A),1);
    arma::colvec Bn =  sum(square(B),1);
    arma::mat C = -2 * (A * B.t());
    C.each_col() += An;
    C.each_row() += Bn.t();
    C = sqrt(C);
    arma::umat ordMat = arma::umat(C.n_rows, C.n_cols);
    for(int a = 0; a < C.n_cols; a = a + 1 ) {
        ordMat.col(a) = arma::sort_index(C.col(a), "ascend") + 1;
    }
    return wrap(ordMat); 
}