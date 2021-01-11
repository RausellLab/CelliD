#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List MCAStep1(NumericMatrix X) {
    arma::mat AM = arma::mat(X.begin(), X.rows(), X.cols(), true);
    arma::colvec rmin = arma::min(AM,1);
    arma::colvec rmax = arma::max(AM,1);
    arma::colvec range = (rmax -rmin);
    AM.each_col() -= rmin;
    AM.each_col() /= range;
    arma::mat FM = join_cols(AM, 1 - AM);
    AM.clear();
    long total = arma::accu(FM);
    arma::rowvec colsum = arma::sum(FM,0);
    arma::colvec rowsum = arma::sum(FM,1);
    FM.each_row() /= sqrt(colsum);
    FM.each_col() /= sqrt(rowsum);
    arma::colvec Dc = 1/(sqrt(rowsum/total));
    return List::create(Named("Z") = wrap(FM),
                        Named("Dc") = wrap(Dc));
}

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List MCAStep2(NumericMatrix Z, NumericMatrix V, NumericVector Dc) {
    arma::mat AV = arma::mat(V.begin(), V.rows(), V.cols(),  true);
    arma::mat AZ = arma::mat(Z.begin(), Z.rows(), Z.cols(),  true);
    arma::colvec ADc = arma::colvec(Dc);
    arma::mat FeaturesCoordinates =  AZ * AV;
    int AZcol = AZ.n_cols;
    AZ.clear();
    FeaturesCoordinates.each_col() %= ADc;
    ADc.clear();
    return List::create(Named("cellsCoordinates") = wrap(std::sqrt(AZcol) * AV),
                        Named("featuresCoordinates") = wrap(FeaturesCoordinates.head_rows(FeaturesCoordinates.n_rows/2)));
}