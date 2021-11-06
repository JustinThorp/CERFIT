//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <vector>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List find_split(const arma::vec &y,const arma::vec &x,const arma::vec &trt, //arma::mat z_mat,
const NumericVector &cutpts,String method,arma::vec &propensity, double minbucket, String response_type) {
  double cutpt;
  double mingroup;
  double top;
  double bottom;
  arma::mat x_t;
  arma::mat bread;
  arma::mat meat;
  NumericVector min_z_split;
  //arma::vec x = arma::sort(x_unsort);
  //arma::uvec x_sort_index = arma::sort_index(x_unsort);
  //arma::vec y = y_unsort(x_sort_index);
  //arma::vec trt = trt_unsort(x_sort_index);
  int n_cut = cutpts.length();
  int n = y.n_rows;//cutpts.length();
  arma::vec res(n);
  arma::vec z(n);
  //arma::mat z_mat = x * cutpts.t();
  std::vector<double> score(n_cut);
  //NumericVector score(n_cut);
  //= arma::zeros<arma::vec>(n);
  NumericVector x_numeric = wrap(x);
  //NumericVector trt_numeric = wrap(trt);
  arma::vec double_trt = 2*trt;
  if (method == "RCT") {
    arma::mat x_lm (n,4);
    arma::mat beta(4,1);
    x_lm.fill(1);
    x_lm.col(2) = trt;
    double se_mat;
    for (int i = 0; i < n_cut ;i++) {
      cutpt = cutpts(i);
      //z = z_mat.row(i).t();
      z = ifelse(x_numeric <= cutpt,1.0,0.0);
      //z = z_mat.col(i);
      NumericVector test = wrap(double_trt + z);
      mingroup = min(table(test));
      min_z_split = {sum(z),sum(1-z)};
      if (min(min_z_split) < .1*y.n_rows || mingroup < minbucket/2) {
        score[i] = -1;//NA_REAL;
      } else {
        if (response_type == "binary") {
          score[i] = 5.0;//lm_rand(y,z,trt); //Need to make a glm function
        } else {
          x_lm.col(1) = z;
          x_lm.col(3) = z % trt;
          if (arma::rank(x_lm) != 4) {
            score[i] = -1;
            continue;
          }
          try {
            x_t = x_lm.t();
            bread = (x_t * x_lm).i();
            //bread = arma::inv_sympd(x_t * x_lm);
            beta = bread*x_t*y;//arma::solve(x_lm,y,solve_opts::fast);
            res = y - x_lm * beta;
            double sigma = arma::as_scalar(res.t() * res) / (n-4);
            //meat = x_t.each_row() % (res % res).t() * x_lm;
            se_mat = sigma*bread(3,3);// * meat * bread;
            score[i] = beta(3,0)*beta(3,0) / se_mat;
          } catch(...) {
            score[i] = -1;
          }
        }
      }
    }
  } else {
    arma::mat x_lm (n,5);
    arma::mat beta(5,1);
    x_lm.fill(1);
    x_lm.col(2) = trt;
    x_lm.col(4) = propensity;
    arma::mat se_mat;
    for (int i = 0; i < n_cut ;i++) {
      cutpt = cutpts(i);
      //z = z_mat.row(i).t();
      z = ifelse(x_numeric <= cutpt,1.0,0.0);
      //z(i) = 1;
      NumericVector test = wrap(double_trt + z);
      mingroup = min(table(test));
      min_z_split = {sum(z),sum(1-z)};
      if (min(min_z_split) < .1*y.n_rows || mingroup < minbucket/2) {
        score[i] = -1;//NA_REAL;
      } else {
        if (response_type == "binary") {
          x_lm.col(1) = z;
          x_lm.col(3) = z % trt;
          if (arma::rank(x_lm) != 5.0) {
            score[i] = -1;
            continue;
          }
          try {
            x_t = x_lm.t();
            bread = (x_t * x_lm).i();
            //bread = arma::inv_sympd(x_t * x_lm);
            beta = bread*x_t*y;//arma::solve(x_lm,y,solve_opts::fast);
            res = y - x_lm * beta;
            //double sigma = var(res);
            //meat = x_t.each_row() % (res % res).t() * x_lm;
            meat = x_t * arma::diagmat(res % res) * x_lm;
            se_mat = bread * meat * bread;
            top = beta(3,0) * beta(3,0);
            bottom = se_mat(3,3);
            if (bottom < .00000001) {
              score[i] = -1;
            } else {
              score[i] = top / bottom;
            }
          } catch(...) {
            score[i] = -1;
          }//lm_rand(y,z,trt); //Need to make a glm function
        } else {
          x_lm.col(1) = z;
          x_lm.col(3) = z % trt;
          if (arma::rank(x_lm) != 5.0) {
            score[i] = -1;
            continue;
          }
          x_t = x_lm.t();
          bread = (x_t * x_lm).i();
          //bread = arma::inv_sympd(x_t * x_lm);
          beta = bread*x_t*y;//arma::solve(x_lm,y,solve_opts::fast);
          res = y - x_lm * beta;
          //double sigma = var(res);
          meat = x_t.each_row() % (res % res).t() * x_lm;
          //meat = x_t * arma::diagmat(res % res) * x_lm;
          se_mat = bread * meat * bread;
          double top = arma::as_scalar(beta(3,0) * beta(3,0));
          double bottom = arma::as_scalar(se_mat(3,3));
          score[i] = top / bottom;
          //if (bottom < .00000001) {
          //  score[i] = -1;
          //} else {
          //  score[i] = top / bottom;
          //}


        }
      }
    }
  }
  NumericVector score_temp = wrap(score);
  double max_score = max(score_temp);
  int max_index = which_max(score_temp);
  double cutoff = cutpts[max_index];
  if (max_score == -1){
    max_score = NA_REAL;
    cutoff = NA_REAL;
  }
  return List::create(Named("cutoff") = cutoff, Named("stat") = max_score);
}


