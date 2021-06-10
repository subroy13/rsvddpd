// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.rSVDdpd_cpp)]]
Rcpp::List rSVDdpd_cpp(arma::mat X, float alpha, arma::mat A, arma::mat B, int nd = NA_INTEGER, 
                double tol = 1e-4, double eps = 1e-4, int maxiter = 100) {
    
    
    // Determine the approximate rank of the matrix
    int rank, max_rank;
    if (X.n_rows < X.n_cols) {
        rank = X.n_rows;
    } else {
        rank = X.n_cols;
    }
    
    // Do some checking of passed arguments
    if ( Rcpp::traits::is_na<INTSXP>(nd)  ) {
       max_rank = rank; 
    } else if ((nd <= rank) && (nd > 0) ) {
        max_rank = nd;
    } else {
        Rcpp::stop("nd must be between 0 and minimum of the number of rows and colums of X");
    }
    
    if (!X.is_finite()) {
        Rcpp::stop("X matrix must not contain NA or Inf values");
    }
    
    if ((alpha < 0) || (alpha > 1)) {
        Rcpp::warning("Robustness parameter alpha should be between 0 and 1. Recommended value is 0.3");
    }
    
    if ((tol < 0) || (eps < 0)) {
        Rcpp::stop("Tolerance levels tol and eps must be small positive real numbers");
    }
    
    
    
    
    // Do some normalization to avoid overflow or underflow
    arma::vec P = {0.25, 0.75}; 
    arma::vec scale_factors = arma::quantile(arma::vectorise(X), P);
    X = X / (scale_factors(1) - scale_factors(0));
    
    arma::vec Lambda = arma::zeros(rank);
    arma::vec maxit_reached = arma::zeros(max_rank);

    // Iteration through the singular values
    for (int r = 0; r < max_rank; r++) {
        double curr_norm = arma::norm(X, 2);
        if (curr_norm > tol) {
            int n_iter = 0;
            bool is_converged = false;
            double sigma = 1;
            
            // some initially declared variables to be used throughout
            arma::mat weight_mat = arma::zeros(X.n_rows, X.n_cols);
            arma::mat numerator = arma::zeros(X.n_rows, X.n_cols);
            arma::mat denominator = arma::zeros(X.n_rows, X.n_cols);
            
            
            // Perform iteration between left and right singular vectors
            do {
                arma::vec curr_a = A.col(r);
                arma::vec curr_b = B.col(r);
                double curr_lambda = arma::as_scalar(Lambda(r));
                
                // FIX RIGHT, OPTIMIZE LEFT
                // Do fixed point iteration to get best left singular vector
                int left_iter = 0;
                bool fixed = false;  // holds convergence for fixed point iteration
                
                arma::vec c = arma::ones(X.n_rows);
                
                do {
                    arma::vec curr_c = c;
                    weight_mat = X - curr_c * curr_b.t();
                    weight_mat = arma::exp(-0.5 * alpha * arma::square(weight_mat) / sigma);
                    numerator = weight_mat % X;
                    numerator.each_row() %= curr_b.t();
                    denominator = weight_mat;
                    denominator.each_row() %= (arma::square(curr_b).t());
                    c = arma::sum(numerator, 1) / arma::sum(denominator, 1);
                    left_iter += 1;
                    
                    // check if fixed point criterion is met
                    fixed = (arma::norm((c - curr_c), 2) < eps) || (left_iter > maxiter);
                    if (left_iter > maxiter) {
                        maxit_reached(r) = 1;
                    }
                    
                    // allow user interruption
                    Rcpp::checkUserInterrupt();
                    
                } while (!fixed);
                
                // Need to apply Gram Schmidt
                if (r > 0) {
                    arma::mat temp = A.cols(0, (r-1));
                    c -= ((temp * temp.t()) * c);
                }
                
                curr_lambda = arma::norm(c, 2);
                curr_a = c / curr_lambda;
                
                
                // FIX LEFT, OPTIMIZE RIGHT
                // Do fixed point iteration to get best left singular vector
                int right_iter = 0;
                fixed = false;  // holds convergence for fixed point iteration
                
                arma::vec d = arma::ones(X.n_cols);
                
                do {
                    arma::vec curr_d = d;
                    weight_mat = X - curr_a * curr_d.t();
                    weight_mat = arma::trunc_exp(-0.5 * alpha * arma::square(weight_mat) / sigma);
                    numerator = weight_mat % X;
                    numerator.each_col() %= curr_a;
                    denominator = weight_mat;
                    denominator.each_col() %= (arma::square(curr_a));
                    d = arma::sum(numerator, 0).t() / arma::sum(denominator, 0).t();
                    right_iter += 1;
                    
                    // check if fixed point criterion is met
                    fixed = (arma::norm((d - curr_d), 2) < eps) || (right_iter > maxiter);
                    
                    if (right_iter > maxiter) {
                        maxit_reached(r) = 1;
                    }
                    
                    // allow user interruption
                    Rcpp::checkUserInterrupt();
            
                } while (!fixed);
                
                
                // Need to apply Gram Schmidt
                if (r > 0) {
                    arma::mat temp = B.cols(0, (r-1));
                    d -= ((temp * temp.t()) * d);
                } 
                
                curr_lambda = arma::norm(d, 2);
                curr_b = d / curr_lambda;
                
                // updating sigma square
                numerator = arma::square(X - curr_lambda * curr_a * curr_b.t());
                weight_mat = arma::trunc_exp(-0.5 * alpha * numerator);
                double denom = arma::as_scalar(arma::accu(weight_mat)) - (alpha / std::pow(alpha+1, 1.5)); 
                sigma = arma::as_scalar(arma::accu(numerator % weight_mat)) / denom;
                // do not let sigma become very small
                if (sigma < 1e-1) {
                    sigma = 1e-1;
                }
                
                n_iter += 1;
                
                // allow user interruption
                Rcpp::checkUserInterrupt();
                
                // check if final convergence criteria is met
                bool is_convl = ( std::abs(curr_lambda - arma::as_scalar(Lambda(r)) ) < eps );
                bool is_conva = ( arma::norm(curr_a - A.col(r), 2) < eps );
                bool is_convb = ( arma::norm(curr_b - B.col(r), 2) < eps );
                
                is_converged = (n_iter > maxiter) || (is_convl && is_conva && is_convb);
                if (n_iter > maxiter) {
                    maxit_reached(r) = 1;
                }
                
                // Update the values
                sigma = 1;
                Lambda(r) = curr_lambda;
                A.col(r) = curr_a;
                B.col(r) = curr_b;
                
            } while (!is_converged);
            
            // One singular value obtain, proceed to the next
            X = X - Lambda(r) * A.col(r) * B.col(r).t();
            curr_norm = arma::norm(X, 2);
            
        }
    }
    
    // Change the singular values as required
    Lambda = (scale_factors(1) - scale_factors(0)) * Lambda;
    
    Rcpp::List L = Rcpp::List::create(Rcpp::Named("d") = Lambda.subvec(0,max_rank-1), 
                                      Rcpp::Named("u") = A.cols(0, max_rank - 1),
                                      Rcpp::Named("v") = B.cols(0, max_rank - 1),
                                      Rcpp::Named("maxitreach") = maxit_reached.subvec(0,max_rank-1) );
    
    return(L);
}












