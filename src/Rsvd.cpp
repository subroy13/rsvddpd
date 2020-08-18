// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' Robust Singular Value Decomposition using Density Power Divergence
//'
//' \code{Rsvd} returns the singular value decomposition of a matrix with robust
//' singular values in presence of outliers
//' 
//' @param X \code{matrix}, whose singular value decomposition is required
//' @param alpha \code{numeric}, robustness parameter between 0 and 1. See details for more.
//' @param nd \code{integer}, must be lower than \code{nrow(X)} and \code{ncol(X)} both. If 
//' NA, defaults to \code{min(nrow(X), ncol(X))} 
//' @param tol \code{numeric}, a tolerance level. If the residual matrix has lower 
//' norm than this, then subsequent singular values will be taken as 0.
//' @param eps \code{numeric}, a tolerance level for the convergence of singular 
//' vectors. If in subsequent iterations the singular vectors do not change its 
//' norm beyond this, then the iteration will stop.
//' @param maxiter \code{integer}, upper limit to the maximum number of iterations.
//' 
//' @return A list containing different components of the decomposition \eqn{X = U D V'}
//' \itemize{
//' \item d - The robust singular values, namely the diagonal entries of \eqn{D}.
//' \item u - The matrix of left singular vectors \eqn{U}. Each column is a singular vector.
//' \item v - The matrix of right singular vectors \eqn{V}. Each column is a singular vector.
//' }
//' 
//' @details The usual singular value decomposition is highly prone to error in 
//' presence of outliers, since it tries to minimize the \eqn{L_2} norm of the errors
//' between the matrix \eqn{X} and its best lower rank approximation. While there is
//' considerable effort to impose robustness using \eqn{L_1} norm of the errors instead
//' of \eqn{L_2} norm, such estimation lacks efficiency. Application of density power
//' divergence bridges the gap.
//' \deqn{DPD(f|g) = \int f^{(1+\alpha)} - (1 + \frac{1}{\alpha}) \int f^{\alpha}g + \frac{1}{\alpha} \int g^{(1 + \alpha)} }
//' The parameter \code{alpha} must be between 0 and 1, lower \code{alpha} means less robustness
//' but more efficiency in estimation, while higher \code{alpha} means high robustness but 
//' less efficiency in estimation. The recommended value of \code{alpha} is 0.3.
//' The function tries to obtain the best rank one approximation of a matrix by minimizing 
//' this density power divergence of the true errors with that of a normal distribution centered
//' at the origin.  
//' 
//' @seealso \code{\link{svd}}
//' @export
//' @examples
//' X = matrix(1:20, nrow = 4, ncol = 5)
//' Rsvd(X, alpha = 0.3)
// [[Rcpp::export]]
Rcpp::List Rsvd(arma::mat X, float alpha, int nd = NA_INTEGER, 
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
        Rcpp::stop("Robustness parameter alpha must be between 0 and 1. Recommended value is 0.3");
    }
    
    if ((tol < 0) || (eps < 0)) {
        Rcpp::stop("Tolerance levels tol and eps must be small positive real numbers");
    }
    
    
    
    
    // Do some normalization to avoid overflow or underflow
    arma::vec P = {0.25, 0.75}; 
    arma::vec scale_factors = arma::quantile(arma::vectorise(X), P);
    X = X / (scale_factors(1) - scale_factors(0));
    
    // Choose random matrix to initialize
    arma::mat A = arma::randu(X.n_rows, rank);
    arma::mat B = arma::randu(X.n_cols, rank);
    
    for (int i = 0; i < rank; i++) {
        A.col(i) /= arma::norm(A.col(i), 2);
        B.col(i) /= arma::norm(B.col(i), 2);
    }
    
    arma::vec Lambda = arma::zeros(rank);

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
                
                sigma = 1;
                
                n_iter += 1;
                // allow user interruption
                if (n_iter % 10 == 0) {
                    Rcpp::checkUserInterrupt();
                }
                
                // check if final convergence criteria is met
                bool is_convl = ( std::abs(curr_lambda - arma::as_scalar(Lambda(r)) ) < eps );
                bool is_conva = ( arma::norm(curr_a - A.col(r), 2) < eps );
                bool is_convb = ( arma::norm(curr_b - B.col(r), 2) < eps );
                
                is_converged = (n_iter > maxiter) || (is_convl && is_conva && is_convb);
                
                // Update the current values
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
                                      Rcpp::Named("v") = B.cols(0, max_rank - 1));
    
    return(L);
}












