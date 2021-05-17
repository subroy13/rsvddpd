# Simulation functions for Package rsvddpd
# Authors:
#       - Subhrajyoty Roy
#       - Ayanendranath Basu
#       - Abhik Ghosh


#' Add outlier to matrix
#'
#' \code{AddOutlier} returns a matrix with outliers randomly added to a matrix
#'  given certain proportion of contamination
#' 
#' @param X \code{matrix}, to which outliers are added
#' @param proportion \code{numeric}, proportion of elements, rows or columns to be contaminated. 
#' Must be between 0 and 1.
#' @param value \code{numeric}, the outlying value to be used for contamination
#' @param seed \code{numeric}, a seed to reproduce the randomization behaviour
#' @param method \code{character}, must be one of the following:
#' \itemize{
#' \item \code{"element"} - For contaminating at random positions of the matrix
#' \item \code{"row"} - For contaminating an entire row of the matrix
#' \item \code{"col"} - For contaminating an entire column of the matrix
#' }
#' 
#' @return A \code{matrix} with elements / rows / columns contaminated.
#' @note Due to randomization, it is possible that the none of the entries of the matrix 
#' become contaminated. In that case, it is recommended to use different seed value.
#' @export
#' @examples 
#' X = matrix(1:20, nrow = 4, ncol = 5)
#' AddOutlier(X, 0.5, 10, seed = 1234)
AddOutlier <- function(X, proportion, value, seed = NULL, method = "element") {
    
    # Consistency check of the arguments
    if (length(dim(X)) != 2) {
        stop("X must be a matrix")
    }
    
    if (!method %in% c("element", "row", "col")) {
        stop("method must be one of element, row or col")
    }
    
    if ((proportion < 0) | (proportion > 1)) {
        stop("proportion must be between 0 and 1")
    }
    
    if (!is.finite(value)) {
        stop("value must be a finite and non-missing quantity")
    }
    
    
    # Setup a seed for reproducibility
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    if (method == "element") {
        indicator <- stats::rbinom(nrow(X) * ncol(X), size = 1, prob = proportion)
    } else if (method == "row") {
        indicator <- rep(stats::rbinom(nrow(X), size = 1, prob = proportion), ncol(X))
    } else {
        indicator <- rep(stats::rbinom(ncol(X), size = 1, prob = proportion), each = nrow(X))
    }
    
    return(X + indicator * value)
}



#' Simulate SVD and measure performances of various algorithms
#' 
#' \code{simSVD} simulates various models for the errors in the data matrix, and summarize
#' performance of a singular value decomposition algorithm under presence or absence of 
#' outlying data introduced through various outlying schemes, using Monte Carlo approach.
#' 
#' @param trueSVD \code{list}, containing three different named components.
#' \itemize{
#' \item d - a \code{vector} containing the singular values.
#' \item u - a \code{matrix} with left singular vectors, each column being a singular vector.
#' \item v - a \code{matrix} with right singular vectors, each column being a singular vector.
#' }
#' @param svdfun \code{function} which takes a \code{numeric} matrix as first argument and
#' returns singular value decomposition of it as a \code{list}, with three components
#'  d, u and v as indicated before.
#' @param B \code{numeric}, denoting the number of Monte Carlo simulation.
#' @param seed \code{numeric}, a seed value used for reproducibility.
#' @param dist \code{character} string, denoting the distribution from which errors will be generated. 
#' It must be equal to one of the following: \code{\link[stats:rnorm]{normal}}, \code{\link[stats:rcauchy]{cauchy}}, 
#' \code{\link[stats:rexp]{exp}}, \code{\link[stats:rlogis]{logis}}, \code{\link[stats:rlnorm]{lognormal}}
#' @param tau \code{numeric}, a value between 0 and 1, see details for more.
#' @param outlier \code{logical}, if \code{TRUE}, simulates the situation by adding outliers.
#' @param out_method \code{character}, the method to add outliers. Must be one of "element", "row" or "col". See \link{AddOutlier} for details.
#' @param out_value \code{numeric}, the outlying observation. See \link{AddOutlier} for details.
#' @param out_prop a \code{numeric}, between 0 and 1 denoting the proportion of contamination. See \link{AddOutlier} for details.
#' @param return_details \code{logical}, whether to return detailed results for each Monte Carlo simulation. See value for details.
#' @param ... extra arguments to be passed to \code{svdfun} function.
#' 
#' @return Based on whether \code{return_details} is \code{TRUE} or \code{FALSE}, returns a list with two or one components.
#' \itemize{
#' \item Simulations :
#' \itemize{
#' \item Lambda - A \code{matrix} containing obtained singular values from all Monte Carlo Simulations.
#' \item Left - A \code{matrix} containing the dissimilarities between left singular vectors of true SVD and obtained SVD.
#' \item Right - A \code{matrix} containing the dissimilarities between right singular vectors of true SVD and obtained SVD.
#' }
#' \item Summary :
#' \itemize{
#' \item Bias - A \code{numeric vector} showing biases of the singular vectors obtained by \code{svdfun} algorithm.
#' \item MSE - A \code{numeric vector} showing MSE of the singular vectors obtained by \code{svdfun} algorithm.
#' \item Variance - A \code{numeric vector} showing variances of the singular vectors obtained by \code{svdfun} algorithm.
#' \item Left - A \code{numeric vector} showing average dissimilarities between true and estimated left singular vectors.
#' \item Right - A \code{numeric vector} showing average dissimilarities between true and estimated right singular vectors.
#' }
#' }
#' If \code{return_details} is \code{FALSE}, only Summary component of the larger list is returned.
#' @export
simSVD <- function(trueSVD, svdfun, B = 100, seed = NULL, dist = "normal", tau = 0.95,
                    outlier = FALSE, out_method = "element", out_value = 10, out_prop = 0.1, 
                    return_details = FALSE, ...) {
    
    # Consistency check of the parameters
    if (is.null(trueSVD$d) | is.null(trueSVD$u) | is.null(trueSVD$v)) {
        stop("trueSVD must be a list with three components, d, u and v representing 
             the singular values, left and right singular vectors")
    } else if ( (length(dim(trueSVD$u))!= 2) | (length(dim(trueSVD$v))!=2) ) {
        stop("u and v components of trueSVD must be matrices")
    } else if ( (ncol(trueSVD$u) != length(trueSVD$d)) | (ncol(trueSVD$v) != length(trueSVD$d)) ) {
        stop("The dimensions of different components of trueSVD do not match")
    } else {
        X <- trueSVD$u %*% diag(trueSVD$d) %*% t(trueSVD$v)
    }
    
    if (B < 0) {
        stop("B must be greater than 0")
    }
    
    if (!dist %in% c("normal", "cauchy", "exp", "logis", "lognorm")) {
        stop("Only normal, cauchy, exponential, logistic and lognormal distribution of errors are supported.")
    }
    
    # Start initializing the memory to store results
    dims <- dim(X)
    lambdas <- matrix(NA, nrow = B, ncol = length(trueSVD$d))
    left_score <- matrix(NA, nrow = B, ncol = ncol(trueSVD$u))
    right_score <- matrix(NA, nrow = B, ncol = ncol(trueSVD$v))
    
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    error_seeds <- sample(1:1e6, size = B)
    if (outlier) {
        out_seeds <- sample(1:1e6, size = B)
    }
    
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (i in 1:B) {
        set.seed(error_seeds[i])
        
        if (dist == "normal") {
            errors <- matrix(stats::rnorm(dims[1] * dims[2]), nrow = dims[1], ncol = dims[2])
        } else if (dist == "cauchy") {
            errors <- matrix(stats::rcauchy(dims[1] * dims[2]), nrow = dims[1], ncol = dims[2])
        } else if (dist == "exp") {
            errors <- matrix(stats::rexp(dims[1] * dims[2]), nrow = dims[1], ncol = dims[2])
        } else if (dist == "logis") {
            errors <- matrix(stats::rlogis(dims[1] * dims[2]), nrow = dims[1], ncol = dims[2])
        } else {
            errors <- matrix(stats::rlnorm(dims[1] * dims[2]), nrow = dims[1], ncol = dims[2])
        }
        
        tmpX <- X + errors
        
        if (outlier) {
           tmpX <- rsvddpd::AddOutlier(X, proportion = out_prop, value = out_value, 
                                       method = out_method, seed = out_seeds[i]) 
        }
        
        tryCatch({
            y <- svdfun(tmpX, ...)
            
            lambdas[i, ] <- y$d[1:length(trueSVD$d)]
            left_score[i, ] <- 1 - abs( colSums(y$u[, 1:ncol(trueSVD$u)] * trueSVD$u) )
            right_score[i, ] <- 1 - abs( colSums(y$v[, 1:ncol(trueSVD$v)] * trueSVD$v) )
        },
        warning = function(cond) {
        },
        error = function(cond) {
        })
        
        utils::setTxtProgressBar(pb, value = i)
    }
    
    close(pb)  # close the connection to the progress bar
    
    
    # create summary measures
    non_na_index <- stats::complete.cases(lambdas)
    lambdas <- lambdas[non_na_index, ]
    left_score <- left_score[non_na_index, ]
    right_score <- right_score[non_na_index, ]
    
    # since some of the simulations might not converge, we use robust estimators to estimate bias and MSE
    if (!is.null(seed)) {
        set.seed(seed)
    } 
    b <- MASS::cov.rob(lambdas, quantile.used = tau * nrow(lambdas))  
    
    good_index <- b$best
    
    bias <- colMeans(lambdas[good_index, ], na.rm = T) - trueSVD$d
    mse <- colMeans( (lambdas[good_index, ] - matrix(trueSVD$d, nrow = length(good_index), ncol = ncol(lambdas), byrow = TRUE))^2, na.rm = T)
    variance <- mse - bias^2
    
    left_score_mean <- colMeans(left_score[good_index, ], na.rm = T)
    right_score_mean <- colMeans(right_score[good_index, ], na.rm = T)
    
    
    if (return_details) {
        return(list("Simulations" = list("Lambda" = lambdas, "Left" = left_score, "Right" = right_score), 
                    "Summary" = list("Bias" = bias, "MSE" = mse, "Variance" = variance, "Left" = left_score_mean, "Right" = right_score_mean) ))
    } else {
        return(list("Bias" = bias, "MSE" = mse, "Variance" = variance, "Left" = left_score_mean, "Right" = right_score_mean) )
    }
    
}
























