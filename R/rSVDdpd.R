#' Robust Singular Value Decomposition using Density Power Divergence
#'
#' \code{rSVDdpd} returns the singular value decomposition of a matrix with robust
#' singular values in presence of outliers
#' 
#' @param X \code{matrix}, whose singular value decomposition is required
#' @param alpha \code{numeric}, robustness parameter between 0 and 1. See details for more.
#' @param nd \code{integer}, must be lower than \code{nrow(X)} and \code{ncol(X)} both. If 
#' NA, defaults to \code{min(nrow(X), ncol(X))} 
#' @param tol \code{numeric}, a tolerance level. If the residual matrix has lower 
#' norm than this, then subsequent singular values will be taken as 0.
#' @param eps \code{numeric}, a tolerance level for the convergence of singular 
#' vectors. If in subsequent iterations the singular vectors do not change its 
#' norm beyond this, then the iteration will stop.
#' @param maxiter \code{integer}, upper limit to the maximum number of iterations.
#' @param initu \code{matrix}, initializing vectors for left singular values. Must be of dimension \code{nrow(X)} \eqn{\times} \code{min(nrow(X), ncol(X))}. If \code{NULL}, defaults to random initialization.
#' @param initv \code{matrix}, initializing vectors for right singular values. Must be of dimension \code{ncol(X)} \eqn{\times} \code{min(nrow(X), ncol(X))}. If \code{NULL}, defaults to random initialization.
#' 
#' @return A list containing different components of the decomposition \eqn{X = U D V'}
#' \itemize{
#' \item d - The robust singular values, namely the diagonal entries of \eqn{D}.
#' \item u - The matrix of left singular vectors \eqn{U}. Each column is a singular vector.
#' \item v - The matrix of right singular vectors \eqn{V}. Each column is a singular vector.
#' }
#' 
#' @details The usual singular value decomposition is highly prone to error in 
#' presence of outliers, since it tries to minimize the \eqn{L_2} norm of the errors
#' between the matrix \eqn{X} and its best lower rank approximation. While there is
#' considerable effort to impose robustness using \eqn{L_1} norm of the errors instead
#' of \eqn{L_2} norm, such estimation lacks efficiency. Application of density power
#' divergence bridges the gap.
#' \deqn{DPD(f|g) = \int f^{(1+\alpha)} - (1 + \frac{1}{\alpha}) \int f^{\alpha}g + \frac{1}{\alpha} \int g^{(1 + \alpha)} }
#' The parameter \code{alpha} should be between 0 and 1, if not, then a warning is shown.
#' Lower \code{alpha} means less robustness
#' but more efficiency in estimation, while higher \code{alpha} means high robustness but 
#' less efficiency in estimation. The recommended value of \code{alpha} is 0.3.
#' The function tries to obtain the best rank one approximation of a matrix by minimizing 
#' this density power divergence of the true errors with that of a normal distribution centered
#' at the origin.  
#' 
#' @seealso \code{\link{svd}}
#' @export
#' @examples
#' X = matrix(1:20, nrow = 4, ncol = 5)
#' rSVDdpd(X, alpha = 0.3)

rSVDdpd <- function(X, alpha, nd = NA, tol = 1e-4, eps = 1e-4, maxiter = 100L, initu = NULL, initv = NULL) {
    # checking of initial values
    n <- nrow(X)
    p <- ncol(X)
    r <- min(n, p)
    if (is.null(initu)) {
        initu <- matrix(stats::runif(n * r), nrow = n, ncol = r)
    }
    if (is.null(initv)) {
        initv <- matrix(stats::runif(p * r), nrow = p, ncol = r)
    }
    if (any(dim(initu) != c(n, r)) | any(dim(initv) != c(p, r))) {
        stop("The dimensionality of initial values does not match!")
    }
    normA <- matrixStats::colSums2(initu^2)
    normB <- matrixStats::colSums2(initv^2)
    initu <- initu / matrix(normA, nrow = n, ncol = r, byrow = T)
    initv <- initv / matrix(normB, nrow = p, ncol = r, byrow = T)
    
    out <- .rSVDdpd_cpp(X, alpha, initu, initv, nd, tol, eps, maxiter)
    
    # check convergence issues
    if (any(as.logical(out$maxitreach))) {
        warning(paste("Maximum iteration reached for obtaining singular values ", paste(which(out$maxitreach==1), collapse = ", "), ". Consider increasing the maximum number of iteration." ))
    }
    if (is.unsorted(rev(out$d))) {
        warning("Singular values are not ordered. Indicates convergence issues. Please retry with another initialization.")
    }
    if (out$d[1] > sqrt(norm(X, "1") * norm(X, "I")) ) {
        warning("First singular value is beyond theoretical bounds. Indicates convergence issues. Please retry with another initialization.")
    }
    
    return(list("d" = out$d, "u" = out$u, "v" = out$v))
}


#' Calculate optimal robustness parameter
#'
#' \code{cv.alpha} returns the optimal robustness parameter
#' @param X \code{matrix}, whose singular value decomposition is required
#' @param alphas \code{numeric vector}, vector of robustness parameters to try.
#'
#' @export
cv.alpha <- function(X, alphas = 10) {
    n <- nrow(X)
    p <- ncol(X)
    if (length(alphas) == 1) {
        alphas <- seq(0, 1, length.out = alphas)
    }
    n_alpha <- length(alphas)
    cvscore <- numeric(n_alpha)
    
    out1 <- suppressWarnings(rSVDdpd(X, alpha = 1, nd = 1))
    trueu <- out1$d[1] * out1$u[, 1]
    truev <- out1$d[1] * out1$v[, 1]
    resids <- X - out1$d[1] * out1$u[, 1] %*% t(out1$v[, 1])
    sigma_hat <- stats::median(abs(resids))/0.675
    
    pb <- utils::txtProgressBar(max = n_alpha, style = 3)
    for (i in 1:n_alpha) {
        out <- suppressWarnings(rSVDdpd(X, alpha = alphas[[i]], nd = 1))
        biasu <- out$d * out$u[,1] - trueu
        biasv <- out$d * out$v[, 1] - truev
        vars <- sigma_hat^2 * (1 + alphas[[i]]^2/(1 + 2*alphas[[i]]))^1.5
        cvscore[i] <- sum(biasu^2) + sum(biasv^2) + (n + p)*vars
        utils::setTxtProgressBar(pb, value = i)
    }
    close(pb)
    
    return(list("alphas" = alphas, "cvscore" = cvscore, "opt.alpha" = alphas[which.min(cvscore)]))
}








