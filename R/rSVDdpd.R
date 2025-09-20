#' Robust Singular Value Decomposition using Density Power Divergence
#'
#' \code{rSVDdpd} returns the singular value decomposition of a matrix with robust
#' singular values in presence of outliers
#' 
#' @param X \code{matrix}, whose singular value decomposition is required
#' @param alpha \code{numeric}, robustness parameter between 0 and 1. See details for more.
#' @param nd \code{integer}, must be lower than \code{nrow(X)} and \code{ncol(X)} both. If NA, determined by \code{rank.rSVDdpd(X, alpha, maxrank)}
#' @param maxrank \code{integer}, maximum rank to be considered if \code{nd} is not specified. If NA, defaults to \code{min(nrow(X), ncol(X))} 
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
#' @references S. Roy, A. Basu and A. Ghosh (2021), A New Robust Scalable Singular Value Decomposition Algorithm for Video Surveillance Background Modelling
#' \url{https://arxiv.org/abs/2109.10680}
#' @seealso \code{\link{rank.rSVDdpd}}, \code{\link{svd}}
#' @export
#' @examples
#' X = matrix(1:20, nrow = 4, ncol = 5)
#' rSVDdpd(X, alpha = 0.3)

rSVDdpd <- function(X, alpha, nd = NA, maxrank = NA, tol = 1e-4, eps = 1e-4, maxiter = 100L, initu = NULL, initv = NULL) {
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

    if (is.na(nd)) {
        if (is.na(maxrank)) {
            maxrank <- r
        }
        message("Estimating rank using rank.rSVDdpd (DICMR criterion).")
        nd <- rank.rSVDdpd(X, alpha = alpha, maxrank = maxrank)["DICMR"]
    }
    
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
#' @references S. Roy, A. Basu and A. Ghosh (2021), A New Robust Scalable Singular Value Decomposition Algorithm for Video Surveillance Background Modelling
#' \url{https://arxiv.org/abs/2109.10680}
#' 
#' @return A list containing
#' \itemize{
#' \item The choices of the robust parameters.
#' \item Corresponding cross validation score.
#' \item Best choice of the robustness parameter.
#' }
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


#' Rank Estimation for Robust Singular Value Decomposition
#'
#' \code{rank.rSVDdpd} estimates the optimal rank of a given matrix under
#' robust SVD using Density Power Divergence (DPD) criteria.
#'
#' The function computes three penalized criteria for rank determination:
#' \itemize{
#'   \item \strong{DIC} — Divergence Information Criterion.
#'   \item \strong{RCC} — Robust Cross-Validation Criterion.
#'   \item \strong{DICMR} — Modified Divergence Information Criterion with Matrix Rank
#'     penalty (recommended).
#' }
#'
#' @param X \code{matrix}, the data matrix for which robust rank estimation is required.
#' @param alpha \code{numeric}, robustness parameter between 0 and 1 (default \code{0.5}).
#'   Controls the trade-off between robustness and efficiency in the DPD measure.
#' @param maxrank \code{integer}, maximum rank to be considered. Defaults to
#'   \code{min(dim(X))}.
#'
#' @return A named integer vector of length 3, giving the estimated ranks according
#'   to each criterion:
#'   \itemize{
#'     \item \code{DIC} — estimated rank from DIC.
#'     \item \code{RCC} — estimated rank from RCC.
#'     \item \code{DICMR} — estimated rank from DICMR (recommended).
#'   }
#'
#' @details The function computes a full robust SVD (up to \code{maxrank}) using
#'   \code{\link{rSVDdpd}}. It then evaluates the DPD divergence at different
#'   candidate ranks and applies penalty adjustments for model complexity.
#'   The final estimated rank minimizes the penalized criterion.
#'
#' @seealso \code{\link{rSVDdpd}}, \code{\link{svd}}
#'
#' @export
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' rank.rSVDdpd(X, alpha = 0.3, maxrank = 5)
rank.rSVDdpd <- function(X, alpha = 0.5, maxrank = NULL) {
	if (is.null(maxrank)) {
        maxrank <- min(dim(X))
    }
    n <- nrow(X)
    p <- ncol(X)
    X.rsvd <- rSVDdpd(X, alpha = alpha, nd = maxrank) # maximum rank decomposition

    # compute the estimates of sigma
    sigmahat <- numeric(maxrank)
    dpd.vals <- numeric(maxrank)

    penalty.dic <- numeric(maxrank)
    penalty.rcc <- numeric(maxrank)
    penalty.dicmr <- numeric(maxrank)

    for (rank in 1:maxrank) {
        E <- X - X.rsvd$u[, 1:rank, drop = FALSE] %*% diag(X.rsvd$d[1:rank], rank, rank) %*% t(X.rsvd$v[, 1:rank, drop = FALSE])
        sigmahat[rank] <- median(abs(E)) * 1.4826   # consistent estimate of sigma

        # compute the DPD values
        C <- (2*pi)^(-alpha/2) * sigmahat[rank]^(-alpha)    
        dpd.vals[rank] <- C * ((1 + alpha)^(-1/2) - (1+1/alpha) * mean(exp(-alpha * E^2/sigmahat[rank]^2))  )

        penalty.dic[rank] <- rank * (alpha + 1) * (2*pi)^(-alpha/2) * ((1 + alpha)/(1+2*alpha))^(rank/2)
        penalty.rcc[rank] <- rank * log(n * p)/(n * p)
        penalty.dicmr[rank] <- C * ((1 + alpha)/(1+2*alpha))^(3/2) * rank * (1/n + 1/p)
    }

    # Return estimated ranks by minimizing penalized criteria
    return(c(
        "DIC"   = which.min(dpd.vals + penalty.dic),
        "RCC"   = which.min(dpd.vals + penalty.rcc),
        "DICMR" = which.min(dpd.vals + penalty.dicmr)  # recommended
    ))
}
