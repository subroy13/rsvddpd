# rsvddpd 1.0

* Initial version of the package. Contains four functions.
    - `rSVDdpd` performs the robust svd using density power divergence. Also provides different convergence rule based warnings.
    - `AddOutlier` adds different types of outliers in a matrix.
    - `simSVD` simulates different scenarios and evaluates performance of SVD under pure model or under contamination.
    - `cv.alpha` function for finding the optimal robustness parameter.