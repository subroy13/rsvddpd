# rsvddpd 1.0

* Initial version of the package. Contains three functions.
    - `rSVDdpd` performs the robust svd using density power divergence.
    - `AddOutlier` adds different types of outliers in a matrix.
    - `simSVD` simulates different scenarios and evaluates performance of SVD under pure model or under contamination.

# rsvddpd 1.0.1
- Minor bug fixes and addition of user interruption during iterations with longer execution time
- Allowing robustness parameters beyond 1, with an issue of warning instead of stopping the execution.
- Added `cv.alpha` function for finding the optimal robustness parameter.
- Added different convergence rule based warnings.