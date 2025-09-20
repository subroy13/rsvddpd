## Test environments
* local macOS-arm64 R 4.5.1

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Resubmission
This is a resubmission to address a NOTE from the CRAN team.

In this version I have:

*   Removed the explicit `CXX_STD = CXX11` specification from `src/Makevars` and `src/Makevars.win` as requested. The package now uses the default C++ standard provided by R, which resolves the compatibility issue with `RcppArmadillo`.

*   Added a new function `rank.rSVDdpd` which is used to determine the rank for matrix factorization. 

*   Introduced a new argument `maxrank` to the `rSVDdpd` function for automatic matrix rank determination.

*   Updated the vignettes accordingly to reflect the changes due to modification of `rSVDdpd` function.

*   Updated the `DESCRIPTION` file:
    *   Incremented the package version to `1.0.1`.
    *   Updated the DOI identifier of arXiv paper to resolve another NOTE.
    *   Added `V8` to `Suggests` to resolve a NOTE about math rendering in vignettes.

*   Updated the `NEWS.md` file to reflect the changes in this version.

Thank you for your consideration.