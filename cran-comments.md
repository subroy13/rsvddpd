## Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.
There was 2 NOTEs. According to [stackoverflow](https://stackoverflow.com/questions/64402688/information-on-o-files-for-x64-is-not-available-note-on-r-package-checks-using), this note is a false positive and can be ignored.

## Resubmission
This is a resubmission. In this version I have:

* Updated the vignette to remove a demo from pcaMethods package.
* Removed dependencies of RobRSVD and pcaMethods packages from DESCRIPTION file.
* Provided a reference to the algorithm used to compute the decomposition. 
* Set up github action R CMD checks.
* Changed the url text appropriately.