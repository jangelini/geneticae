## New submission to CRAN
This is a new release to CRAN. In this version, we have:

* Included new imputation methods specifically developed for multi-environment trials based on Angelini et al. (2024).
* Renamed several functions to use more descriptive and appropriate names for the user (e.g., rAMMI to rAMMIModel).
* Fixed an issue in `rAMMIModel.R` where datasets with replicates were not being handled correctly when calculating genotype-by-environment means.
* Added a comprehensive tutorial vignette ("Tutorial") covering all new features and SREG robust model.

## Test environments
* local Linux Mint 22.2 Cinnamon, R 4.3.3
* win-builder (devel, release and oldrelease)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs, other than the NEW SUBMISSION note.

## Downstream dependencies
There are currently no downstream dependencies for this package.
