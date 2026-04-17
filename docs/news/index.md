# Changelog

## geneticae 1.0.1

- **Documentation Update**:
  - Major improvements to
    [`rAMMIModel()`](https://github.com/jangelini/geneticae/reference/rAMMIModel.md)
    and
    [`rAMMIPlot()`](https://github.com/jangelini/geneticae/reference/rAMMIPlot.md)
    documentation to clarify robust estimation methods and return
    objects.
  - Updated package `Description` with latest references for robust
    SREG (2022) and MET data imputation (2024).

## geneticae 1.0.0

CRAN release: 2026-04-16

**New features:**

- Included new imputation methods specifically developed for
  multi-environment trials based on Angelini et al. (2024).

**Renamed functions:**

- Several functions were renamed to be more descriptive and consistent:
  rAMMI is now rAMMIModel, rSREG is now rSREGModel, and GGE is now
  GGEModel.

**Bug fixes:**

- Fixed an issue in rAMMIModel.R where datasets with replicates were not
  being handled correctly when calculating genotype-by-environment
  means.

**Documentation:**

- Added a comprehensive tutorial vignette (“Tutorial”) covering all new
  features and SREG robust model.

## geneticae 0.4.0

CRAN release: 2022-07-20

Added SREG robust model in GGEmodel() function.

Dependecy on archived bcv package removed.

## geneticae 0.3.0

CRAN release: 2022-02-09

Dependecy on archived GGEBiplotGUI package removed.

## geneticae 0.2.0

CRAN release: 2021-12-20

Changed default option to “symmetrical” in the SVP argument of the
GGEmodel() function.

Changed agridat package from Imports to Suggests in DESCRIPTION file.

## geneticae 0.1.0

CRAN release: 2021-09-16

First submission to CRAN.

## geneticae 0.0.9000

This is the first development version of the package.
