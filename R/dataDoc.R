#' Clones from the PLRV population
#'
#' @description  Resistance study to PLRV (Patato Leaf Roll Virus), the disease
#'   is leaf curl due to the effect of the virus. 28 genotypes were experimented
#'   in 6 locations in Peru, being 504 the total number of records. This dataset
#'   was obtained from agricolae Package.
#'
#' @name plrv
#' @docType data
#' @usage data(plrv)
#'
#' @format A data frame with 504 observations on the following 6 variables
#'
#'   Genotype a factor with levels 102.18 104.22 121.31 141.28 157.26 163.9
#'   221.19 233.11 235.6 241.2 255.7 314.12 317.6 319.20 320.16 342.15 346.2
#'   351.26 364.21 402.7 405.2 406.12 427.7 450.3 506.2 Canchan Desiree Unica.
#'
#'   Locality a factor with levels Ayacucho (Ayac), La Molina 02 (LM-02), San
#'   Ramon 02 (SR-02), Huancayo (Hyo-02), La Molina 03 (LM-03), San Ramon 03
#'   (SR-03).
#'
#'   Rep a numeric vector
#'
#'   WeightPlant, WeightPlot and Yield are numeric vectors
#'
#' @keywords datasets
#'
#' @references Felipe de Mendiburu (2020). agricolae: Statistical Procedures for
#'   Agricultural Research. R package version 1.3-2.
#'   \url{https://CRAN.R-project.org/package=agricolae}
#'
#' @examples
#' library(geneticae)
#' data(plrv)
#' str(plrv)
#'
NULL

#' Winter wheat varieties from Ontario
#'
#' @description Yield of 18 winter wheat varieties grown at 9 environments in
#'   Ontario in 1993. This dataset was obtained from agridat Package.
#'
#' @name yan.winterwheat
#' @docType data
#' @usage data(yan.winterwheat)
#'
#' @format A data frame with 162 observations on the following 3 variables.
#'
#'   gen genotype
#'
#'   env environment
#'
#'   yield yield in metric tons per hectare
#'
#' @keywords datasets
#' @references Kevin Wright (2018). agridat: Agricultural Datasets. R package
#'   version 1.16.\url{https://CRAN.R-project.org/package=agridat}
#'
#' @examples
#' library(geneticae)
#' data(yan.winterwheat)
#' str(yan.winterwheat)
NULL
