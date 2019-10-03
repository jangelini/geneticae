#' Data clones from the PLRV population
#'
#' Six environments: Ayacucho, La Molina 02, San Ramon 02, Huancayo,
#' La Molina 03, San Ramon 03. This data.frame was obtained from agricolae Package.
#' @name plrv
#' @docType data
#' @usage data(plrv)
#'
#' @format A data frame with 504 observations on the following 6 variables.
#' Genotype a factor with levels 102.18 104.22 121.31 141.28 157.26 163.9 221.19 233.11
#' 235.6 241.2 255.7 314.12 317.6 319.20 320.16 342.15 346.2 351.26 364.21 402.7
#' 405.2 406.12 427.7 450.3 506.2 Canchan Desiree Unica
#' Locality a factor with levels Ayac Hyo-02 LM-02 LM-03 SR-02 SR-03
#' Rep a numeric vector
#' WeightPlant a numeric vector
#' WeightPlot a numeric vector
#' Yield a numeric vector
#'
#' @keywords datasets
#'
#' @references
#' AAAAA
#'
#' @examples
#' library(geneticae)
#' data(plrv)
#' str(plrv)
#'
NULL

#' Multi-environment trial of winter wheat in Ontario
#'
#' Yield of 18 varieties of winter wheat grown at 9 environments in Ontario in 1993.
#' This data.frame was obtained from agridat Package.
#'
#' @name yan.winterwheat
#' @docType data
#' @usage data(yan.winterwheat)
#'
#' @format A data frame with 162 observations on the following 3 variables.
#' gen genotype
#' env environment
#' yield yield in metric tons per hectare
#'
#' @keywords datasets
#' @references
#' Weikai Yan and Manjit S. Kang and Baoluo Ma and Sheila Woods, 2007, GGE Biplot vs. AMMI
#' Analysis of Genotype-by-Environment Data, Crop Science, 2007, 47, 641–653. \url{http://doi.org/10.2135/cropsci2006.06.0374}
#' @examples
#' library(geneticae)
#' data(yan.winterwheat)
#' str(yan.winterwheat)
NULL
