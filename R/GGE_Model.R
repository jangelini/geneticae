#'Site regression model
#'
#'@description The site regression (SREG) model also named genotype plus
#'genotype-by-environment (GGE) model is powerful a tools for effective
#'analysis and interpretation of multienvironment data structure in breeding
#'programs.This function is a modification of \emph{GGEModel} of GGEBiplots
#'package, where the data frame or matrix contains genotype by environment
#'means with the genotypes in rows and environments in columns. By contrast,
#'repetitions are allowed in this function.
#'
#' @param Data a data frame or matrix that contains the genotypes in the rows
#'   and the environments in the columns when there no replications in the
#'   experiment. In case of replications, the data frame must have the columns
#'   in the following order: genotypes, environments, replications and the
#'   response variable.
#' @param rep logical. If TRUE the genotype by environment means is calculated.
#' @param centering centering method. Either "tester" for tester centered
#'   (G+GE), "global" for global centered (E+G+GE), "double" for double centred
#'   (GE) or "none" for no centering. If a centering method is not used, the
#'   \code{\link[geneticae]{GGEPlot}} function can not be used.
#' @param scaling scaling method. Either "sd" for standard deviation or "none"
#'   for no scaling.
#' @param SVP method for singular value partitioning. Either "row","column",
#'   "dual" or "symmetrical".
#' @return A list of class \code{GGE_Model} containing:
#'   \item{coordgenotype}{plotting coordinates for genotypes from all
#'   components} \item{coordenviroment}{plotting coordinates for environments
#'   from all components} \item{eigenvalues}{vector of eigenvalues from each
#'   component} \item{vartotal}{overall variance} \item{varexpl}{percentage of
#'   variance explained by each component} \item{labelgen}{genotype names}
#'   \item{labelenv}{environment names} \item{axes}{axis labels}
#'   \item{Data}{scaled and centered input data} \item{centering}{name of
#'   centering method} \item{scaling}{name of scaling method} \item{SVP}{name
#'   of SVP method}
#' @references Yan W, Kang M (2003). \emph{GGE Biplot Analysis: A Graphical Tool
#'   for Breeders, Geneticists, and Agronomists}. CRC Press.
#' @references Yan W, Kang M (2002). \emph{Singular-Value Partitioning in Biplot
#'   Analysis of Multienvironment Trial Data}. Agronomy Journal, 94, 990-996.
#'   \url{http://dx.doi.org/10.2134/agronj2002.0990}
#' @export
#' @examples
#'   library(geneticae)
#'   data(Ontario)
#'   GGE1<-GGEmodel(Ontario, centering="tester", rep=FALSE)
#'   GGEPlot(GGE1)
#' @importFrom stats var
#' @importFrom GGEBiplots stattable GGEModel
#'
GGEmodel <- function(Data,rep=FALSE,centering="tester",scaling="none",SVP="column"){

  if (missing(Data)) stop("Need to provide Data data frame or matrix")

  if(any(is.na(Data))){stop("Missing data in input data frame, run the imputation function first to complete the data set")}

  # stopifnot(
  #   class(Data) %in% c("matrix", "data.frame"),
  #   class(rep) == "logical",
  #   class(centering) %in% c("tester", "global","double","none"),
  #   class(scaling) %in% c("sd", "none"),
  #   class(SVP) %in% c("row", "column","dual","symmetrical")
  # )


  if(rep==TRUE){
    Data<-stattable(Data[,1],Data[,2],Data[,4],FUN=mean)
    model<-GGEModel(Data,centering=centering,scaling=scaling,SVP=SVP)
  }

  if(rep==FALSE){
    model<-GGEModel(Data,centering=centering,scaling=scaling,SVP=SVP)
  }

  class(model)<-"GGEModel"
  return(model)
}
