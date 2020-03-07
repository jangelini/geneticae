#'Site regression model
#'
#'@description The site regression (SREG) model also named genotype plus
#'  genotype-by-environment (GGE) model is powerful a tools for effective
#'  analysis and interpretation of multienvironment data structure in breeding
#'  programs.This function is a modification of \emph{GGEModel} of GGEBiplots
#'  package, where the data frame or matrix contains genotype by environment
#'  means with the genotypes in rows and environments in columns. By contrast,
#'  repetitions are allowed in this function.
#'
#'@param Data a data frame
#'@param genotype name of the column that contains the genotypes
#'@param environment name of the column that contains the environments
#'@param response name of the column that contains the response
#'@param rep name of the column that contains the replications.If this argument
#'  is NULL, there is no replications in the data.
#'@param centering centering method. Either "tester" for tester centered (G+GE),
#'  "global" for global centered (E+G+GE), "double" for double centred (GE) or
#'  "none" for no centering. If a centering method is not used, the
#'  \code{\link[geneticae]{GGEPlot}} function can not be used.
#'@param scaling scaling method. Either "sd" for standard deviation or "none"
#'  for no scaling.
#'@param SVP method for singular value partitioning. Either "row","column",
#'  "dual" or "symmetrical".
#'@return A list of class \code{GGE_Model} containing:
#'  \item{coordgenotype}{plotting coordinates for genotypes from all components}
#'  \item{coordenviroment}{plotting coordinates for environments from all
#'  components} \item{eigenvalues}{vector of eigenvalues from each component}
#'  \item{vartotal}{overall variance} \item{varexpl}{percentage of variance
#'  explained by each component} \item{labelgen}{genotype names}
#'  \item{labelenv}{environment names} \item{axes}{axis labels}
#'  \item{Data}{scaled and centered input data} \item{centering}{name of
#'  centering method} \item{scaling}{name of scaling method} \item{SVP}{name of
#'  SVP method}
#'@references Yan W, Kang M (2003). \emph{GGE Biplot Analysis: A Graphical Tool
#'  for Breeders, Geneticists, and Agronomists}. CRC Press.
#'@references Yan W, Kang M (2002). \emph{Singular-Value Partitioning in Biplot
#'  Analysis of Multienvironment Trial Data}. Agronomy Journal, 94, 990-996.
#'  \url{http://dx.doi.org/10.2134/agronj2002.0990}
#'@export
#' @examples
#'
#'  library(geneticae)
#'  # Data without replication
#'  GGE1 <- GGEmodel(yan.winterwheat, genotype="gen",environment="env", response="yield",
#'  centering = "tester")
#'  GGE1
#'
#'  # Data with replication
#'  data(plrv)
#'  GGE2 <- GGEmodel(plrv, genotype="Genotype",environment="Locality",
#'  response="Yield", rep="Rep", centering = "tester")
#'  GGE2
#'
#'@importFrom stats var
#'@importFrom GGEBiplots GGEModel
#'@importFrom tidyr spread
#'@importFrom dplyr group_by summarise rename
#'
GGEmodel <- function(Data, genotype="gen", environment="env", response="yield",
                     rep=NULL, centering="tester",scaling="none",SVP="column"){

  if (missing(Data)) stop("Need to provide Data data frame or matrix")
  if(any(is.na(Data))){stop("Missing data in input data frame, run the imputation function first to complete the data set")}
  stopifnot(
     class(Data) %in% c("data.frame"),
     class(rep)%in% c("character", "NULL"),
     class(genotype) == "character",
     class(environment) == "character",
     class(response) == "character",
     centering %in% c("tester", "global","double","none"),
     scaling %in% c("sd", "none"),
     SVP %in% c("row", "column","dual","symmetrical")
  )

  if(!is.null(rep)){
  Data <-
    Data %>%
    group_by(!!sym(genotype), !!sym(environment)) %>%
    summarise(mean_resp=mean(!!sym(response)))%>%
    spread(!!sym(environment), mean_resp) %>%
    as.data.frame()

  } else{
    Data <-
      Data %>%
      spread(!!sym(environment), !!sym(response)) %>%
      as.data.frame()

  }


  rownames(Data) <- Data[,{{genotype}}]
  Data[,{{genotype}}] <- NULL

  # GGEModel is a function of GGEBiplots package
  model<-GGEModel(Data,centering=centering,scaling=scaling,SVP=SVP)


  class(model)<-"GGEModel"
  return(model)
}
