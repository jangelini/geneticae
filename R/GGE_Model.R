#'Site regression model
#'
#'@description The site regression model also named genotype plus
#'  genotype-by-environment (GGE) model is a powerful tool for effective
#'  analysis and interpretation of multienvironment data structure in breeding
#'  programs.This function is a modification of
#'  \code{\link[GGEBiplots]{GGEModel}} of GGEBiplots package, where the input
#'  dataset contains genotype by environment means with the genotypes in rows
#'  and environments in columns. By contrast, this function allows a less
#'  restrictive format of dataset. Repetitions and other variables that will not
#'  be used in the analysis may be present in the dataset.
#'
#'@param Data dataframe with genotypes, environments, repetitions (if any) and
#'  the phenotypic trait of interest. The order of the variables is indistinct,
#'  even additional variables that will not be used in the model may be present
#'  in the data.
#'@param genotype column name containing genotypes.
#'@param environment column name containing environments.
#'@param response column name containing phenotypic trait.
#'@param rep column name containing replications,if this argument is NULL, there
#'  is no replication available on the data. Defaults to NULL.
#'@param centering centering method. Either "tester" for tester centered (G+GE),
#'  "global" for global centered (E+G+GE), "double" for double centred (GE) or
#'  "none" for no centering. Default to "tester".
#'@param scaling scaling method. Either "sd" for standard deviation or "none"
#'  for no scaling. Default to "none".
#'@param SVP method for singular value partitioning. Either "row","column",
#'  "dual" or "symmetrical". Default to "column".
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
#'@references Sam Dumble (2017). GGEBiplots: GGE Biplots with 'ggplot2'. R
#'  package version 0.1.1. \url{https://CRAN.R-project.org/package=GGEBiplots}
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
#'  data(yan.winterwheat)
#'  GGE1 <- GGEmodel(yan.winterwheat, genotype="gen",environment="env", response="yield",
#'  centering = "tester")
#'
#'  # Data with replication
#'  data(plrv)
#'  GGE2 <- GGEmodel(plrv, genotype="Genotype",environment="Locality",
#'  response="Yield", rep="Rep", centering = "tester")
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
