#' Produces genotype plus genotype-by-environment model from a 2-way table of
#' means
#'
#'
#' @param Data a data frame or matrix that contains the genotypes in the rows and the
#' environments in the columns when there are no replications of the experiment. In case
#' of replications, the data frame must have the columns in the following order: genotypes,
#' environments, replications and the response variable.
#' @param rep logical. If TRUE the genotype by environment means is calculated.
#' @param type imputation method. Either "EM-AMMI", "EM-SVD", "Gabriel","WGabriel"
#' @export
#'
#' @importFrom stats var
#' @importFrom GGEBiplots stattable
#' @importFrom bcv impute.svd
#' @importFrom missMDA imputePCA
#'
imputation <- function(Data,rep=FALSE,type="EM-AMMI") {
  if(rep==TRUE){
    Data<-stattable(Data[,1],Data[,2],Data[,4],FUN=mean)
    if(type=="EM-AMMI"){
      matrix<-EM.AMMI(Data)
    }
    if(type=="EM-SVD"){
      matrix<-impute.svd(Data)
    }
    if(type=="Gabriel"){
      matrix<-Gabriel.Calinski(Data)
    }
    if(type=="WGabriel"){
      matrix<-WGabriel(Data,0.8,1)
    }
    if(type=="EM-PCA"){
      matrix<-imputePCA(Data)
    }
  }

  if(rep==FALSE){
    if(type=="EM-AMMI"){
        matrix<-EM.AMMI(Data)
    }
    if(type=="EM-SVD"){
      matrix<-impute.svd(Data)
    }
    if(type=="Gabriel"){
      matrix<-Gabriel.Calinski(Data)
    }
    if(type=="WGabriel"){
      matrix<-WGabriel(Data,0.8,1)
    }
    if(type=="EM-PCA"){
      matrix<-imputePCA(Data)
    }

  }

return(matrix)

}
