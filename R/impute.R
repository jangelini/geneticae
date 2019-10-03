#' Produces genotype plus genotype-by-environment model from a 2-way table of
#' means
#'
#' @param Data a data frame or matrix that contains the genotypes in the rows and the
#'   environments in the columns when there are no replications of the experiment. In case
#'   of replications, the data frame must have the columns in the following order: genotypes,
#'   environments, replications and the response variable.
#' @param rep logical. If TRUE the genotype by environment means is calculated.
#' @param type imputation method. Either "EM-AMMI", "EM-SVD", "Gabriel","WGabriel"
#' @return matrix with imputed data
#' @export
#'
#' @examples
#' library(geneticae)
#' library(reshape2)
#' data(yan.winterwheat)
#' dat1 <- yan.winterwheat
#' dat <-t(round(acast(dat1, env~gen, value.var='yield'),2))
#' # generates missing data
#' dat[1,2]<-NA
#' dat[3,4]<-NA
#' dat[2,7]<-NA
#' imputation(dat, rep=FALSE, type="EM-SVD")
#' imputation(dat, rep=FALSE, type="EM-AMMI")
#' imputation(dat, rep=FALSE, type="Gabriel")
#' imputation(dat, rep=FALSE, type="WGabriel")
#' imputation(dat, rep=FALSE, type="EM-PCA")
#'
#' @importFrom stats var
#' @importFrom GGEBiplots stattable
#' @importFrom bcv impute.svd
#' @importFrom missMDA imputePCA
#'
imputation <- function(Data,rep=FALSE,type="EM-AMMI") {

  if (missing(Data)) stop("Need to provide Data data frame")
  if (!any(is.na(Data))) stop("There are not missing data in input data frame")
    stopifnot(
    class(Data) %in%  c("matrix", "data.frame"),
    type %in% c("EM-AMMI", "EM-SVD","Gabriel","WGabriel","EM-PCA"),
    class(rep) == "logical"
  )


  if(rep==TRUE){
    Data<-stattable(Data[,1],Data[,2],Data[,4],FUN=mean)
    if(type=="EM-AMMI"){
      matrix<-EM.AMMI(Data)$X
    }
    if(type=="EM-SVD"){
      matrix<-impute.svd(Data)$x
    }
    if(type=="Gabriel"){
      matrix<-Gabriel.Calinski(Data)$GabrielImput
    }
    if(type=="WGabriel"){
      matrix<-WGabriel(Data,0.8,1)$GabrielWImput
    }
    if(type=="EM-PCA"){
      matrix<-imputePCA(Data)$completeObs
    }
  }

  if(rep==FALSE){
    if(type=="EM-AMMI"){
        matrix<-EM.AMMI(Data)$X
    }
    if(type=="EM-SVD"){
      matrix<-impute.svd(Data)$x
    }
    if(type=="Gabriel"){
      matrix<-Gabriel.Calinski(Data)$GabrielImput
    }
    if(type=="WGabriel"){
      matrix<-WGabriel(Data,0.8,1)$GabrielWImput
    }
    if(type=="EM-PCA"){
      matrix<-imputePCA(Data)$completeObs
    }

  }

return(matrix)
}
