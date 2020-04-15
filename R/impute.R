#'Imputation of missing cells in two-way data sets
#'
#'AMMI or GGE methods require the data to be complete, that is, no missing cells
#'are allowed. This function offers several methods to impute the missing cells,
#'so that the previously mentioned models can be adjusted later.
#'
#'@param Data a data frame
#'@param genotype name of the column that contains the genotypes
#'@param environment name of the column that contains the environments
#'@param response name of the column that contains the response
#'@param rep name of the column that contains the replications.If this argument
#'  is NULL, there is no replications in the data.
#'@param type imputation method. Either "EM-AMMI", "EM-SVD",
#'  "Gabriel","WGabriel","EM-PCA".
#'@inheritParams EM.AMMI
#'@inheritParams bcv::impute.svd
#'@inheritParams missMDA::imputePCA
#'@inheritParams Gabriel.Calinski
#'@inheritParams WGabriel
#'@return matrix with imputed data
#'@references Paderewski, J. (2013). An R function for imputation of missing
#'  cells in two-way data sets by EM-AMMI algorithm. Communications in Biometry
#'  and Crop Science 8 (2), 60–69.
#'@references FALTAN CITAS....GABRIEL EMPCA...
#'
#'
#'@export
#'
#' @examples
#' library(geneticae)
#' # Data without replication
#' data(yan.winterwheat)
#' dat <- yan.winterwheat
#' # generates missing data
#' dat[1,3]<-NA
#' dat[3,3]<-NA
#' dat[2,3]<-NA
#' \dontrun{
#' imputation(dat, genotype="gen",environment="env", response="yield",
#' type="EM-SVD") imputation(dat, genotype="gen",environment="env",
#' response="yield", type="EM-AMMI") imputation(dat,
#' genotype="gen",environment="env", response="yield", type="Gabriel")
#' imputation(dat, genotype="gen",environment="env", response="yield",
#' type="WGabriel") imputation(dat, genotype="gen",environment="env",
#' response="yield", type="EM-PCA")
#'}
#'# Data with replication
#'  data(plrv)
#'  dat2<-plrv
#'  dat2[1,3]<-NA
#' dat2[3,3]<-NA
#' dat2[2,3]<-NA
#' \dontrun{
#' imputation(dat2, genotype="Genotype",environment="Locality",
#' response="Yield", rep= "Rep", type="EM-SVD") imputation(dat2,
#' genotype="Genotype",environment="Locality", response="Yield", rep= "Rep",
#' type="EM-AMMI") imputation(dat2, genotype="Genotype",environment="Locality",
#' response="Yield", rep= "Rep", type="Gabriel") imputation(dat2,
#' genotype="Genotype",environment="Locality", response="Yield", rep= "Rep",
#' type="WGabriel") imputation(dat2, genotype="Genotype",environment="Locality",
#' response="Yield", rep= "Rep", type="EM-PCA")
#'}
#'
#'@importFrom stats var
#'@importFrom bcv impute.svd
#'@importFrom missMDA imputePCA
#'@importFrom dplyr group_by summarise rename
#'
imputation <- function(Data, genotype="gen",environment="env", response="yield", rep=NULL,type="EM-AMMI",
                       PC.nb=1, initial.values=NA, precision=0.01, max.iter=1000, change.factor=1, simplified.model=FALSE,
                       k = min(nrow(Data), ncol(Data)), tol = max(nrow(Data), ncol(Data)) * 1e-10, maxiter = 100,
                       Winf=0.8,Wsup=1,
                       ncp = 2, scale = TRUE, method = c("Regularized","EM"),
                       row.w = NULL, coeff.ridge = 1, threshold = 1e-06, seed = NULL, nb.init = 1) {

  if (missing(Data)) stop("Need to provide Data data frame")
  if (!any(is.na(Data))) stop("There are not missing data in input data frame")
    stopifnot(
    class(Data) %in%  c("matrix", "data.frame"),
    type %in% c("EM-AMMI", "EM-SVD","Gabriel","WGabriel","EM-PCA"),
    class(rep)%in% c("character", "NULL")
  )

  # if(!is.null(rep)){
  #   Data <-
  #     Data %>%
  #     group_by({{genotype}}, {{environment}}) %>%
  #     summarise(mean_resp=mean({{response}}))%>%
  #     plyr::rename(gen={{genotype}}, env= {{environment}}, mean_resp) %>%
  #     spread(env, mean_resp) %>%
  #     as.data.frame()
  #
  # } else{
  #   Data <-
  #     Data %>%
  #     spread({{environment}}, {{response}})%>%
  #     as.data.frame()
  # }
  #

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

  rownames(Data) <- pull(Data, genotype)
  Data <- dplyr::select(Data, -!!sym(genotype))
  # Data[,{{genotype}}] <- NULL


  if(type=="EM-AMMI"){
     matrix<-EM.AMMI(Data, PC.nb=PC.nb, initial.values=initial.values, precision=precision,
                     max.iter=max.iter, change.factor=change.factor, simplified.model=simplified.model)$X
   }
   if(type=="EM-SVD"){
     matrix<-impute.svd(Data, k = k, tol = tol, maxiter = maxiter)$x
   }
   if(type=="Gabriel"){
     matrix<-Gabriel.Calinski(Data)$GabrielImput
   }
   if(type=="WGabriel"){
     matrix<-WGabriel(Data,Winf,Wsup)$GabrielWImput
   }
   if(type=="EM-PCA"){
     matrix<-imputePCA(Data, ncp = ncp, scale = scale, method = method,
                       row.w = row.w, coeff.ridge = coeff.ridge, threshold = threshold, seed = seed, nb.init = nb.init,
                       maxiter = maxiter)$completeObs
   }

return(matrix)
}
