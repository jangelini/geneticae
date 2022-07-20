#'Site Regression model
#'
#'@description The Site Regression model (also called genotype +
#'  genotype-by-environment (GGE) model) is a powerful tool for effective
#'  analysis and interpretation of data from multi-environment trials in
#'  breeding programs. There are different functions in R to fit the SREG model,
#'  such as the \code{\link[GGEBiplots]{GGEModel}} from the
#'  \href{https://CRAN.R-project.org/package=GGEBiplots}{GGEBiplots package}.
#'  However, this function has the following improvements: \itemize{ \item
#'  Includes recently published robust versions of the SREG model (Angelini et
#'  al., 2022). \item It can be used for data from trials with repetitions
#'  (there is no need to calculate means beforehand). \item Other variables not
#'  used in the analysis can be present in the dataset.}
#'
#'@param Data dataframe with genotypes, environments, repetitions (if any) and
#'  the phenotypic trait of interest. Additional variables that will not be used
#'  in the model may be present in the data.
#'@param genotype column name for genotypes.
#'@param environment column name for environments.
#'@param response column name for the phenotypic trait.
#'@param rep column name for replications. If this argument is NULL, there are
#'  no replications in the data. Defaults to NULL.
#'@param model method for fitting the SREG model: `"SREG"`,`"CovSREG"`,`"hSREG"`
#'  or `"ppSREG"` (see References). Defaults to `"SREG"`.
#'@param SVP method for singular value partitioning. Either `"row"`, `"column"`,
#'  or `"symmetrical"`. Defaults to `"symmetrical"`.
#'@return A list of class \code{GGE_Model} containing: \item{model}{SREG model
#'  version.} \item{coordgenotype}{plotting coordinates for each genotype in
#'  every component.} \item{coordenviroment}{plotting coordinates for each
#'  environment in every component.} \item{eigenvalues}{vector of eigenvalues
#'  for each component.} \item{vartotal}{overall variance.}
#'  \item{varexpl}{percentage of variance explained by each component.}
#'  \item{labelgen}{genotype names.} \item{labelenv}{environment names.}
#'  \item{axes}{axis labels.} \item{Data}{scaled and centered input data.}
#'  \item{SVP}{name of SVP method.}
#'
#'@return A biplot of class \code{ggplot}
#'
#'@details A linear model by robust regression using an M estimator proposed by
#'  Huber (1964, 1973) fitted by iterated re-weighted least squares, in
#'  combination with three robust SVD/PCA procedures, resulted in a total of
#'  three robust SREG alternatives. The robust SVD/PCA considered were:
#'  \itemize{ \item CovSREG: robust PCA that is obtained by replacing the
#'  classical estimates of location and covariance by their robust analogues
#'  using Minimum Regularized Covariance Determinant (MRCD) approach; \item
#'  hSREG: robust PCA method that tries to combine the advantages of both
#'  approaches, PCA based on a robust covariance matrix and based on projection
#'  pursuit; \item ppSREG: robust PCA that uses the projection pursuit and
#'  directly calculates the robust estimates of the eigenvalues and eigenvectors
#'  without going through robust covariance estimation. It is a very attractive
#'  method for bigdata situations, which are very common in METs (a few
#'  genotypes tested in a large number of environments), as the principal
#'  components can be calculated sequentially. }
#'@references  Julia Angelini, Gabriela Faviere, Eugenia Bortolotto, Gerardo
#'  Domingo Lucio Cervigni & Marta Beatriz Quaglino (2022) Handling outliers in
#'  multi-environment trial data analysis: in the direction of robust SREG
#'  model, Journal of Crop Improvement, DOI: 10.1080/15427528.2022.2051217
#'@export
#'@examples
#'  library(geneticae)
#'
#'  # Data without replication
#'  library(agridat)
#'  data(yan.winterwheat)
#'  GGE1 <- GGEmodel(yan.winterwheat, genotype="gen", environment="env", response="yield")
#'
#'  # Data with replication
#'  data(plrv)
#'  GGE2 <- GGEmodel(plrv, genotype = "Genotype", environment = "Locality",
#'                   response = "Yield", rep = "Rep")
#'
#'@importFrom MASS rlm
#'@importFrom pcaMethods robustSvd
#'@importFrom rrcov PcaHubert PcaCov PcaProj CovControlMrcd getEigenvalues
#'  getLoadings getScores
#'@importFrom stats lm residuals
#'@importFrom dplyr group_by summarise rename pull %>%
#'@importFrom rlang sym

GGEmodel<- function(Data, genotype = "gen", environment = "env", response = "yield",
                      rep=NULL, model = "SREG", SVP="symmetrical")
{

  if (missing(Data)) stop("Need to provide Data data frame")
  if(any(is.na(Data))){stop("Missing data in input data frame, run the imputation function first to complete the data set")}
  stopifnot(
    class(Data) == "data.frame",
    class(genotype) == "character",
    class(environment) == "character",
    class(response) == "character",
    class(rep)%in% c("character", "NULL"),
    model %in% c("SREG", "hSREG", "ppSREG", "CovSREG"),
    SVP %in% c("symmetrical", "row", "column")
  )

  if(!is.null(rep)){
    Data <-
      Data %>%
      group_by(!!sym(environment), !!sym(genotype)) %>%
      summarise(y=mean(!!sym(response)))
    response = "y"
  }

  gen <- as.factor(pull(Data, genotype))
  env <- as.factor(pull(Data, environment))
  y <- pull(Data, response)

  Ngen <- nlevels(gen)
  Nenv <- nlevels(env)
  labelgen <- unique(gen)
  labelenv <- unique(env)
  Ngen <- nlevels(gen)
  Nenv <- nlevels(env)
  Max <- min(Ngen - 1, Nenv)

  if (model == "SREG"){
    # classic regression
    lm.x<- lm(y ~ env, data=Data, contrasts=list(env="contr.sum"))
    residuals.x<- matrix(residuals(lm.x), nrow = Ngen , ncol = Nenv)
    svd.x<- svd(residuals.x)
    singlevalue<-svd.x$d
    eigenvalues = singlevalue^2
    loading<-(svd.x$v)
    score<-(svd.x$u%*%diag(svd.x$d))
  }
  #
  else if (model == "CovSREG"){
    # robust regression
    rlm.x<- rlm(y~env, data=Data, contrasts=list(env="contr.sum"),maxit = 50)
    residuals.x<- matrix(residuals(rlm.x), nrow = Ngen, ncol = Nenv)

    PcaCov <- PcaCov(residuals.x, cov.control=CovControlMrcd())

    loading<- getLoadings(PcaCov)
    score<- getScores(PcaCov)

    singlevalue<- sqrt(getEigenvalues(PcaCov))
    eigenvalues = singlevalue^2
    n <- NROW(score)
    singlevalue <- singlevalue * sqrt(n)
  }
  else if (model == "hSREG"){
    # robust regression
    rlm.x<- rlm(y~env, data=Data, contrasts=list(env="contr.sum"),maxit = 50)
    residuals.x<- matrix(residuals(rlm.x), nrow = Ngen, ncol = Nenv)

    modeloHubert<- PcaHubert(residuals.x, mad=FALSE)

    loading<- getLoadings(modeloHubert)
    score<- getScores(modeloHubert)

    singlevalue<- sqrt(getEigenvalues(modeloHubert))
    eigenvalues = singlevalue^2
    n <- NROW(score)
    singlevalue <- singlevalue * sqrt(n)
  }
  else if (model == "ppSREG"){
    # robust regression
    rlm.x<- rlm(y~env, data=Data, contrasts=list(env="contr.sum"),maxit = 50)
    residuals.x<- matrix(residuals(rlm.x), nrow = Ngen, ncol = Nenv)

    modeloProj<- PcaProj(residuals.x)
    loading<- getLoadings(modeloProj)
    score<- getScores(modeloProj)

    singlevalue<- sqrt(getEigenvalues(modeloProj))
    eigenvalues = singlevalue^2
    n <- NROW(score)
    singlevalue <- singlevalue * sqrt(n)
  }

  coordgenotype <- score
  coordenviroment <- loading
  singlevalue_scale <- singlevalue


  if(SVP=="row"){
    scale=0
  }
  if(SVP=="column"){
    scale=1
  }
  if(SVP=="symmetrical"){
    scale=0.5
  }

  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) singlevalue_scale <- singlevalue_scale^scale else singlevalue_scale <- 1

  coordgenotype <- t(t(score) / singlevalue_scale)
  coordenviroment <- t(t(loading) * singlevalue_scale)

  vartotal = round(as.numeric(sum(eigenvalues)),2)
  varexpl = round(as.numeric((eigenvalues/vartotal) *100), 2)

  centering="tester"
  scaling = "none"
  labelaxes <- paste("Component ",1:ncol(diag(singlevalue)), sep = "")

  GGEmodel=list(model=model,
                  coordgenotype=coordgenotype,
                  coordenviroment=coordenviroment,
                  eigenvalues=eigenvalues,
                  vartotal=vartotal,
                  varexpl=varexpl,
                  labelgen=labelgen,
                  labelenv=labelenv,
                  labelaxes=labelaxes,
                  Data=residuals.x,
                  SVP=SVP)
  class(GGEmodel)<-"GGEModel"
  return(GGEmodel)
}

