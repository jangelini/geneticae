#'Imputation of missing cells in two-way data sets
#'
#'Missing values are not allowed by the AMMI or GGE methods. This function
#'provides several methods to impute missing observations in data from
#'multi-environment trials and to subsequently adjust the mentioned methods.
#'
#'@param Data dataframe containing genotypes, environments, repetitions (if any)
#'  and the phenotypic trait of interest. Columns can be in any order in the
#'  dataset and other variables that will not be used in the analysis can be
#'  present.
#'@param genotype column name containing genotypes.
#'@param environment column name containing environments.
#'@param response column name containing the phenotypic trait.
#'@param rep column name containing replications. If this argument is NULL,
#'  there are no replications available in the data. Defaults to NULL.
#'@param type imputation method. Either "EM-AMMI", "EM-SVD",
#'  "Gabriel","WGabriel","EM-PCA". Defaults to "EM-AMMI".
#'@param nPC number of components used to predict the missing values.
#'  Default to 2.
#'@param initial.values initial values of the missing cells. It can be a single
#'  value or a vector of length equal to the number of missing cells (starting
#'  from the missing values in the first column). If omitted, the initial values
#'  will be obtained by the main effects from the corresponding model, that is,
#'  by the grand mean of the observed data increased (or decreased) by row and
#'  column main effects.
#'@param precision threshold for assessing convergence.
#'@param maxiter maximum number of iteration for the algorithm.
#'@param change.factor  When `change.factor` is equal to 1, the previous
#'  approximation is changed with the new values of missing cells (standard
#'  EM-AMMI algorithm). However, when `change.factor` less than 1, then the new
#'  approximations are computed and the values of missing cells are changed in
#'  the direction of this new approximation but the change is smaller. It could
#'  be useful if the changes are cyclic and thus convergence could not be
#'  reached. Usually, this argument should not affect the final outcome (that
#'  is, the imputed values) as compared to the default value of `change.factor`
#'  = 1.
#'@param simplified.model the AMMI model contains the general mean, effects of
#'  rows, columns and interaction terms. So the EM-AMMI algorithm in step 2
#'  calculates the current effects of rows and columns; these effects change
#'  from iteration to iteration because the empty (at the outset) cells in each
#'  iteration are filled with different values. In step 3 EM-AMMI uses those
#'  effects to re-estimate cells marked as missed (as default,
#'  simplified.model=FALSE). It is, however, possible that this procedure will
#'  not converge. Thus the user is offered a simplified EM-AMMI procedure that
#'  calculates the general mean and effects of rows and columns only in the
#'  first iteration and in next iterations uses these values
#'  (simplified.model=TRUE). In this simplified procedure the initial values
#'  affect the outcome (whilst EM-AMMI results usually do not depend on initial
#'  values). For the simplified procedure the number of iterations to
#'  convergence is usually smaller and, furthermore, convergence will be reached
#'  even in some cases where the regular procedure fails. If the regular
#'  procedure does not converge for the standard initial values, the simplified
#'  model can be used to determine a better set of initial values.
#'@param k rank of the SVD approximation.
#'@param scale boolean. By default TRUE leading to a same weight for each
#'  variable
#'@param method "Regularized" by default or "EM"
#'@param row.w row weights (by default, a vector of 1 for uniform row weights)
#'@param coeff.ridge 1 by default to perform the regularized imputePCA
#'  algorithm; useful only if method="Regularized". Other regularization terms
#'  can be implemented by setting the value to less than 1 in order to
#'  regularized less (to get closer to the results of the EM method
#'@param seed integer, by default seed = NULL implies that missing values are
#'  initially imputed by the mean of each variable. Other values leads to a
#'  random initialization
#'@param nb.init integer corresponding to the number of random initializations;
#'  the first initialization is the initialization with the mean imputation
#'@param Winf peso inferior
#'@param Wsup peso superior
#'
#'@return imputed data matrix
#'@references Paderewski, J. (2013). An R function for imputation of missing
#'  cells in two-way data sets by EM-AMMI algorithm. Communications in Biometry
#'  and Crop Science 8, 60–69.
#'@references Julie Josse, Francois Husson (2016). missMDA: A Package for
#'  Handling Missing Values in Multivariate Data Analysis. Journal of
#'  Statistical Software 70, 1-31.
#'@references Arciniegas-Alarcón S., García-Peña M., Dias C.T.S., Krzanowski
#'  W.J. (2010). \emph{An alternative methodology for imputing missing data in
#'  trials with genotype-by-environment interaction}. Biometrical Letters 47,
#'  1–14.
#'@references  Perry P.O. (2015). bcv: Cross-Validation for the SVD
#'  (Bi-Cross-Validation). R package version 1.0.1.
#'@references Arciniegas-Alarcón S., García-Peña M., Krzanowski W.J., Dias
#'  C.T.S. (2014). An alternative methodology for imputing missing data in
#'  trials with genotype-byenvironment interaction: some new aspects.
#'  Biometrical Letters 51, 75-88.
#'
#'@export
#'
#' @examples
#' library(geneticae)
#' # Data without replications
#' data(yan.winterwheat)
#'
#' # generating missing values
#' yan.winterwheat[1,3]<-NA
#' yan.winterwheat[3,3]<-NA
#' yan.winterwheat[2,3]<-NA
#'
#' \dontrun{
#' imputation(yan.winterwheat, genotype = "gen", environment = "env",
#'            response = "yield", type = "EM-AMMI")
#' }
#' # Data with replications
#' data(plrv)
#' plrv[1,3] <- NA
#' plrv[3,3] <- NA
#' plrv[2,3] <- NA
#' \dontrun{
#' imputation(plrv, genotype = "Genotype", environment = "Locality",
#'            response = "Yield", rep = "Rep", type = "EM-SVD")
#' }
#'
#'@importFrom stats var
#'@importFrom bcv impute.svd
#'@importFrom missMDA imputePCA
#'@importFrom dplyr group_by summarise rename
#'
imputation <- function(Data, genotype="gen",environment="env", response="yield", rep=NULL,type="EM-AMMI",
                       nPC=2, initial.values=NA, precision=0.01, maxiter=1000, change.factor=1, simplified.model=FALSE,
                       k = min(nrow(Data), ncol(Data)), scale = TRUE, method = "EM",
                       row.w = NULL, coeff.ridge = 1, seed = NULL, nb.init = 1, Winf=0.8,Wsup=1) {

  if (missing(Data)) stop("Need to provide Data data frame")
  if (!any(is.na(Data))) stop("There are not missing data in input data frame")
    stopifnot(
    class(Data) %in% c("data.frame"),
    class(rep)%in% c("character", "NULL"),
    class(genotype) == "character",
    class(environment) == "character",
    class(response) == "character",
    type %in% c("EM-AMMI", "EM-SVD","Gabriel","WGabriel","EM-PCA"),
    class(nPC) == "numeric",
    # class(initial.values)  %in% c(NA, "vector", "numeric"),
    class(precision) == "numeric",
    class(maxiter) == "numeric",
    class(change.factor) == "numeric",
    class(simplified.model) == "logical",
    # class(k) == "numeric",
    class(scale) == "logical",
    method %in% c("Regularized", "EM"),
    class(row.w)  %in% c("NULL", "vector", "numeric"),
    class(coeff.ridge) == "numeric",
    class(seed) %in% c("NULL", "numeric"),
    class(nb.init) == "numeric",
    class(Winf) == "numeric",
    class(Wsup) == "numeric"
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

  rownames(Data) <- pull(Data, genotype)
  Data <- dplyr::select(Data, -!!sym(genotype))

  if(type=="EM-AMMI"){
     matrix<-EM.AMMI(Data, PC.nb=nPC, initial.values=initial.values, precision=precision,
                     max.iter=maxiter, change.factor=change.factor, simplified.model=simplified.model)$X
   }
   if(type=="EM-SVD"){
     matrix<-impute.svd(Data, k = k, tol = precision, maxiter = maxiter)$x
   }
   if(type=="Gabriel"){
     matrix<-Gabriel.Calinski(Data)$GabrielImput
   }
   if(type=="WGabriel"){
     matrix<-WGabriel(Data,Winf,Wsup)$GabrielWImput
   }
   if(type=="EM-PCA"){
     matrix<-imputePCA(Data, ncp = nPC, scale = scale, method = method,
                       row.w = row.w, coeff.ridge = coeff.ridge, threshold = precision, seed = seed, nb.init = nb.init,
                       maxiter = maxiter)$completeObs
   }

return(matrix)
}
