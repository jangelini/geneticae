#' Imputation of missing cells in two-way data sets
#'
#' Missing values are not allowed by the AMMI, GGE or SREG methods. This function
#' provides several methods to impute missing observations in data from
#' multi-environment trials and to subsequently adjust the mentioned methods.
#'
#' @param Data dataframe containing genotypes, environments, repetitions (if any)
#'   and the phenotypic trait of interest. Other variables that will not be used
#'   in the analysis can be present.
#' @param genotype column name containing genotypes.
#' @param environment column name containing environments.
#' @param response column name containing the phenotypic trait.
#' @param rep column name containing replications. If this argument is NULL,
#'   there are no replications available in the data. Defaults to NULL.
#' @param type imputation method. Either "EM-AMMI", "EM-GGE", "EM-SREG", "EM-bSREG",
#'   "Gabriel", "Eigenvector", "WGabriel", "EM-PCA". Defaults to "EM-AMMI".
#' @param nPC number of components used to predict the missing values.
#'   Default to 2.
#' @param initial.values initial values of the missing cells. It can be a single
#'   value or a vector of length equal to the number of missing cells.
#' @param precision threshold for assessing convergence.
#' @param maxiter maximum number of iteration for the algorithm.
#' @param change.factor When `change.factor` is equal to 1, the previous
#'   approximation is changed with the new values (standard EM). Smaller values
#'   can help convergence if changes are cyclic.
#' @param simplified.model logical. If TRUE, calculates effects only in the first 
#'   iteration to speed up convergence or help in cases where the regular 
#'   procedure fails.
#' @param scale boolean. By default TRUE for "EM-PCA".
#' @param method "Regularized" or "EM" for "EM-PCA".
#' @param row.w row weights for "EM-PCA".
#' @param coeff.ridge ridge coefficient for "EM-PCA".
#' @param seed integer for random initialization in "EM-PCA".
#' @param nb.init number of random initializations for "EM-PCA".
#' @param Winf lower weight for WGabriel.
#' @param Wsup upper weight for WGabriel.
#'
#'@return A matrix of the imputed data.
#'
#'@details
#'
#'Often, multi-environment experiments are unbalanced because several genotypes
#'are not tested in some environments. Several methodologies have been proposed
#'in order to solve this lack of balance caused by missing values, some of which
#'are included in this function:
#'
#'\itemize{
#'\item EM-AMMI: an iterative scheme built round the above procedure is used to
#'obtain AMMI imputations from the EM algorithm. The additive parameters are
#'initially set by computing the grand mean, genotype means and environment
#'means obtained from the observed data. The residuals for the observed cells
#'are initialized as the cell mean minus the genotype mean minus the environment
#'mean plus the grand mean, and interactions for the missing positions are
#'initially set to zero. The initial multiplicative parameters are obtained from
#'the SVD of this matrix of residuals, and the missing values are filled by the
#'appropriate AMMI estimates. In subsequent iterations, the usual AMMI procedure
#'is applied to the completed matrix and the missing values are updated by the
#'corresponding AMMI estimates. The arguments used for this method
#'are:initial.values, precision, maxiter, change.factor and simplified.model
#'
#'\item EM-GGE: Iterative SVD-based imputation focusing on G+GE.
#' 
#' \item EM-SREG: Iterative algorithm using the Sites Regression model. Supports variants 
#'  like standard SVD and Bayesian PCA (EM-bSREG). 
#' 
#'\item Gabriel: combines regression and lower-rank approximation using SVD.
#'This method initially replaces the missing cells by arbitrary values, and
#'subsequently the imputations are refined through an iterative scheme that
#'defines a partition of the matrix for each missing value in turn and uses a
#'linear regression of columns (or rows) to obtain the new imputation. The
#'arguments used for this method is only the dataframe.
#'
#'\item WGabriel: is a a modification of Gabriel method that uses weights chosen
#'by cross-validation. The arguments used for this method are Winf and Wsup.
#'
#'\item EM-PCA: impute the missing entries of a mixed data using the iterative
#'PCA algorithm. The algorithm first consists imputing missing values with
#'initial values. The second step of the iterative PCA algorithm is to perform
#'PCA on the completed dataset to estimate the parameters. Then, it imputes the
#'missing values with the reconstruction formulae of order nPC (the fitted
#'matrix computed with nPC components for the scores and loadings). These steps
#'of estimation of the parameters via PCA and imputation of the missing values
#'using the fitted matrix are iterate until convergence. The arguments used for
#'this methods are: nPC, scale, method, row.w, coeff.ridge, precision, seed,
#'nb.init and maxiter
#'
#'}
#'
#'
#'@references Paderewski, J. (2013). \emph{An R function for imputation of missing
#'  cells in two-way data sets by EM-AMMI algorithm}. Communications in Biometry
#'  and Crop Science 8, 60–69.
#'  
#' @references 
#' Yan, W. (2013). \emph{Biplot analysis of incomplete two-way data}. 
#' Crop Science, 53(1), 48-57. \doi{10.2135/cropsci2012.05.0301}
#' 
#' Arciniegas-Alarcón, S., García-Peña, M., Krzanowski, W., & Dias, C. T. S. 
#' (2014b). \emph{An alternative methodology for imputing missing data in 
#' trials with genotype-by-environment interaction: some new aspects}. 
#' Biometrical Letters, 51(2), 75-88. \doi{10.2478/bile-2014-0006}
#' 
#' @references 
#' Angelini, J., Cervigni, G. D. L., & Quaglino, M. B. (2024). \emph{New 
#' imputation methodologies for genotype-by-environment data: an extensive 
#' study of properties of estimators}. Euphytica, 220(6), 92. 
#' \doi{10.1007/s10681-024-03344-z}
#' 
#'@references Julie Josse, Francois Husson (2016). missMDA: A Package for
#'  Handling Missing Values in Multivariate Data Analysis. Journal of
#'  Statistical Software 70, 1-31.
#'@references Arciniegas-Alarcón S., García-Peña M., Dias C.T.S., Krzanowski
#'  W.J. (2010). \emph{An alternative methodology for imputing missing data in
#'  trials with genotype-by-environment interaction}. Biometrical Letters 47,
#'  1–14.
#'@references Arciniegas-Alarcón S., García-Peña M., Krzanowski W.J., Dias
#'  C.T.S. (2014). \emph{An alternative methodology for imputing missing data in
#'  trials with genotype-byenvironment interaction: some new aspects.}
#'  Biometrical Letters 51, 75-88.
#'
#'@export
#'
#' @examples
#' library(geneticae)
#' # Data without replications
#' library(agridat)
#' data(yan.winterwheat)
#'
#' # generating missing values
#' yan.winterwheat[1,3]<-NA
#' yan.winterwheat[3,3]<-NA
#' yan.winterwheat[2,3]<-NA
#'
#' imputation(yan.winterwheat, genotype = "gen", environment = "env",
#'            response = "yield", type = "EM-AMMI")
#'
#' # Data with replications
#' data(plrv)
#' head(plrv)
#' plrv$Yield[plrv$Locality == "Ayac" & plrv$Rep %in% c(1, 2, 3) & plrv$Genotype == '102.18'] <- NA
#' 
#' imputation(plrv, nPC = 2,genotype = "Genotype", environment = "Locality", 
#'            response = "Yield", rep ='Rep', type = "EM-AMMI")
#'            
#' imputation(plrv, genotype = "Genotype", environment = "Locality", 
#'            response = "Yield", rep ='Rep', type = "EM-SREG")
#'
#'@importFrom stats var
#'@importFrom missMDA imputePCA
#'@importFrom dplyr group_by summarise rename %>% pull select all_of
#'@importFrom rlang sym
#'@importFrom tidyr pivot_wider
#'
imputation <- function(Data, genotype="gen", environment="env", response="yield", rep=NULL, 
                       type="EM-AMMI", nPC=2, initial.values=NA, precision=0.01, 
                       maxiter=1000, change.factor=1, simplified.model=FALSE,
                       scale = TRUE, method = "EM", row.w = NULL, coeff.ridge = 1, 
                       seed = NULL, nb.init = 1, Winf=0.8, Wsup=1) {
  
  if (missing(Data)) stop("Need to provide Data data frame")
  if (!any(is.na(Data))) stop("There are no missing data in input data frame")
  
  stopifnot(
    is.data.frame(Data),
    type %in% c("EM-AMMI", "EM-GGE", "EM-SREG", "EM-bSREG", "Gabriel", "Eigenvector", "WGabriel", "EM-PCA")
  )

  if(!is.null(rep)){
    Data_wide <- Data %>%
      dplyr::group_by(!!rlang::sym(genotype), !!rlang::sym(environment)) %>%
      dplyr::summarise(mean_resp = mean(!!rlang::sym(response), na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(environment), values_from = mean_resp) %>%
      as.data.frame()
  } else {
    Data_wide <- Data %>%
      tidyr::pivot_wider(names_from = !!rlang::sym(environment), values_from = !!rlang::sym(response)) %>%
      as.data.frame()
  }

  rownames(Data_wide) <- Data_wide[[genotype]]

  Data_mat <- Data_wide %>% dplyr::select(-dplyr::all_of(genotype))

  Data_mat <- as.matrix(Data_mat)
  Data_mat[is.na(Data_mat)] <- NA

  Data_mat <- as.data.frame(Data_mat)
  
  
  if(type == "EM-AMMI"){
    matrix_res <- EM.AMMI(Data_mat, PC.nb = nPC, initial.values = initial.values, 
                          precision = precision, max.iter = maxiter, 
                          change.factor = change.factor, 
                          simplified.model = simplified.model)$X
  }
  
  if(type == "EM-GGE"){
    matrix_res <- BiplotImputfun(Data_mat, precision = precision, 
                                 max.iter = maxiter, n_pc = nPC)$Biplot.imput
  }
  
  if(type %in% c("EM-SREG", "EM-bSREG")){
    matrix_res <- EM.SREG(Data_mat, PC.nb = nPC, initial.values = initial.values,
                          precision = precision, max.iter = maxiter, 
                          change.factor = change.factor, 
                          simplified.model = simplified.model, 
                          type = type)$X
  }
  
  if(type == "Gabriel"){
    matrix_res <- Gabriel.Calinski(Data_mat)$GabrielImput
  }
  
  else if(type == "Eigenvector"){
    matrix_res <- Eigenvectorfun(Data_mat, f = nPC)$EigenImputaciones
  }
  
  if(type == "WGabriel"){
    matrix_res <- WGabriel(Data_mat, Winf, Wsup)$GabrielWImput
  }
  
  if(type == "EM-PCA"){
    matrix_res <- missMDA::imputePCA(Data_mat, ncp = nPC, scale = scale, method = method,
                                     row.w = row.w, coeff.ridge = coeff.ridge, 
                                     threshold = precision, seed = seed, 
                                     nb.init = nb.init, maxiter = maxiter)$completeObs
  }
  
  return(as.matrix(matrix_res))
}
