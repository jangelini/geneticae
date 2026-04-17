#' Robust AMMI Model Fitting
#'
#' @description Fits a classical or robust Additive Main effects and Multiplicative 
#' Interaction (AMMI) model for genotype-by-environment data.
#'
#' @param Data a dataframe with genotypes, environments, repetitions (if any) and 
#'   the phenotypic trait of interest. Other variables that will not be used in 
#'   the analysis can be included.
#' @param genotype column name containing genotypes. Defaults to `"gen"`.
#' @param environment column name containing environments. Defaults to `"env"`.
#' @param response column name containing the phenotypic trait of interest. Defaults to `"Y"`.
#' @param rep column name containing replications. If this argument is `NULL` 
#'   (default), it is assumed that the data already contains means per genotype 
#'   in each environment. If provided, means are calculated automatically.
#' @param Ncomp number of principal components to retain for the interaction part. 
#'   Defaults to 2.
#' @param type method for fitting the AMMI model: `"AMMI"` (classical), `"rAMMI"`, 
#'   `"hAMMI"`, `"gAMMI"`, `"lAMMI"` or `"ppAMMI"` (robust variants). 
#'   Defaults to `"AMMI"`.
#'
#' @details 
#' To overcome the problem of data contamination with outlying observations, 
#' Rodrigues, Monteiro and Lourenco (2015) propose a robust AMMI model based on 
#' the M-Huber estimator and robust SVD/PCA procedures. 
#' 
#' The `type` argument allows choosing between several robust strategies:
#' \itemize{
#'   \item \bold{AMMI}: Classical AMMI model using Least Squares and standard SVD.
#'   \item \bold{rAMMI}: Uses the L1 norm instead of the L2 norm to compute a 
#'     robust approximation to the SVD (via \pkg{pcaMethods}).
#'   \item \bold{hAMMI}: Uses the Hubert's approach (PcaHubert) combining 
#'     projection-pursuit and robust covariance estimation.
#'   \item \bold{gAMMI}: Uses the Grid search algorithm for PCA (PcaGrid).
#'   \item \bold{lAMMI}: Performs PCA on the data projected onto a unit sphere (PcaLocantore).
#'   \item \bold{ppAMMI}: Uses projection-pursuit (PcaProj) to calculate robust 
#'     eigenvalues and eigenvectors.
#' }
#' 
#' @return A list of class \code{rAMMI} containing:
#' \item{gen_scores}{Matrix of genotype scores (U * D).}
#' \item{env_scores}{Matrix of environment loadings (V).}
#' \item{eigenvalues}{Vector of singular values for the retained components.}
#' \item{gen_labels}{Names of the genotypes.}
#' \item{env_labels}{Names of the environments.}
#' \item{Ncomp}{Number of principal components used.}
#' \item{type}{The fitting method used.}
#' \item{vartotal}{Total variance explained by the multiplicative terms.}
#'
#' @references Rodrigues P.C., Monteiro A., Lourenco V.M. (2015). \emph{A robust 
#'   AMMI model for the analysis of genotype-by-environment data}. Bioinformatics 
#'   32, 58–66.
#'
#' @importFrom MASS rlm
#' @importFrom pcaMethods robustSvd
#' @importFrom rrcov PcaHubert PcaGrid PcaLocantore PcaProj
#' @importFrom stats lm residuals as.formula
#' @importFrom dplyr group_by summarise pull arrange %>%
#' @importFrom rlang sym
#' @export
#' 
#' @examples
#' library(agridat)
#' data(yan.winterwheat)
#' 
#' # Classical AMMI
#' mod_ammi <- rAMMIModel(yan.winterwheat, genotype = "gen", 
#'                        environment = "env", response = "yield", type = "AMMI")
#' 
#' # Robust AMMI (using Hubert's method)
#' mod_rammi <- rAMMIModel(yan.winterwheat, genotype = "gen", 
#'                         environment = "env", response = "yield", type = "hAMMI")
#' 
#' 
rAMMIModel <- function(Data, genotype = "gen", environment = "env", response = "Y", rep = NULL, Ncomp = 2, type = "AMMI") {
  
  if (missing(Data)) stop("Need to provide Data data frame")
  if (any(is.na(Data))) stop("Missing data in input data frame")
  
  # --- PROCESAMIENTO DE DATOS ---
  Data_avg <- Data %>%
    group_by(!!sym(genotype), !!sym(environment)) %>%
    summarise(y = mean(!!sym(response), na.rm = TRUE), .groups = "drop") %>%
    arrange(!!sym(environment), !!sym(genotype))
  
  gen <- as.factor(Data_avg[[genotype]])
  env <- as.factor(Data_avg[[environment]])
  y   <- Data_avg$y
  
  Ngen <- nlevels(gen)
  Nenv <- nlevels(env)

  
  if (type == "AMMI") {
    lm.x <- lm(y ~ gen + env, data = Data_avg)
    residuals.x <- matrix(residuals(lm.x), nrow = Ngen, ncol = Nenv)
  }

  if (type != "AMMI") {
    rlm.x <- MASS::rlm(y ~ gen + env, data = Data_avg, maxit = 100) # Aumentamos maxit por seguridad
    rresiduals.x <- matrix(residuals(rlm.x), nrow = Ngen, ncol = Nenv)
  }
  
  if (type == "AMMI") {
    svd.x <- svd(residuals.x)
    singlevalue <- svd.x$d
    loading <- svd.x$v
    score <- svd.x$u %*% diag(svd.x$d)
  } else if (type == "rAMMI") {
    svd_result <- tryCatch({
      robustSvd(rresiduals.x)
    }, error = function(e) {
      stop(paste("The robustSvd algorithm did not converge with this data.",
                 "Try using type='AMMI' or rrcov methods like 'hAMMI'."))
    })
    singlevalue <- svd_result$d
    loading <- svd_result$v
    score <- svd_result$u %*% diag(svd_result$d)
  } else {
    res_pca <- switch(type,
                      "hAMMI"  = PcaHubert(rresiduals.x, mcd = FALSE),
                      "gAMMI"  = PcaGrid(rresiduals.x),
                      "lAMMI"  = PcaLocantore(rresiduals.x),
                      "ppAMMI" = PcaProj(rresiduals.x))
    
    singlevalue <- if(type == "lAMMI") res_pca@eigenvalues else sqrt(res_pca@eigenvalues)
    loading     <- res_pca@loadings
    score       <- res_pca@scores
    singlevalue <- singlevalue * sqrt(nrow(score))
  }
  
  res <- list(
    gen_scores = score,
    env_scores = loading,
    eigenvalues = singlevalue,
    gen_labels = levels(gen),
    env_labels = levels(env),
    Ncomp = Ncomp,
    type = type,
    vartotal = sum(singlevalue^2)
  )
  class(res) <- "rAMMI"
  return(res)
}