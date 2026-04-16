#' Robust AMMI Model
#'
#' @description Fits a classical or robust Additive Main effects and Multiplicative Interaction (AMMI) model.
#'
#' @param Data a dataframe with genotypes, environments and the phenotypic trait.
#' @param genotype column name containing genotypes.
#' @param environment column name containing environments.
#' @param response column name containing the phenotypic trait.
#' @param rep column name containing replications. If provided, means are calculated.
#' @param Ncomp number of principal components to retain.
#' @param type method for fitting: `"AMMI"`, `"rAMMI"`, `"hAMMI"`, `"gAMMI"`, `"lAMMI"` or `"ppAMMI"`.
#'
#' @importFrom MASS rlm
#' @importFrom pcaMethods robustSvd
#' @importFrom rrcov PcaHubert PcaGrid PcaLocantore PcaProj
#' @importFrom stats lm residuals as.formula
#' @importFrom dplyr group_by summarise pull arrange %>%
#' @importFrom rlang sym
#' @export
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