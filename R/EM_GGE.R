#' Biplot Imputation Method
#'
#' @description This function implements the Biplot imputation method as 
#' proposed by Yan (2013). It is an iterative algorithm that uses the singular 
#' value decomposition (SVD) to impute missing values in a genotype by 
#' environment matrix.
#'
#' @param X A data frame or matrix with genotypes in rows and 
#'   environments in columns.
#' @param precision (optional) Convergence threshold. The algorithm stops 
#'   when the relative change in imputed values is less than this value. 
#'   Default is 0.01.
#' @param max.iter (optional) Maximum number of iterations. Default is 1000.
#' @param n_pc Number of principal components to use for imputation. 
#'   Default is 2.
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
#' @return A list containing:
#' \itemize{
#'   \item \code{X_imputed}: The final matrix with missing values filled.
#'   \item \code{iteration}: Total number of iterations performed.
#'   \item \code{convergence}: Final relative change reached.
#'   \item \code{fitted}: The fitted values from the AMMI/Biplot model.
#' }
#' 
#' @importFrom stats sd
#' @keywords internal
#' @export


BiplotImputfun <- function(X, precision = 0.01, max.iter = 1000, n_pc = 2) {
  
  X_mat <- as.matrix(X)
  indicamissing <- is.na(X_mat)
  indicaobservEM <- !indicamissing
  
  col_means <- colMeans(X_mat, na.rm = TRUE)
  X_original <- X_mat
  X_original[indicamissing] <- col_means[col(X_mat)[indicamissing]]
  
  iterconverg <- 1
  stabilitycrit <- 1 + precision
  
  RSS_N1 <- sqrt(sum(X_original[indicaobservEM]^2) / sum(indicaobservEM))
  
  while (stabilitycrit > precision && iterconverg <= max.iter) {
    
    C_original <- matrix(colMeans(X_original), nrow = nrow(X_original), ncol = ncol(X_original), byrow = TRUE)
    
    desvest_original <- matrix(apply(X_original, 2, stats::sd), 
                               nrow = nrow(X_original), ncol = ncol(X_original), byrow = TRUE)
    
    pij_original <- (X_original - C_original) / desvest_original
    
    dvs_pij <- svd(pij_original)
    
    U <- dvs_pij$u[, 1:n_pc, drop = FALSE]
    D <- diag(dvs_pij$d[1:n_pc], nrow = n_pc)
    V <- dvs_pij$v[, 1:n_pc, drop = FALSE]
    
    YanImputation <- U %*% D %*% t(V)
    
    Y_EMammi1 <- C_original + (YanImputation * desvest_original)
    
    RSS_N <- sqrt(sum((X_original[indicamissing] - Y_EMammi1[indicamissing])^2) / sum(indicamissing))
    stabilitycrit <- RSS_N / RSS_N1
    
    X_original[indicamissing] <- Y_EMammi1[indicamissing]
    
    iterconverg <- iterconverg + 1
  }
  
  return(list(
    Biplot.iter = iterconverg - 1,
    Biplot.converg = stabilitycrit,
    Biplot.imput = X_original,
    X2.Chapeu = Y_EMammi1
  ))
}