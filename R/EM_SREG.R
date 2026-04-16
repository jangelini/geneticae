#' EM-SREG Imputation Method
#'
#' @description Iterative algorithm to impute missing values in two-way tables 
#'   using the Sites Regression (SREG) model. It supports several variants 
#'   including standard SVD and Bayesian PCA.
#'
#' @param X A data frame or matrix with genotypes in rows and environments in columns.
#' @param PC.nb Number of principal components to be used. Default is 1.
#' @param initial.values (optional) Initial values for missing cells. If NA, 
#'   initial values are obtained from column means (environment effects).
#' @param precision Convergence threshold. Default is 0.01.
#' @param max.iter Maximum number of iterations. Default is 1000.
#' @param change.factor Step size for updating missing values (standard is 1).
#' @param simplified.model Logical. If TRUE, effects are only calculated in the first iteration.
#' @param type Method type: "EM-SREG" (Standard), "EM-bSREG" (Bayesian).
#'
#' @references 
#' Angelini, J., Cervigni, G. D. L., & Quaglino, M. B. (2024). \emph{New 
#' imputation methodologies for genotype-by-environment data: an extensive 
#' study of properties of estimators}. Euphytica, 220(6), 92. 
#' \doi{10.1007/s10681-024-03344-z}
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{X}: The imputed matrix.
#'   \item \code{iter}: The number of iterations until convergence.
#' }
#' 
#' @importFrom stats sd
#' @importFrom pcaMethods pca scores loadings
#' @keywords internal
#' @export
EM.SREG <- function(X, PC.nb = 1, initial.values = NA, precision = 0.01,
                    max.iter = 1000, change.factor = 1, simplified.model = FALSE, 
                    type = c("EM-SREG", "EM-bSREG")) {
  
  type <- match.arg(type)
  X <- as.matrix(X)
  indicamissing <- is.na(X)
  
  # Determinar número máximo de componentes posibles
  X_missing_mask <- !indicamissing
  max.IPC <- min(rowSums(X_missing_mask), colSums(X_missing_mask) - 1)
  PC.nb.used <- if (max.IPC < PC.nb) max.IPC else PC.nb
  
  # --- 1. Inicialización ---
  X.new <- X
  if (length(initial.values) == 1) {
    if (!is.na(initial.values)) {
      X.new[indicamissing] <- initial.values
    } else {
      # Efecto de columna (Ambiente) como inicialización para SREG
      col_m <- colMeans(X, na.rm = TRUE)
      X.new[indicamissing] <- col_m[col(X)[indicamissing]]
    }
  } else {
    X.new[indicamissing] <- initial.values[indicamissing]
  }
  
  iteration <- 1
  change <- precision + 1
  
  # --- 2. Algoritmo Iterativo ---
  while ((change > precision) && (iteration < max.iter)) {
    
    if (iteration == 1 || !simplified.model) {
      x.mean <- mean(X.new)
      # Centrado por columna (SREG se centra en el ambiente)
      X_centered <- scale(X.new - x.mean, center = TRUE, scale = FALSE)
      col_eff <- attr(X_centered, "scaled:center")
    } else {
      X_centered <- X.new - x.mean - matrix(col_eff, nrow(X.new), ncol(X.new), byrow = TRUE)
    }
    
    # --- Ajuste según tipo de SREG ---
    if (PC.nb.used >= 1) {
      if (type == "EM-SREG") {
        # SVD estándar
        SVD <- svd(X_centered)
        interaction.adj <- SVD$u[, 1:PC.nb.used, drop = FALSE] %*% 
          diag(SVD$d[1:PC.nb.used], nrow = PC.nb.used) %*% 
          t(SVD$v[, 1:PC.nb.used, drop = FALSE])
      } else if (type == "EM-bSREG") {
        # Requiere paquete pcaMethods
        if (!requireNamespace("pcaMethods", quietly = TRUE)) stop("Package 'pcaMethods' needed for EM-bSREG.")
        Xpca <- pcaMethods::pca(X_centered, method = "bpca", nPcs = PC.nb.used, center = FALSE)
        interaction.adj <- pcaMethods::scores(Xpca) %*% t(pcaMethods::loadings(Xpca))
      } else {
        # Placeholder para los otros tipos (gSREG, ppSREG, etc.)
        interaction.adj <- 0 
      }
    } else {
      interaction.adj <- 0
    }
    
    # Reconstrucción del modelo SREG: Media + Efecto Ambiente + Interacción
    X.fitted <- x.mean + matrix(col_eff, nrow(X), ncol(X), byrow = TRUE) + interaction.adj
    
    # Calcular cambio solo en los faltantes
    X.next <- X.new
    X.next[indicamissing] <- X.fitted[indicamissing]
    
    change <- max(abs(X.next[indicamissing] - X.new[indicamissing]))
    
    # Aplicar factor de cambio y actualizar
    X.new <- change.factor * X.next + (1 - change.factor) * X.new
    iteration <- iteration + 1
  }
  
  return(list(X = X.new, iter = iteration))
}