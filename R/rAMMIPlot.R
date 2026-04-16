#' AMMI Biplots with ggplot2
#'
#' @description Produces a biplot for objects of class 'AMMI'.
#'
#' @param model_res an object of class 'AMMI' from AMMIModel.
#' @param colGen genotype colour. Defaults to "gray47".
#' @param colEnv environment colour. Defaults to "darkred".
#' @param sizeGen genotype text size.
#' @param sizeEnv environment text size.
#' @param titles logical, show plot title.
#' @param footnote logical, show footnote with explained variance.
#' @param axis_expand expansion factor for axis limits.
#' @param limits logical. If `TRUE` axes are automatically rescaled. Defaults to
#'  `TRUE`.
#' @param axes logical, if this argument is `TRUE` axes passing through the
#'  origin are drawn. Defaults to `TRUE`.
#' @param axislabels logical, if this argument is `TRUE` labels axes are included.
#'  Defaults to `TRUE`
#'
#' @import ggplot2
#' @export
rAMMIPlot <- function(model_res, colGen = "gray47", colEnv = "darkred", sizeGen = 6, sizeEnv = 6, 
                      titles = TRUE, footnote = TRUE, axis_expand = 1.2, limits = TRUE, 
                      axes = TRUE, axislabels = TRUE) {
  
  if (!inherits(model_res, "rAMMI")) stop("Se requiere un objeto de clase 'rAMMI'")
  type <- label <- d1 <- d2 <- x0 <- y0 <- radio <- NULL

  scaling <- 0.5
  idx <- 1:model_res$Ncomp
  lambda_scaling <- model_res$eigenvalues[idx]^scaling
  
  coordgenotype <- t(t(model_res$gen_scores[, idx]) / lambda_scaling)
  coordenviroment <- t(t(model_res$env_scores[, idx]) * lambda_scaling)
  
  varexpl <- round((model_res$eigenvalues[idx]^2 / model_res$vartotal) * 100, 2)
  
  plotdata <- data.frame(
    rbind(
      data.frame(coordgenotype, type = "genotype", label = model_res$gen_labels),
      data.frame(coordenviroment, type = "environment", label = model_res$env_labels)
    )
  )
  colnames(plotdata)[1:2] <- c("Component1", "Component2")
  
  p <- ggplot(data = plotdata, aes(x = Component1, y = Component2)) +
    theme_classic()
  
  if (axes) {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
  }
  
  p <- p + 
    geom_segment(data = subset(plotdata, type == "environment"),
                 aes(x = 0, y = 0, xend = Component1, yend = Component2),
                 arrow = arrow(length = unit(0.2, "cm")), col = colEnv, alpha = 0.7) +
    geom_text(data = subset(plotdata, type == "environment"),
              aes(label = label), col = colEnv, size = sizeEnv, vjust = -0.5) +
    geom_text(data = subset(plotdata, type == "genotype"),
              aes(label = label), col = colGen, size = sizeGen)
  
  if (titles) p <- p + ggtitle(paste("AMMI Biplot - Model:", model_res$type))
  
  if (axislabels) {
    p <- p + xlab(paste0("PC1 (", varexpl[1], "%)")) +
      ylab(paste0("PC2 (", varexpl[2], "%)"))
  }
  
  if (footnote) {
    foot_txt <- paste(model_res$type, "biplot explaining", sum(varexpl), "% of GxE variation")
    p <- p + labs(caption = foot_txt)
  }
  
  if (limits) p <- p + coord_fixed(ratio = 1)
  
  return(p)
}