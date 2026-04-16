library(testthat)
library(geneticae)
library(ggplot2)

# --- FUNCIÓN AUXILIAR PARA UN MODELO ESTABLE ---
generar_modelo_gge <- function() {
  set.seed(123)
  df <- expand.grid(Genotype = paste0("G", 1:6), 
                    Locality = paste0("L", 1:4)) %>%
    mutate(Yield = 10 + rnorm(n()))
  # rSREGModel es el motor para GGE
  return(rSREGModel(df, "Genotype", "Locality", "Yield"))
}

test_that("rSREGPlot genera todos los tipos de biplot sin errores", {
  gge_obj <- generar_modelo_gge()
  
  # Tipos que no requieren parámetros adicionales
  tipos_simples <- c("Biplot", "Relationship Among Environments", 
                     "Which Won Where/What", "Discrimination vs. representativeness",
                     "Ranking Environments", "Mean vs. Stability", "Ranking Genotypes")
  
  for (t in tipos_simples) {
    p <- rSREGPlot(gge_obj, type = t)
    expect_s3_class(p, "ggplot")
    expect_true(length(p$layers) > 0)
  }
}) #