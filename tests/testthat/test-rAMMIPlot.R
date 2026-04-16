library(testthat)
library(ggplot2)
library(geneticae)

# --- FUNCIÓN AUXILIAR ESTABLE ---
generar_modelo_test <- function() {
  # Datos con estructura aditiva clara para facilitar convergencia
  set.seed(123)
  df <- expand.grid(Genotype = paste0("G", 1:6), 
                    Locality = paste0("L", 1:4)) %>%
    mutate(Yield = 10 + 
             as.numeric(as.factor(Genotype)) * 0.5 + 
             as.numeric(as.factor(Locality)) * 1.5 + 
             rnorm(n(), sd = 0.1))
  
  # Usamos AMMI para los tests de interfaz (es el más rápido y estable)
  return(rAMMIModel(df, "Genotype", "Locality", "Yield", type = "AMMI"))
}

# --- TESTS ---

test_that("rAMMIPlot genera un objeto ggplot válido", {
  modelo <- generar_modelo_test()
  p <- rAMMIPlot(modelo)
  expect_s3_class(p, "ggplot")
})

test_that("rAMMIPlot respeta los argumentos de personalización", {
  modelo <- generar_modelo_test()
  p_custom <- rAMMIPlot(modelo, titles = FALSE, axislabels = FALSE, footnote = FALSE)
  
  expect_null(p_custom$labels$title)
  expect_null(p_custom$labels$caption)
  
  # Verificamos que no haya porcentajes si axislabels = FALSE
  label_x <- p_custom$labels$x
  if (!is.null(label_x)) expect_false(grepl("%", label_x))
})

test_that("rAMMIPlot muestra porcentajes cuando axislabels es TRUE", {
  modelo <- generar_modelo_test()
  p <- rAMMIPlot(modelo, axislabels = TRUE)
  
  expect_true(grepl("%", p$labels$x))
  expect_true(grepl("%", p$labels$y))
})

test_that("rAMMIPlot contiene las capas esenciales", {
  modelo <- generar_modelo_test()
  p <- rAMMIPlot(modelo)
  
  capas <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_true(any(capas == "GeomSegment"))
  expect_true(any(capas == "GeomText"))
})