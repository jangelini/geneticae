library(testthat)
library(dplyr)
library(geneticae)

# Función auxiliar para generar datos limpios
get_simple_test_data <- function() {
  set.seed(123)
  expand.grid(
    Genotype = paste0("G", 1:8), 
    Locality = paste0("L", 1:4), 
    Rep = 1:2
  ) %>%
    mutate(Yield = 10 + rnorm(n()))
}

test_that("rAMMIModel procesa correctamente los datos y promedios", {
  df <- get_simple_test_data()
  
  # Probamos el modelo base (AMMI clásico)
  # Esto valida que dplyr, group_by, summarise y arrange funcionen bien
  res <- rAMMIModel(df, "Genotype", "Locality", "Yield", rep = "Rep", type = "AMMI")
  
  expect_s3_class(res, "rAMMI")
  expect_equal(nrow(res$gen_scores), 8)
  expect_equal(nrow(res$env_scores), 4)
  expect_equal(res$type, "AMMI")
})

test_that("rAMMIModel funciona con métodos estables de rrcov", {
  df <- get_simple_test_data()
  
  # Testeamos hAMMI y gAMMI que suelen ser más estables que robustSvd
  estables <- c("hAMMI", "gAMMI")
  
  for (t in estables) {
    expect_error(
      rAMMIModel(df, "Genotype", "Locality", "Yield", type = t),
      NA,
      info = paste("Falló el método estable:", t)
    )
  }
})

test_that("rAMMIModel lanza error con NAs", {
  df <- get_simple_test_data()
  df$Yield[1] <- NA
  
  expect_error(
    rAMMIModel(df, "Genotype", "Locality", "Yield"),
    "Missing data in input data frame"
  )
})