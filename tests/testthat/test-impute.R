library(geneticae)
library(agridat)
library(dplyr)

# 1. Preparación de datos
data(yan.winterwheat)
data(plrv)

# Dataset con NAs
dat_na <- yan.winterwheat
dat_na[1, 3] <- NA
dat_na[5, 3] <- NA

test_that("Validación de errores de entrada en imputación", {
  expect_error(imputation(yan.winterwheat), "There are no missing data in input data frame")
  expect_error(imputation(genotype="gen"), "Need to provide Data data frame")
  expect_error(imputation(dat_na, type = "METODO_FALSO"))
})

test_that("Consistencia de los métodos de imputación principales", {
  res_ammi <- imputation(dat_na, genotype="gen", environment="env", response="yield")
  expect_true(is.matrix(res_ammi))
  expect_false(any(is.na(res_ammi)))
  expect_equal(nrow(res_ammi), length(unique(yan.winterwheat$gen)))
})

test_that("Prueba de métodos específicos (EM-PCA, Gabriel, GGE)", {
  # EM-PCA
  expect_false(any(is.na(imputation(dat_na, type = "EM-PCA", nPC = 1))))
  
  # Gabriel
  expect_false(any(is.na(imputation(dat_na, type = "Gabriel"))))
  
  # EM-GGE
  expect_false(any(is.na(imputation(dat_na, type = "EM-GGE"))))
})

test_that("Manejo de repeticiones y balanceo", {
  # Usamos suppressWarnings para el aviso interno de max() -Inf en la convergencia
  res_rep <- suppressWarnings(
    imputation(plrv %>% mutate(Yield = ifelse(row_number() == 2, NA, Yield)), 
               genotype = "Genotype", environment = "Locality", 
               response = "Yield", rep = "Rep", type = "EM-AMMI")
  )
  
  expect_true(is.matrix(res_rep))
  expect_equal(nrow(res_rep), length(unique(plrv$Genotype)))
  expect_false(any(is.na(res_rep)))
})

test_that("Argumentos técnicos e integridad", {
  expect_silent(suppressWarnings(imputation(dat_na, nPC = 1, precision = 0.1)))
  
  # Igualdad por defecto
  expect_equal(imputation(dat_na), imputation(dat_na, type = "EM-AMMI"))
})