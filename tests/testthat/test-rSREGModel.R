test_that("Validación de entradas y errores iniciales", {
  # Error por falta de Data
  expect_error(rSREGModel(genotype="gen"), "Need to provide Data data frame")
  
  # Error por presencia de NAs
  dat_na <- yan.winterwheat
  dat_na[1, "yield"] <- NA
  expect_error(rSREGModel(dat_na, "gen", "env", "yield"), "Missing data in input data frame")
  
  # Error por argumentos de clase incorrecta
  expect_error(rSREGModel(yan.winterwheat, genotype = 123, environment="env", response="yield"), "genotype")
  expect_error(rSREGModel(yan.winterwheat, "gen", "env", "yield", model = "Inexistente"), "model")
})

test_that("Estructura y clase del objeto de salida", {
  model_sreg <- rSREGModel(yan.winterwheat, "gen", "env", "yield")
  expect_s3_class(model_sreg, "rSREGModel")
  expect_named(model_sreg, c("model", "coordgenotype", "coordenviroment", 
                             "eigenvalues", "vartotal", "varexpl", "labelgen", 
                             "labelenv", "labelaxes", "Data", "SVP"))
})

test_that("Partición de Valores Singulares (SVP) y Varianza", {
  m_sym <- rSREGModel(yan.winterwheat, "gen", "env", "yield", SVP = "symmetrical")
  m_row <- rSREGModel(yan.winterwheat, "gen", "env", "yield", SVP = "row")
  
  expect_equal(m_sym$eigenvalues, m_row$eigenvalues)
  expect_false(identical(m_sym$coordgenotype, m_row$coordgenotype))
})