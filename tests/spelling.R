if(requireNamespace('spelling', quietly = TRUE))
  spelling::spell_check_test(vignettes = "/media/julia/KINGSTON/Tesis doctorado en estadistica/Actualizacion geneticae/geneticae/", error = FALSE,
                             skip_on_cran = TRUE)
