test_that("PCR and PLSR choose components by CV and give valid metrics", {
  d <- sim_data(n = 40)
  xy <- merge_of_lab_and_spectrum(d$soil, d$spectra)
  for (f in list(pcr_preprocess, plsr_preprocess)) {
    fit <- f(xy$soil, xy$vnir.matrix, xy$j, preprocess = "snv",
             type_of_soil = "SOC", ncomp_max = 10)
    expect_length(fit, 4)
    lv <- vapply(fit, attr, numeric(1), "LV")
    expect_true(all(lv >= 1 & lv <= 10))
    res <- results(fit)
    expect_true(res[["cv-r2"]] < 1)          # no longer a perfect, invalid fit
    expect_gt(res[["val-r2"]], 0.8)          # the SOC feature is recovered
  }
})

test_that("predictions have one value per sample", {
  d <- sim_data(n = 30)
  xy <- merge_of_lab_and_spectrum(d$soil, d$spectra)
  X <- as.matrix(spectral_preprocess(xy$vnir.matrix, "raw"))
  cal <- xy$j[[1]]
  for (m in c("pcr", "plsr", "lasso", "cubist", "rf")) {
    fp <- MLSP:::.fit_predict(m, X[cal, ], xy$soil$SOC[cal], X[-cal, ],
                              ncomp_max = 5, rf_ntree = 50)
    expect_length(fp$cal, length(cal))
    expect_length(fp$val, nrow(X) - length(cal))
  }
})

test_that("msd.comp handles constant predictions without warnings", {
  expect_silent(m <- msd.comp(c(1, 2, 3, 4), rep(2.5, 4)))
  expect_equal(unname(m["r2"]), 0)
})

test_that("ml_f ranks combinations and runs all 26 methods", {
  d <- sim_data(n = 30)
  top <- ml_f(d$soil, d$spectra, smoother_selection = c("raw", "snv"),
              type_of_soil = "SOC", models = "plsr", top_n = 1)
  expect_equal(nrow(top), 1)
  expect_equal(top$preprocessing, "snv")

  all26 <- ml_f(d$soil, d$spectra, type_of_soil = 2, models = "plsr",
                model_selection = FALSE, ncomp_max = 5)
  expect_equal(nrow(all26), 26)
  expect_setequal(all26$preprocessing, preprocess_methods()$method)
  expect_false(anyNA(all26[["val-rpiq"]]))
  expect_error(ml_f(d$soil, d$spectra, smoother_selection = "sg.c.11.0",
                    type_of_soil = 2, models = "plsr"), "sg.c.11.0 -> sg.q.11.0")
})
