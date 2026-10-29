test_that("the catalogue has raw, SNV, log(1/R) and 23 Savitzky-Golay filters", {
  pm <- preprocess_methods()
  expect_equal(nrow(pm), 26)
  expect_false(anyDuplicated(pm$method) > 0)
  sg <- pm[pm$family == "Savitzky-Golay", ]
  expect_equal(nrow(sg), 23)
  expect_equal(as.vector(table(sg$derivative)), c(9, 9, 5))
  expect_true(all(c("sg.q.11.0", "sg.q.11.1", "sg.c.11.1", "sg.q.11.2") %in% sg$method))
})

test_that("all 26 preprocessing methods run and return 202 finite bands", {
  d <- sim_data(n = 5)
  for (m in preprocess_methods()$method) {
    out <- spectral_preprocess(d$vnir, m)
    expect_equal(dim(out), c(5L, 202L), info = m)
    expect_true(all(is.finite(as.matrix(out))), info = m)
  }
})

test_that("all 26 methods give distinct spectra", {
  d <- sim_data(n = 3)
  pm <- preprocess_methods()$method
  outs <- lapply(pm, function(m) as.matrix(spectral_preprocess(d$vnir, m)))
  for (a in 1:(length(pm) - 1)) for (b in (a + 1):length(pm)) {
    expect_gt(max(abs(outs[[a]] - outs[[b]])), 1e-8, label = paste(pm[a], pm[b]))
  }
})

test_that("excluded settings are identical to the method that replaces them", {
  d <- sim_data(n = 3)
  X <- unname(as.matrix(d$vnir))
  agg <- function(Y) unname(colsMean(Y[, 7:2146], 10)[, -c(1:4, 64:66, 144:146, 213:214)])
  sgf <- function(p, n, m) agg(t(gsignal::sgolayfilt(t(X), p = p, n = n, m = m)))
  for (n in c(5, 7, 9, 11)) {
    expect_equal(sgf(3, n, 0), sgf(2, n, 0), tolerance = 1e-10)   # cubic = quadratic smoother
    expect_equal(sgf(3, n, 2), sgf(2, n, 2), tolerance = 1e-10)   # cubic = quadratic 2nd derivative
  }
  for (n in c(3, 5, 7, 9, 11)) {
    expect_equal(sgf(1, n, 1), sgf(2, n, 1), tolerance = 1e-10)   # linear = quadratic 1st derivative
  }
  expect_equal(sgf(2, 3, 0), agg(X), tolerance = 1e-10)           # 3-point quadratic = raw
  expect_error(spectral_preprocess(d$vnir, "sg.c.11.0"), "sg.q.11.0")
  expect_error(spectral_preprocess(d$vnir, "sg.l.7.1"), "sg.q.7.1")
  expect_error(spectral_preprocess(d$vnir, "sg.q.3.0"), "raw")
})

test_that("Savitzky-Golay works on named spectra (gsignal names bug)", {
  d <- sim_data(n = 3)
  expect_false(is.null(colnames(d$vnir)))
  expect_silent(spectral_preprocess(d$vnir, "sg.q.11.1"))
})

test_that("9- and 11-point quadratic smoothers differ from the 7-point one", {
  d <- sim_data(n = 3)
  q7 <- as.matrix(spectral_preprocess(d$vnir, "sg.q.7.0"))
  q9 <- as.matrix(spectral_preprocess(d$vnir, "sg.q.9.0"))
  q11 <- as.matrix(spectral_preprocess(d$vnir, "sg.q.11.0"))
  expect_gt(max(abs(q9 - q7)), 0)
  expect_gt(max(abs(q11 - q9)), 0)
  X <- unname(as.matrix(d$vnir))
  direct <- t(gsignal::sgolayfilt(t(X), p = 2, n = 11))
  expect_equal(unname(q11), unname(colsMean(direct[, 7:2146], 10)[,
    -c(1:4, 64:66, 144:146, 213:214)]))
})

test_that("raw() equals spectral_preprocess(, 'raw') and input is validated", {
  d <- sim_data(n = 3)
  expect_equal(raw(d$vnir), spectral_preprocess(d$vnir, "raw"))
  expect_error(spectral_preprocess(d$vnir[, 1:100], "raw"), "2151")
  expect_error(spectral_preprocess(d$vnir, "sg.x.3.0"), "Unknown")
})
