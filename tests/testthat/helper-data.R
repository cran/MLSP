# Simulated VNIR data: additive and multiplicative scatter plus an absorption
# feature at 1700 nm whose depth depends on soil organic carbon (SOC).
sim_data <- function(n = 40, seed = 1) {
  set.seed(seed)
  wl <- 350:2500
  soc <- runif(n, 0.5, 4)
  offset <- runif(n, -0.05, 0.05)
  mult <- runif(n, 0.9, 1.1)
  shape <- 0.35 + 0.15 * (wl - 350) / 2150 - 0.04 * exp(-((wl - 1930) / 40)^2)
  spec <- t(sapply(seq_len(n), function(i)
    offset[i] + mult[i] * (shape - 0.02 * soc[i] * exp(-((wl - 1700) / 50)^2)) +
      stats::rnorm(length(wl), 0, 5e-4)))
  colnames(spec) <- wl
  list(spectra = data.frame(Wavelength = seq_len(n), spec, check.names = FALSE),
       soil = data.frame(LAB_NUM = seq_len(n), SOC = soc, clay = runif(n, 5, 40)),
       vnir = as.data.frame(spec, check.names = FALSE))
}
