#' Merge Soil Laboratory Data with Spectral Data
#'
#' Averages replicate scans of each sample, matches the spectra to the soil
#' laboratory data, and creates repeated random calibration/validation splits
#' used by every model in \code{\link{ml_f}}.
#'
#' @param soil_data A data frame of soil laboratory data with a sample
#'   identifier column (default \code{LAB_NUM}) and one column per soil property.
#' @param data_NaturaSpec_cleaned A data frame of spectra with a sample
#'   identifier column (default \code{Wavelength}, matching \code{LAB_NUM})
#'   and 2151 reflectance columns (350-2500 nm at 1 nm). Several rows (scans)
#'   per sample are averaged.
#' @param soil_id Name of the identifier column in \code{soil_data}.
#' @param spec_id Name of the identifier column in \code{data_NaturaSpec_cleaned}.
#' @param n_splits Number of random calibration/validation splits.
#' @param cal_prop Proportion of samples in each calibration set.
#'
#' @details
#' Splits are drawn with \code{\link[base]{sample}}; call
#' \code{\link[base]{set.seed}} beforehand for reproducible results. All
#' preprocessing and calibration methods are evaluated on the same splits, so
#' their accuracy metrics are directly comparable.
#'
#' @return A named list:
#' \describe{
#'   \item{soil}{Soil laboratory data for the matched samples (identifier first).}
#'   \item{vnir.matrix}{Matrix of averaged reflectance (2151 columns).}
#'   \item{j}{List of \code{n_splits} vectors of calibration-row indices.}
#'   \item{rm1, rm2, rm3, rm4}{Indices of the 1-nm bands excluded as edges or
#'     detector artifacts (350-395, 986-1015, 1786-1815, and 2476-2500 nm).}
#'   \item{ind}{Indices of the 1-nm bands averaged into 10-nm bands.}
#'   \item{rmove}{Indices of the 10-nm bands removed after averaging.}
#'   \item{vars}{Names of the 202 retained 10-nm bands.}
#' }
#'
#' @examples
#' # Simulated data: 50 soil samples scanned at 350-2500 nm (2151 bands).
#' set.seed(1)
#' n <- 50
#' wl <- 350:2500
#' soc <- runif(n, 0.5, 4)
#' spec <- t(sapply(soc, function(s)
#'   0.5 - 0.05 * s * exp(-((wl - 1400) / 600)^2) + rnorm(length(wl), 0, 0.01)))
#' colnames(spec) <- wl
#' spectra <- data.frame(Wavelength = 1:n, spec, check.names = FALSE)
#' soil_data <- data.frame(LAB_NUM = 1:n, SOC = soc, clay = runif(n, 5, 40),
#'   pH = runif(n, 5, 8))
#' merged <- merge_of_lab_and_spectrum(soil_data, spectra)
#' str(merged[1:3], max.level = 1)
#'
#' @importFrom stats aggregate
#'
#' @export
merge_of_lab_and_spectrum <- function(soil_data, data_NaturaSpec_cleaned,
                                      soil_id = "LAB_NUM", spec_id = "Wavelength",
                                      n_splits = 4, cal_prop = 0.7) {
  spectra <- data_NaturaSpec_cleaned
  if (!soil_id %in% names(soil_data)) {
    stop(sprintf("'soil_data' has no identifier column '%s'.", soil_id), call. = FALSE)
  }
  if (!spec_id %in% names(spectra)) {
    stop(sprintf("The spectral data have no identifier column '%s'.", spec_id), call. = FALSE)
  }
  band_cols <- setdiff(names(spectra), spec_id)
  if (length(band_cols) != 2151) {
    stop(sprintf(paste0("The spectral data must have 2151 reflectance columns ",
                        "(350-2500 nm) besides '%s'; found %d."),
                 spec_id, length(band_cols)), call. = FALSE)
  }

  # Average replicate scans per sample
  spec_mean <- aggregate(spectra[, band_cols, drop = FALSE],
                         by = list(.id = spectra[[spec_id]]), FUN = mean)
  names(spec_mean)[1] <- soil_id

  merged <- merge(soil_data, spec_mean, by = soil_id, sort = TRUE)
  if (nrow(merged) < 10) {
    stop(sprintf("Only %d samples matched between the soil and spectral data.",
                 nrow(merged)), call. = FALSE)
  }
  soil <- merged[, names(soil_data), drop = FALSE]
  vnir.matrix <- merged[, band_cols, drop = FALSE]

  n.smpl <- nrow(merged)
  ncal <- round(n.smpl * cal_prop)
  j <- lapply(seq_len(n_splits), function(k) sort(sample(seq_len(n.smpl), ncal)))

  rm1 <- 1:46; rm2 <- 637:666; rm3 <- 1437:1466; rm4 <- 2127:2151
  ind <- 7:2146
  rmove <- c(1:4, 64:66, 144:146, 213:214)
  vars <- band_cols[seq(11, 2141, 10)][-rmove]

  list(soil = soil, vnir.matrix = vnir.matrix, j = j,
       rm1 = rm1, rm2 = rm2, rm3 = rm3, rm4 = rm4,
       ind = ind, rmove = rmove, vars = vars)
}
