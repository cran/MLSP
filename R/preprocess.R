## Spectral preprocessing for 350-2500 nm VNIR spectra recorded at 1 nm
## (2151 bands). Every method returns 202 band-averaged (10 nm) predictors
## after removing the detector-artifact regions near 1000 nm and 1800 nm and
## the noisy spectral edges.

# Band layout shared by all preprocessing methods ----------------------------
.n_bands <- 2151L
.agg_index <- 7:2146                                    # columns averaged in 10-band groups
.agg_remove <- c(1:4, 64:66, 144:146, 213:214)          # 10-nm groups dropped after averaging
.artifact_bands <- c(1:46, 637:666, 1437:1466, 2127:2151) # 1-nm bands excluded (edges, detector joins)

.check_spectra <- function(vnir.matrix) {
  if (!(is.data.frame(vnir.matrix) || is.matrix(vnir.matrix))) {
    stop("'vnir.matrix' must be a data frame or matrix of spectra (samples in rows).",
         call. = FALSE)
  }
  if (ncol(vnir.matrix) != .n_bands) {
    stop(sprintf(paste0("'vnir.matrix' must have %d columns (350-2500 nm at 1 nm); ",
                        "it has %d."), .n_bands, ncol(vnir.matrix)), call. = FALSE)
  }
  X <- as.matrix(vnir.matrix)
  if (!is.numeric(X)) stop("'vnir.matrix' must contain only numeric values.", call. = FALSE)
  if (anyNA(X)) stop("'vnir.matrix' contains missing values.", call. = FALSE)
  band_names <- colnames(vnir.matrix)
  if (is.null(band_names)) band_names <- as.character(350:2500)
  dimnames(X) <- NULL
  list(X = X, band_names = band_names)
}

# Average 1-nm bands into 10-nm bands and drop artifact regions.
.aggregate_bands <- function(X, band_names) {
  out <- data.frame(colsMean(X[, .agg_index, drop = FALSE], 10))
  out <- out[, -.agg_remove, drop = FALSE]
  names(out) <- band_names[seq(11, 2141, 10)][-.agg_remove]
  out
}

# Savitzky-Golay filter applied to every spectrum (row) of X.
# The matrix is passed unnamed: gsignal::sgolayfilt() fails on named vectors.
.sg_filter <- function(X, p, n, m) {
  t(gsignal::sgolayfilt(t(X), p = p, n = n, m = m))
}

#' Catalogue of Spectral Preprocessing Methods
#'
#' Lists the 26 spectral preprocessing methods in the package: raw
#' reflectance, standard normal variate (SNV), pseudo-absorbance
#' (\eqn{\log_{10}(1/R)}), and 23 Savitzky-Golay (SG) filters.
#'
#' @section Savitzky-Golay filters:
#' An SG filter fits a polynomial by least squares to the points inside a
#' moving window centred on each wavelength and replaces the centre value with
#' the fitted value (smoothing) or with a derivative of the fitted polynomial
#' (Savitzky and Golay, 1964). Three settings define each filter, following the
#' \code{SG.P.N.D} naming used by Clingensmith et al. (2019), who applied SG
#' smoothers and first and second derivatives with 11-point windows to soil
#' VNIR spectra:
#' \describe{
#'   \item{Polynomial order (P)}{\code{l} = linear (1), \code{q} = quadratic
#'     (2), \code{c} = cubic (3). Higher orders follow narrow absorption
#'     features more closely; lower orders smooth more strongly. A linear SG
#'     smoother is a simple moving average.}
#'   \item{Window size (N)}{3, 5, 7, 9, or 11 points (nm). Wider windows
#'     remove more noise but can flatten narrow absorption features. The
#'     window must exceed the polynomial order.}
#'   \item{Derivative (D)}{0 = smoothing, 1 = first derivative, which removes
#'     additive baseline offsets, and 2 = second derivative, which also removes
#'     linear baseline slopes and separates overlapping absorption bands.}
#' }
#' The 23 filters are linear smoothers (3- to 11-point windows), quadratic
#' smoothers (5 to 11), quadratic first derivatives (3 to 11), cubic first
#' derivatives (5 to 11), and quadratic second derivatives (3 to 11). For
#' example, \code{sg.q.11.1} is the first derivative from a quadratic
#' polynomial fitted in an 11-point window.
#'
#' @section Equivalent settings that are not included:
#' Because the window is symmetric, raising the polynomial order from even to
#' odd does not change the fitted value, slope, or curvature at the window
#' centre (Savitzky and Golay, 1964). Fourteen further settings therefore
#' reproduce a method above exactly and are not included:
#' \itemize{
#'   \item cubic smoothers equal quadratic smoothers (\code{sg.c.N.0} =
#'     \code{sg.q.N.0}, N = 5 to 11);
#'   \item a 3-point quadratic smoother returns the spectrum unchanged
#'     (\code{sg.q.3.0} = \code{raw});
#'   \item linear first derivatives equal quadratic first derivatives
#'     (\code{sg.l.N.1} = \code{sg.q.N.1}, N = 3 to 11);
#'   \item cubic second derivatives equal quadratic second derivatives
#'     (\code{sg.c.N.2} = \code{sg.q.N.2}, N = 5 to 11).
#' }
#' Results reported under these names, such as \code{SG.C.11.0} in
#' Clingensmith et al. (2019), are reproduced by the equivalent method.
#'
#' @return A data frame with one row per method and columns
#'   \code{method} (the name used in \code{\link{spectral_preprocess}} and
#'   \code{\link{ml_f}}), \code{family}, \code{poly_order}, \code{window},
#'   \code{derivative}, \code{description}, and \code{equivalent_to} (other
#'   SG settings that give the same filter).
#'
#' @references
#' Savitzky A, Golay MJE (1964). Smoothing and differentiation of data by
#' simplified least squares procedures. Analytical Chemistry 36, 1627-1639.
#' \doi{10.1021/ac60214a047}
#'
#' Barnes RJ, Dhanoa MS, Lister SJ (1989). Standard normal variate
#' transformation and de-trending of near-infrared diffuse reflectance
#' spectra. Applied Spectroscopy 43, 772-777. \doi{10.1366/0003702894202201}
#'
#' Rinnan A, van den Berg F, Engelsen SB (2009). Review of the most common
#' pre-processing techniques for near-infrared spectra. TrAC Trends in
#' Analytical Chemistry 28, 1201-1222. \doi{10.1016/j.trac.2009.07.007}
#'
#' Clingensmith CM, Grunwald S, Wani SP (2019). Evaluation of calibration
#' subsetting and new chemometric methods on the spectral prediction of key
#' soil properties in a data-limited environment. European Journal of Soil
#' Science 70(1). \doi{10.1111/ejss.12753}
#'
#' @seealso \code{\link{spectral_preprocess}}, \code{\link{ml_f}}
#'
#' @examples
#' methods <- preprocess_methods()
#' table(methods$family)
#' # Savitzky-Golay filters by polynomial order and derivative
#' sg <- methods[methods$family == "Savitzky-Golay", ]
#' table(sg$poly_order, sg$derivative)
#' methods[, c("method", "description", "equivalent_to")]
#'
#' @export
preprocess_methods <- function() {
  poly_label <- c("1" = "linear", "2" = "quadratic", "3" = "cubic")
  deriv_label <- c("0" = "smoothing", "1" = "1st derivative", "2" = "2nd derivative")
  sg <- rbind(
    expand.grid(window = c(3, 5, 7, 9, 11), poly_order = 1, derivative = 0),
    expand.grid(window = c(5, 7, 9, 11), poly_order = 2, derivative = 0),
    expand.grid(window = c(3, 5, 7, 9, 11), poly_order = 2, derivative = 1),
    expand.grid(window = c(5, 7, 9, 11), poly_order = 3, derivative = 1),
    expand.grid(window = c(3, 5, 7, 9, 11), poly_order = 2, derivative = 2)
  )
  sg_method <- sprintf("sg.%s.%d.%d", c("l", "q", "c")[sg$poly_order], sg$window, sg$derivative)
  equiv <- rep(NA_character_, nrow(sg))
  q_smooth <- sg$poly_order == 2 & sg$derivative == 0
  equiv[q_smooth] <- sprintf("sg.c.%d.0", sg$window[q_smooth])
  q_d1 <- sg$poly_order == 2 & sg$derivative == 1
  equiv[q_d1] <- sprintf("sg.l.%d.1", sg$window[q_d1])
  q_d2 <- sg$poly_order == 2 & sg$derivative == 2 & sg$window >= 5
  equiv[q_d2] <- sprintf("sg.c.%d.2", sg$window[q_d2])

  out <- data.frame(
    method = c("raw", "snv", "log1_R", sg_method),
    family = c("none", "scatter correction", "transformation",
               rep("Savitzky-Golay", nrow(sg))),
    poly_order = c(NA, NA, NA, sg$poly_order),
    window = c(NA, NA, NA, sg$window),
    derivative = c(NA, NA, NA, sg$derivative),
    description = c(
      "Raw reflectance",
      "Standard normal variate",
      "Pseudo-absorbance, log10(1/R)",
      sprintf("SG %s %s, %d-point window", poly_label[as.character(sg$poly_order)],
              deriv_label[as.character(sg$derivative)], sg$window)
    ),
    equivalent_to = c("sg.q.3.0", NA, NA, equiv),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

# Map an excluded, equivalent SG name to the method that reproduces it.
.equivalent_method <- function(method) {
  m <- regmatches(method, regexec("^sg\\.([lqc])\\.(\\d+)\\.([012])$", method))[[1]]
  if (length(m) == 0) return(NA_character_)
  p <- m[2]; n <- as.integer(m[3]); d <- m[4]
  if (p == "q" && n == 3 && d == "0") return("raw")
  if (p == "c" && d == "0" && n >= 5) return(sprintf("sg.q.%d.0", n))
  if (p == "l" && d == "1") return(sprintf("sg.q.%d.1", n))
  if (p == "c" && d == "2" && n >= 5) return(sprintf("sg.q.%d.2", n))
  NA_character_
}

#' Preprocess VNIR Spectra
#'
#' Applies one of the 26 spectral preprocessing methods listed by
#' \code{\link{preprocess_methods}} to VNIR spectra measured from 350 to
#' 2500 nm at 1-nm intervals (2151 bands), averages the result into 10-nm
#' bands, and removes the noisy spectral edges and the detector-artifact
#' regions near 1000 nm and 1800 nm (350-395, 986-1015, 1786-1815, and
#' 2476-2500 nm), leaving 202 predictor bands.
#'
#' @param vnir.matrix A numeric data frame or matrix with one reflectance
#'   spectrum per row and 2151 columns (350-2500 nm).
#' @param method Name of the preprocessing method: \code{"raw"},
#'   \code{"snv"}, \code{"log1_R"}, or a Savitzky-Golay filter named
#'   \code{sg.<polynomial>.<window>.<derivative>} (see
#'   \code{\link{preprocess_methods}}).
#'
#' @details
#' Filters and transformations are applied to the full 1-nm spectrum before
#' band averaging, so Savitzky-Golay windows are in nanometres and derivatives
#' are per nanometre. SNV centres and scales each spectrum with the mean and
#' standard deviation of the retained (non-artifact) bands (Barnes et al.,
#' 1989). The removed regions follow the detector edge effects described by
#' Clingensmith et al. (2019).
#'
#' @return A data frame with one row per sample and 202 columns named by
#'   the wavelength at the centre of each 10-nm band.
#'
#' @references
#' Barnes RJ, Dhanoa MS, Lister SJ (1989). Standard normal variate
#' transformation and de-trending of near-infrared diffuse reflectance
#' spectra. Applied Spectroscopy 43, 772-777. \doi{10.1366/0003702894202201}
#'
#' Clingensmith CM, Grunwald S, Wani SP (2019). Evaluation of calibration
#' subsetting and new chemometric methods on the spectral prediction of key
#' soil properties in a data-limited environment. European Journal of Soil
#' Science 70(1). \doi{10.1111/ejss.12753}
#'
#' @seealso \code{\link{preprocess_methods}}, \code{\link{raw}}
#'
#' @examples
#' wl <- 350:2500
#' vnir <- data.frame(t(sapply(c(0.30, 0.35, 0.40), function(b)
#'   b + 0.1 * (wl - 350) / 2150 - 0.03 * exp(-((wl - 1900) / 40)^2))))
#' names(vnir) <- wl
#' smooth <- spectral_preprocess(vnir, "sg.l.11.0")  # moving average, 11 nm
#' d1 <- spectral_preprocess(vnir, "sg.q.11.1")      # 1st derivative, quadratic
#' d2 <- spectral_preprocess(vnir, "sg.q.11.2")      # 2nd derivative, quadratic
#' dim(d1)
#' # The first derivative removes the offsets between the three spectra
#' round(apply(d1, 2, sd)[1:5], 6)
#'
#' @export
spectral_preprocess <- function(vnir.matrix, method = "raw") {
  if (!is.character(method) || length(method) != 1L) {
    stop("'method' must be a single method name; see preprocess_methods().", call. = FALSE)
  }
  catalogue <- preprocess_methods()
  if (!method %in% catalogue$method) {
    eq <- .equivalent_method(method)
    if (!is.na(eq)) {
      stop(sprintf("'%s' gives exactly the same filter as '%s'; use '%s'.",
                   method, eq, eq), call. = FALSE)
    }
    stop(sprintf("Unknown preprocessing method '%s'; see preprocess_methods().", method),
         call. = FALSE)
  }

  chk <- .check_spectra(vnir.matrix)
  X <- chk$X
  info <- catalogue[catalogue$method == method, ]
  Y <- switch(
    info$family,
    "none" = X,
    "transformation" = {
      if (any(X <= 0)) stop("log10(1/R) requires strictly positive reflectance.", call. = FALSE)
      log10(1 / X)
    },
    "scatter correction" = {
      keep <- setdiff(seq_len(.n_bands), .artifact_bands)
      mu <- rowMeans(X[, keep, drop = FALSE])
      s <- apply(X[, keep, drop = FALSE], 1, stats::sd)
      (X - mu) / s
    },
    "Savitzky-Golay" = .sg_filter(X, p = info$poly_order, n = info$window,
                                  m = info$derivative)
  )
  .aggregate_bands(Y, chk$band_names)
}

#' Raw Reflectance Predictors
#'
#' Averages raw VNIR reflectance (350-2500 nm, 1-nm intervals) into 10-nm
#' bands and removes the noisy spectral edges and the detector-artifact regions
#' near 1000 nm and 1800 nm. Equivalent to
#' \code{spectral_preprocess(vnir.matrix, "raw")}.
#'
#' @inheritParams spectral_preprocess
#'
#' @return A data frame with one row per sample and 202 columns.
#'
#' @seealso \code{\link{spectral_preprocess}}
#'
#' @examples
#' vnir <- data.frame(matrix(runif(3 * 2151, 0.2, 0.6), nrow = 3))
#' names(vnir) <- 350:2500
#' dim(raw(vnir))
#'
#' @export
raw <- function(vnir.matrix) {
  spectral_preprocess(vnir.matrix, "raw")
}
