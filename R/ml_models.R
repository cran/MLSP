#' Fit and Evaluate One Calibration Method
#'
#' Each function fits one calibration method to preprocessed spectra on every
#' calibration/validation split and computes accuracy metrics with
#' \code{\link{msd.comp}}:
#' \itemize{
#'   \item \code{pcr_preprocess()}: principal component regression (PCR)
#'   \item \code{plsr_preprocess()}: partial least squares regression (PLSR)
#'   \item \code{lasso_preprocess()}: least absolute shrinkage and selection
#'     operator (LASSO)
#'   \item \code{rf_preprocess()}: random forest
#'   \item \code{cubist_preprocess()}: 'Cubist' rule-based regression
#' }
#'
#' @details
#' Within each calibration set, model complexity is tuned and accuracy is
#' estimated by cross-validation (CV): the number of PCR/PLSR components is
#' chosen with the one-standard-error rule of
#' \code{\link[pls]{selectNcomp}} (at most \code{ncomp_max} components), the
#' LASSO penalty is the CV minimum from \code{\link[glmnet]{cv.glmnet}},
#' random forest uses out-of-bag predictions, and 'Cubist' uses k-fold CV
#' predictions. The CV predictions give the \code{cal} metrics; the fitted
#' model is then applied to the held-out validation set (\code{val} metrics).
#' Predictors are scaled for PCR and PLSR.
#'
#' @param soil A data frame of soil properties, e.g. the \code{soil} element
#'   returned by \code{\link{merge_of_lab_and_spectrum}}.
#' @param vnir.matrix A numeric matrix or data frame of spectra (2151 columns,
#'   350-2500 nm).
#' @param j A list of calibration-row index vectors (e.g., from
#'   \code{\link{merge_of_lab_and_spectrum}}).
#' @param preprocess A preprocessing method: a method name from
#'   \code{\link{preprocess_methods}} (e.g., \code{"sg.q.11.1"}) or a function
#'   that takes the spectra and returns a data frame of predictors (e.g.,
#'   \code{\link{raw}}).
#' @param type_of_soil The soil property to model: a column name or index of
#'   \code{soil}.
#' @param ncomp_max Maximum number of PCR/PLSR components considered.
#' @param rf_ntree Number of trees for random forest.
#' @param cv_folds Number of CV folds used within each calibration set.
#'
#' @return A list with one element per split. Each element is a matrix of
#'   the \code{\link{msd.comp}} statistics with columns \code{cal} (CV) and
#'   \code{val} (validation), and attribute \code{LV} (the number of latent
#'   variables for PCR/PLSR, \code{NA} otherwise).
#'
#' @references
#' Mevik B-H, Wehrens R (2007). The pls package: principal component and
#' partial least squares regression in R. Journal of Statistical Software
#' 18(2), 1-23. \doi{10.18637/jss.v018.i02}
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
#' soil_data <- data.frame(LAB_NUM = 1:n, SOC = soc, clay = runif(n, 5, 40))
#' xy <- merge_of_lab_and_spectrum(soil_data, spectra)
#'
#' # PLSR on first-derivative spectra (quadratic, 11-point window)
#' fit_plsr <- plsr_preprocess(xy$soil, xy$vnir.matrix, xy$j,
#'                             preprocess = "sg.q.11.1", type_of_soil = "SOC")
#' results(fit_plsr)
#' \donttest{
#' # Random forest on raw reflectance
#' fit_rf <- rf_preprocess(xy$soil, xy$vnir.matrix, xy$j,
#'                         preprocess = "raw", type_of_soil = "SOC")
#' results(fit_rf)
#' }
#'
#' @seealso \code{\link{ml_f}}, \code{\link{results}}
#'
#' @name soil_preprocess
NULL

# Internal helpers -------------------------------------------------------------

.resolve_target <- function(soil, type_of_soil) {
  if (is.character(type_of_soil)) {
    if (!type_of_soil %in% names(soil)) {
      stop(sprintf("Soil property '%s' not found.", type_of_soil), call. = FALSE)
    }
    return(type_of_soil)
  }
  if (!is.numeric(type_of_soil) || type_of_soil < 1 || type_of_soil > ncol(soil)) {
    stop("'type_of_soil' must be a column name or a valid column index of 'soil'.",
         call. = FALSE)
  }
  names(soil)[type_of_soil]
}

.apply_preprocess <- function(vnir.matrix, preprocess) {
  if (is.character(preprocess)) return(spectral_preprocess(vnir.matrix, preprocess))
  if (is.function(preprocess)) return(as.data.frame(preprocess(vnir.matrix)))
  stop("'preprocess' must be a method name or a function.", call. = FALSE)
}

.cv_folds <- function(n, k) sample(rep_len(seq_len(k), n))

# Fit one model on (Xc, yc); return CV predictions on the calibration set,
# predictions on the validation set, and the number of latent variables.
.fit_predict <- function(model, Xc, yc, Xv, ncomp_max = 20, rf_ntree = 500,
                         cv_folds = 10) {
  n <- nrow(Xc)
  k <- max(3, min(cv_folds, floor(n / 3)))
  if (model %in% c("pcr", "plsr")) {
    nc_max <- max(1, min(ncomp_max, n - ceiling(n / k) - 1, ncol(Xc)))
    dat <- data.frame(y = yc)
    dat$X <- Xc
    fitfun <- if (model == "pcr") pls::pcr else pls::plsr
    mod <- fitfun(y ~ X, data = dat, ncomp = nc_max, scale = TRUE,
                  validation = "CV", segments = k)
    nc <- suppressWarnings(pls::selectNcomp(mod, method = "onesigma"))
    nc <- max(1, nc)
    newd <- data.frame(row = seq_len(nrow(Xv)))
    newd$X <- Xv
    return(list(cal = as.numeric(mod$validation$pred[, 1, nc]),
                val = as.numeric(predict(mod, newdata = newd, ncomp = nc)[, 1, 1]),
                lv = nc))
  }
  if (model == "lasso") {
    cvm <- glmnet::cv.glmnet(Xc, yc, alpha = 1, nfolds = k, keep = TRUE)
    idx <- which(cvm$lambda == cvm$lambda.min)
    return(list(cal = as.numeric(cvm$fit.preval[, idx]),
                val = as.numeric(predict(cvm, newx = Xv, s = "lambda.min")),
                lv = NA_real_))
  }
  if (model == "rf") {
    mod <- randomForest::randomForest(x = Xc, y = yc, ntree = rf_ntree)
    return(list(cal = as.numeric(mod$predicted),
                val = as.numeric(predict(mod, newdata = Xv)),
                lv = NA_real_))
  }
  if (model == "cubist") {
    folds <- .cv_folds(n, k)
    cal <- numeric(n)
    for (f in seq_len(k)) {
      tr <- folds != f
      m_f <- Cubist::cubist(x = as.data.frame(Xc[tr, , drop = FALSE]), y = yc[tr])
      cal[!tr] <- predict(m_f, newdata = as.data.frame(Xc[!tr, , drop = FALSE]))
    }
    mod <- Cubist::cubist(x = as.data.frame(Xc), y = yc)
    return(list(cal = cal,
                val = as.numeric(predict(mod, newdata = as.data.frame(Xv))),
                lv = NA_real_))
  }
  stop(sprintf("Unknown model '%s'.", model), call. = FALSE)
}

# Evaluate one model on precomputed predictors over all splits.
.evaluate_model <- function(model, y, X, j, ...) {
  out <- lapply(seq_along(j), function(s) {
    cal_i <- j[[s]]
    fp <- .fit_predict(model, X[cal_i, , drop = FALSE], y[cal_i],
                       X[-cal_i, , drop = FALSE], ...)
    msd <- cbind(cal = msd.comp(y[cal_i], fp$cal), val = msd.comp(y[-cal_i], fp$val))
    attr(msd, "LV") <- fp$lv
    msd
  })
  names(out) <- paste0("split", seq_along(j))
  out
}

.model_preprocess <- function(model, soil, vnir.matrix, j, preprocess, type_of_soil, ...) {
  target <- .resolve_target(soil, type_of_soil)
  X <- as.matrix(.apply_preprocess(vnir.matrix, preprocess))
  .evaluate_model(model, soil[[target]], X, j, ...)
}

#' @rdname soil_preprocess
#' @export
pcr_preprocess <- function(soil, vnir.matrix, j, preprocess, type_of_soil,
                           ncomp_max = 20, cv_folds = 10) {
  .model_preprocess("pcr", soil, vnir.matrix, j, preprocess, type_of_soil,
                    ncomp_max = ncomp_max, cv_folds = cv_folds)
}

#' @rdname soil_preprocess
#' @export
plsr_preprocess <- function(soil, vnir.matrix, j, preprocess, type_of_soil,
                            ncomp_max = 20, cv_folds = 10) {
  .model_preprocess("plsr", soil, vnir.matrix, j, preprocess, type_of_soil,
                    ncomp_max = ncomp_max, cv_folds = cv_folds)
}

#' @rdname soil_preprocess
#' @export
lasso_preprocess <- function(soil, vnir.matrix, j, preprocess, type_of_soil,
                             cv_folds = 10) {
  .model_preprocess("lasso", soil, vnir.matrix, j, preprocess, type_of_soil,
                    cv_folds = cv_folds)
}

#' @rdname soil_preprocess
#' @export
rf_preprocess <- function(soil, vnir.matrix, j, preprocess, type_of_soil,
                          rf_ntree = 500) {
  .model_preprocess("rf", soil, vnir.matrix, j, preprocess, type_of_soil,
                    rf_ntree = rf_ntree)
}

#' @rdname soil_preprocess
#' @export
cubist_preprocess <- function(soil, vnir.matrix, j, preprocess, type_of_soil,
                              cv_folds = 10) {
  .model_preprocess("cubist", soil, vnir.matrix, j, preprocess, type_of_soil,
                    cv_folds = cv_folds)
}
