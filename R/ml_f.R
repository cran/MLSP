#' Compare Spectral Preprocessing and Calibration Methods
#'
#' Runs every combination of the selected spectral preprocessing methods and
#' calibration methods on the same calibration/validation splits and ranks
#' the combinations by accuracy.
#'
#' @param x A data frame of soil laboratory data with a \code{LAB_NUM}
#'   identifier column and one column per soil property (passed to
#'   \code{\link{merge_of_lab_and_spectrum}}).
#' @param y A data frame of spectra with a sample identifier column
#'   (\code{Wavelength}, matching \code{LAB_NUM}) and 2151 reflectance columns
#'   (350-2500 nm at 1 nm).
#' @param smoother_selection Preprocessing methods to compare: \code{"all"}
#'   (default; all 26 methods in \code{\link{preprocess_methods}}: raw
#'   reflectance, SNV, \eqn{\log_{10}(1/R)}, and 23 Savitzky-Golay filters), a
#'   character vector of method names, or a single preprocessing function such
#'   as \code{\link{raw}}.
#' @param type_of_soil The soil property to model: a column name or index of
#'   \code{x} (index 1 is the identifier, so the first property is 2).
#' @param model_selection Logical; if \code{TRUE} (default), return the
#'   \code{top_n} best combinations; if \code{FALSE}, return all combinations.
#' @param models Calibration methods to compare: any of \code{"pcr"},
#'   \code{"plsr"}, \code{"lasso"}, \code{"cubist"}, and \code{"rf"}.
#' @param top_n Number of combinations returned when
#'   \code{model_selection = TRUE}.
#' @param rank_by Metric used for ranking, one of the column names of the
#'   output (e.g., \code{"val-rpiq"}, \code{"val-rmse"}, \code{"val-r2"}).
#'   Higher is better for R2 and RPIQ, lower for RMSE and MSE, and closer to
#'   zero for bias.
#' @param verbose Logical; print progress messages.
#' @param ... Further arguments passed to the model functions
#'   (\code{ncomp_max}, \code{rf_ntree}, \code{cv_folds}); see
#'   \code{\link{soil_preprocess}}.
#'
#' @details
#' Spectra are matched to the soil data and split once (see
#' \code{\link{merge_of_lab_and_spectrum}}), so every combination is evaluated
#' on identical calibration and validation sets. Each preprocessing method is
#' computed once and reused by all calibration methods. Metrics are averaged
#' over the splits (see \code{\link{results}}). A combination that fails to fit
#' is reported with a warning and \code{NA} metrics instead of stopping the run.
#' Call \code{\link[base]{set.seed}} first for reproducible splits.
#'
#' Running all 26 distinct preprocessing methods with all five calibration
#' methods fits 130 combinations on each split and may take several minutes
#' or longer for large spectral libraries. Restricting \code{models}, e.g. to
#' \code{c("plsr", "lasso")}, is much faster.
#'
#' The preprocessing methods and the Savitzky-Golay naming
#' (\code{sg.<polynomial>.<window>.<derivative>}) follow Clingensmith et al.
#' (2019); see \code{\link{preprocess_methods}} for how window size,
#' polynomial order, and derivative order affect the spectra.
#'
#' @return A data frame ranked by \code{rank_by}, with columns \code{rank},
#'   \code{preprocessing}, \code{model}, \code{LV}, the cross-validation
#'   metrics (\code{cv-r2}, \code{cv-bias}, \code{cv-rmse}, \code{cv-mse},
#'   \code{cv-rpiq}), and the validation metrics (\code{val-r2},
#'   \code{val-bias}, \code{val-rmse}, \code{val-mse}, \code{val-rpiq}).
#'
#' @references
#' Clingensmith CM, Grunwald S, Wani SP (2019). Evaluation of calibration
#' subsetting and new chemometric methods on the spectral prediction of key
#' soil properties in a data-limited environment. European Journal of Soil
#' Science 70(1). \doi{10.1111/ejss.12753}
#'
#' @seealso \code{\link{preprocess_methods}}, \code{\link{soil_preprocess}}
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
#'
#' # Three preprocessing methods x two calibration methods
#' ml_f(soil_data, spectra,
#'      smoother_selection = c("raw", "snv", "sg.q.11.1"),
#'      type_of_soil = "SOC", models = c("pcr", "plsr"), top_n = 5)
#' \donttest{
#' # All 26 distinct preprocessing methods with PLSR
#' top10 <- ml_f(soil_data, spectra, type_of_soil = "SOC", models = "plsr")
#' top10
#' }
#'
#' @export
ml_f <- function(x, y, smoother_selection = "all", type_of_soil = 2,
                 model_selection = TRUE,
                 models = c("pcr", "plsr", "lasso", "cubist", "rf"),
                 top_n = 10, rank_by = "val-rpiq", verbose = FALSE, ...) {
  models <- match.arg(models, c("pcr", "plsr", "lasso", "cubist", "rf"),
                      several.ok = TRUE)
  metric_names <- c("LV", paste0("cv-", c("r2", "bias", "rmse", "mse", "rpiq")),
                    paste0("val-", c("r2", "bias", "rmse", "mse", "rpiq")))
  if (!rank_by %in% metric_names[-1]) {
    stop(sprintf("'rank_by' must be one of: %s.",
                 paste(metric_names[-1], collapse = ", ")), call. = FALSE)
  }

  # Preprocessing methods
  if (is.function(smoother_selection)) {
    pre_list <- list(smoother_selection)
    names(pre_list) <- deparse(substitute(smoother_selection))[1]
  } else if (is.character(smoother_selection)) {
    pre_names <- if (identical(smoother_selection, "all")) {
      preprocess_methods()$method
    } else {
      unknown <- setdiff(smoother_selection, preprocess_methods()$method)
      equiv <- vapply(unknown, .equivalent_method, character(1))
      if (any(!is.na(equiv))) {
        stop(sprintf(paste0("These settings give exactly the same filter as an ",
                            "included method; use the method shown instead: %s."),
                     paste(sprintf("%s -> %s", unknown[!is.na(equiv)],
                                   equiv[!is.na(equiv)]), collapse = ", ")),
             call. = FALSE)
      }
      if (length(unknown)) {
        stop(sprintf("Unknown preprocessing method(s): %s; see preprocess_methods().",
                     paste(unknown, collapse = ", ")), call. = FALSE)
      }
      unique(smoother_selection)
    }
    pre_list <- as.list(pre_names)
    names(pre_list) <- pre_names
  } else {
    stop("'smoother_selection' must be \"all\", method names, or a function.",
         call. = FALSE)
  }

  xy <- merge_of_lab_and_spectrum(x, y)
  target <- .resolve_target(xy$soil, type_of_soil)
  yv <- xy$soil[[target]]
  if (!is.numeric(yv) || anyNA(yv)) {
    stop(sprintf("Soil property '%s' must be numeric without missing values.", target),
         call. = FALSE)
  }

  rows <- list()
  for (p in names(pre_list)) {
    if (verbose) message("Preprocessing: ", p)
    X <- tryCatch(as.matrix(.apply_preprocess(xy$vnir.matrix, pre_list[[p]])),
                  error = function(e) {
                    warning(sprintf("Preprocessing '%s' failed: %s", p,
                                    conditionMessage(e)), call. = FALSE)
                    NULL
                  })
    for (m in models) {
      if (verbose) message("  model: ", m)
      res <- if (is.null(X)) NULL else tryCatch(
        results(.evaluate_model(m, yv, X, xy$j, ...)),
        error = function(e) {
          warning(sprintf("'%s' with '%s' failed: %s", m, p, conditionMessage(e)),
                  call. = FALSE)
          NULL
        })
      if (is.null(res)) res <- stats::setNames(rep(NA_real_, 11), metric_names)
      rows[[length(rows) + 1]] <- data.frame(
        preprocessing = p, model = m, t(res), check.names = FALSE,
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)

  key <- out[[rank_by]]
  if (grepl("r2$|rpiq$", rank_by)) key <- -key
  if (grepl("bias$", rank_by)) key <- abs(key)
  out <- out[order(key, na.last = TRUE), , drop = FALSE]
  out <- cbind(rank = seq_len(nrow(out)), out)
  rownames(out) <- NULL
  if (model_selection) out <- utils::head(out, top_n)
  out
}
