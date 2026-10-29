#' Summarize Model Performance Metrics
#'
#' Averages accuracy metrics over the calibration/validation splits returned
#' by \code{\link{pcr_preprocess}}, \code{\link{plsr_preprocess}},
#' \code{\link{lasso_preprocess}}, \code{\link{rf_preprocess}}, or
#' \code{\link{cubist_preprocess}}.
#'
#' @param metric.list A list of \code{\link{msd.comp}} matrices (one per
#'   split) with columns \code{cal} and \code{val}.
#' @param soil_type Not used; kept for backward compatibility.
#'
#' @return A named numeric vector of metrics averaged over splits:
#'   \describe{
#'     \item{LV}{Number of latent variables (PCR/PLSR; \code{NA} otherwise).}
#'     \item{cv-r2, cv-bias, cv-rmse, cv-mse, cv-rpiq}{Coefficient of
#'       determination, bias, root mean squared error, mean squared error, and
#'       ratio of performance to interquartile distance from cross-validation
#'       within the calibration set.}
#'     \item{val-r2, val-bias, val-rmse, val-mse, val-rpiq}{The same metrics
#'       on the independent validation set.}
#'   }
#'
#' @references
#' Bellon-Maurel V, Fernandez-Ahumada E, Palagos B, Roger J-M, McBratney A
#' (2010). Critical review of chemometric indicators commonly used for
#' assessing the quality of the prediction of soil attributes by NIR
#' spectroscopy. TrAC Trends in Analytical Chemistry 29, 1073-1081.
#' \doi{10.1016/j.trac.2010.05.006}
#'
#' @seealso \code{\link{soil_preprocess}}, \code{\link{ml_f}}
#'
#' @examples
#' # Two splits of observed and predicted values
#' split_metrics <- function(y, yc, yv) cbind(cal = msd.comp(y, yc), val = msd.comp(y, yv))
#' y <- c(1.2, 3.4, 2.5, 4.1, 2.2)
#' metric.list <- list(
#'   split1 = split_metrics(y, y + 0.1, y - 0.2),
#'   split2 = split_metrics(y, y - 0.1, y + 0.3))
#' results(metric.list)
#'
#' @export
results <- function(metric.list, soil_type) {
  stats_rows <- c("r2", "bias", "rmse", "mse", "rpiq")
  tab <- t(vapply(metric.list, function(m) {
    lv <- attr(m, "LV")
    c(LV = if (is.null(lv)) NA_real_ else as.numeric(lv),
      stats::setNames(m[stats_rows, "cal"], paste0("cv-", stats_rows)),
      stats::setNames(m[stats_rows, "val"], paste0("val-", stats_rows)))
  }, numeric(11)))
  colMeans(tab)
}
