#' Compute Model Evaluation Metrics
#'
#' This function computes various statistics for comparing observed values `y` with predicted values `yhat`.
#' It includes correlation, regression coefficients, bias, RMSE, MSE, and predictive performance metrics like RPD and RPIQ.
#'
#' @param y Numeric vector of observed values.
#' @param yhat Numeric vector of predicted values (same length as `y`).
#'
#' @return A named numeric vector with the following components:
#' \describe{
#'   \item{r}{Pearson correlation between `y` and `yhat`}
#'   \item{int}{Intercept of regression of `y` on `yhat`}
#'   \item{slope}{Slope of regression of `y` on `yhat`}
#'   \item{r2}{Coefficient of determination (R-squared)}
#'   \item{bias}{Mean bias: mean(yhat) - mean(y)}
#'   \item{rmse}{Root mean squared error}
#'   \item{mse}{Mean squared error}
#'   \item{sb}{Systematic bias component of MSE}
#'   \item{nu}{Non-unity slope component of MSE}
#'   \item{lc}{Lack-of-correlation component of MSE}
#'   \item{rmse.c}{Corrected RMSE after removing bias}
#'   \item{mse.c}{Corrected MSE after removing bias}
#'   \item{rpd}{Ratio of standard deviation to RMSE (RPD)}
#'   \item{rpiq}{Ratio of interquartile range to RMSE (RPIQ)}
#' }
#'
#' @examples
#' \donttest{
#' y_obs <- c(1.2, 3.4, 2.5, 4.1)
#' y_pred <- c(1.1, 3.5, 2.4, 4.0)
#' msd.comp(y_obs, y_pred)
#' }
#'
#' @export
msd.comp <- function(y, yhat){
  n <- length(y)
  r <- cor(y, yhat)
  
  lmy <- lm(y ~ yhat)
  a <- coefficients(lmy)[1]
  b <- coefficients(lmy)[2]
  r2 <- summary(lmy)$r.squared[1]
  
  bias <- mean(yhat) - mean(y)
  mse <- sum((yhat - y)^2) / n
  rmse <- sqrt(mse)
  mse.c <- sum((yhat - bias - y)^2) / n
  rmse.c <- sqrt(mse.c)
  
  sb <- (mean(yhat) - mean(y))^2
  nu <- ((1 - b)^2) * (var(yhat) * ((n - 1) / n))
  lc <- (1 - r^2) * (var(y) * ((n - 1) / n))
  
  rpd <- sd(y) / rmse
  q1 <- quantile(y)[2]
  q3 <- quantile(y)[4]
  rpiq <- (q3 - q1) / rmse
  
  msd.vec <- round(c(r, a, b, r2, bias, rmse, mse, sb, nu, lc, rmse.c, mse.c, rpd, rpiq), 5)
  names(msd.vec) <- c("r", "int", "slope", "r2", "bias", "rmse", "mse", "sb", "nu", "lc", 
                      "rmse.c", "mse.c", "rpd", "rpiq")
  return(msd.vec)
}
