#' Machine Learning Function for Soil Spectral Data
#'
#' This function applies several machine learning models (PCR, PLSR, Random Forest, LASSO, Cubist) 
#' to soil spectral data and compares their performance. Optionally, it can return the best-performing model.
#'
#' @param x A data frame or matrix containing spectral data. Default is `spectral_data`.
#' @param y A vector containing corresponding soil laboratory measurements. Default is `lab_data_in_vector`.
#' @param smoother_selection A parameter specifying the smoothing method to be applied during preprocessing.
#' @param type_of_soil A character string indicating the soil type for model calibration.
#' @param model_selection Logical; if `TRUE` (default), the function returns only the best-performing model.
#'   If `FALSE`, it returns the results from all models.
#'
#' @details
#' The function merges spectral and laboratory data, preprocesses the data, and evaluates the following models:
#' \itemize{
#'   \item PCR (Principal Component Regression)
#'   \item PLSR (Partial Least Squares Regression)
#'   \item RF (Random Forest)
#'   \item LASSO regression
#'   \item Cubist regression
#' }
#' 
#' Each model's performance results are combined into a single results object. If \code{model_selection = TRUE}, 
#' the function returns the model with the highest performance metric (based on the 11th column of the results table).
#'
#' @return A data frame:
#' \describe{
#'   \item{If \code{model_selection = FALSE}}{Returns results for all models.}
#'   \item{If \code{model_selection = TRUE}}{Returns only the best-performing model result.}
#' }
#'
#' @examples
#' \donttest{
#' # Example usage:
#' results <- ml_f(
#'   x,
#'   y,
#'   smoother_selection = "savitzky",
#'   type_of_soil = "loam",
#'   model_selection = TRUE
#' )
#' }
#'
#' @importFrom stats predict
#' @importFrom pls pcr plsr
#' @importFrom randomForest randomForest
#' @importFrom glmnet glmnet
#' @importFrom Cubist cubist
#'
#' @export
ml_f <- function(x, y, smoother_selection, type_of_soil, model_selection = TRUE) {
  xy <- merge_of_lab_and_spectrum(x, y)
  soil <- xy[[1]]
  vnir.matrix <- xy[[2]]
  j <- xy[[3]]
  
  pcr_results <- results(pcr_preprocess(soil, vnir.matrix, j, smoother_selection, type_of_soil))
  plsr_results <- results(plsr_preprocess(soil, vnir.matrix, j, smoother_selection, type_of_soil))
  rf_results <- results(rf_preprocess(soil, vnir.matrix, j, smoother_selection, type_of_soil))
  lasso_results <- results(lasso_preprocess(soil, vnir.matrix, j, smoother_selection, type_of_soil))
  cubist_results <- results(cubist_preprocess(soil, vnir.matrix, j, smoother_selection, type_of_soil))
  
  all_results <- rbind(pcr_results, plsr_results, rf_results, lasso_results, cubist_results)
  
  if (model_selection == FALSE) {
    return(all_results)
  } else {
    return(all_results[which.max(all_results[, 11]), ])
  }
}
