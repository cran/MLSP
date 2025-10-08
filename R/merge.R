#' Merge Soil Laboratory Data with Spectral Data
#'
#' This function merges soil laboratory data with cleaned spectral (VNIR) data, 
#' performs preprocessing, and prepares inputs for calibration and model building.  
#'
#' @param soil_data A data frame containing soil laboratory measurements (must include a column named \code{LAB_NUM}).
#' @param data_NaturaSpec_cleaned A data frame containing cleaned spectral data 
#'   with columns \code{Wavelength}, \code{LAB_NUM}, and reflectance values.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Aggregates spectral data by wavelength and computes mean reflectance values.
#'   \item Merges the soil and spectral datasets by \code{LAB_NUM}.
#'   \item Separates soil variables and VNIR spectral matrix.
#'   \item Creates calibration sample indices using random sampling.
#'   \item Defines spectral bands to remove (detector artifact areas) and 
#'         indices to be used in modeling.
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{soil}{Data frame of soil laboratory data (first 8 columns of merged dataset).}
#'   \item{vnir.matrix}{Matrix of VNIR spectral reflectance values (without metadata columns).}
#'   \item{j}{List of calibration sample indices for cross-validation (4 sets).}
#'   \item{rm1, rm2, rm3, rm4}{Vectors of indices corresponding to spectral bands to be removed 
#'         (detector artifact regions around 1000 nm and 1800 nm).}
#'   \item{ind}{Indices of spectral bands used for aggregation (columns 7–2146).}
#'   \item{rmove}{Indices of bands to be excluded from analysis.}
#'   \item{vars}{Vector of spectral band names retained after removal.}
#' }
#'
#' @examples
#' \donttest{
#' merged <- merge_of_lab_and_spectrum(soil_data, data_NaturaSpec_cleaned)
#' str(merged)
#' }
#'
#' @importFrom stats aggregate merge
#'
#' @export
merge_of_lab_and_spectrum <- function(soil_data, data_NaturaSpec_cleaned) {
  
  NS_Spectral <- aggregate(
    x = data_NaturaSpec_cleaned,
    by = list(data_NaturaSpec_cleaned$Wavelength),
    FUN = mean
  )
  
  colnames(NS_Spectral)[1] <- "LAB_NUM"
  
  Soil_NS_Spectral <- merge(soil_data, NS_Spectral, by = "LAB_NUM")
  Soil_NS_Spectral_master <- Soil_NS_Spectral[, -9]
  soil <- Soil_NS_Spectral_master[, c(1:8)]
  
  soil.var <- names(soil) # include lab_num
  
  # VNIR data source file
  vnir <- Soil_NS_Spectral_master[, -c(2:8)]
  vnir.matrix <- vnir[, -1]
  
  ## sample and variables
  n.smpl <- nrow(vnir)         # number of samples
  ncal <- round(n.smpl * .7, 0) # calibration set size
  
  j <- list()
  for (jj in 1:4) {
    j[[jj]] <- sample(1:n.smpl, ncal)
  }
  
  rm1 <- c(seq.int(1, 46, 1))
  rm2 <- c(seq.int(637, 666, 1))
  rm3 <- c(seq.int(1437, 1466, 1))
  rm4 <- c(seq.int(2127, 2151, 1)) # detector artifact areas
  
  ind <- 7:2146   # indices used in aggregation
  rmove <- c(1:4, 64:66, 144:146, 213:214) # remove these bands
  
  vars <- names(vnir.matrix[, seq(11, 2141, 10)])
  vars <- vars[-rmove]  # retained band names
  
  return(list(soil, vnir.matrix, j, rm1, rm2, rm3, rm4, ind, rmove, vars))
}
