#' Compute Row-wise Means for Groups of Columns
#'
#' This function computes the row-wise mean for consecutive groups of columns in a matrix or data frame.
#'
#' @param x A numeric matrix or data frame.
#' @param ncols An integer specifying the number of consecutive columns to group together.
#'
#' @return A numeric matrix with the same number of rows as \code{x} and \code{ncol(x) / ncols} columns,
#'   where each column is the row-wise mean of a group of \code{ncols} columns.
#'
#' @details
#' The function splits the columns of \code{x} into consecutive groups of \code{ncols} columns
#' and calculates the mean of each row for each group. The number of columns in \code{x} must
#' be divisible by \code{ncols}.
#'
#' @examples
#' mat <- matrix(1:12, nrow = 3)
#' colsMean(mat, 2)
#'
#' @export
colsMean <- function(x, ncols){
  # Convert input to data frame
  x <- data.frame(x)
  
  # Calculate how many groups of columns we will have
  times <- length(x)/ncols
  
  # Check if ncols divides evenly into number of columns
  if(round(times) != times) stop(paste("Number of columns not divisible by", ncols))
  
  # Initialize an empty matrix to store results
  temp <- matrix(NA, ncol = 0, nrow = nrow(x))
  
  # Sequence of starting indices for each group of columns
  index <- seq(1, length(x), by = ncols)
  
  # Loop through each group and compute row-wise mean
  for(i in 1:times){
    temp <- cbind(temp, apply(x[, index[i]:(index[i] + ncols - 1)], MARGIN = 1, FUN = mean))
  }
  
  # Return result matrix
  return(temp)
}
