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
  x <- as.matrix(x)
  times <- ncol(x) / ncols
  if (round(times) != times) stop(paste("Number of columns not divisible by", ncols))
  starts <- seq(1, ncol(x), by = ncols)
  out <- vapply(starts, function(s) rowMeans(x[, s:(s + ncols - 1), drop = FALSE]),
                numeric(nrow(x)))
  matrix(out, nrow = nrow(x))
}
