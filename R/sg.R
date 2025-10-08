#' Aggregate VNIR Spectra by Columns
#'
#' This function aggregates VNIR (Visible and Near-Infrared) spectral data by calculating the mean of every 10 columns 
#' while removing specific detector artifact regions (~1000nm and ~1800nm) and unwanted spectral bands.
#'
#' @param vnir.matrix A numeric matrix or data frame containing VNIR spectral data. Each row corresponds to a sample, 
#' and each column corresponds to a spectral band.
#'
#' @return A data frame containing the aggregated VNIR spectra with cleaned band names. 
#' The number of columns is reduced by averaging over every 10 bands and removing artifact-prone regions.
#'
#' @details 
#' The function removes columns corresponding to detector artifacts:
#' - rm1: bands 1–46
#' - rm2: bands 637–666
#' - rm3: bands 1437–1466
#' - rm4: bands 2127–2151
#' Additionally, columns 1:4, 64:66, 144:146, and 213:214 (after averaging) are removed.
#'
#' @examples
#' \donttest{
#'   raw_spectra <- raw(vnir_matrix)
#' }
#' 
#' @export
raw<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
raw <- data.frame(colsMean(vnir.matrix[,ind], 10)); raw <- raw[,-rmove]; names(raw) <- vars
return(raw)
}
snv<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
snv <- t(apply(vnir.matrix[-c(rm1,rm2,rm3,rm4)], MARGIN=1, FUN=function(x){(x-mean(x))/sd(x)}))
snv <- data.frame(colsMean(snv, 10)); names(snv) <- vars
return(snv)
}
# log(1/R) "Absorbance"
log1_R<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
log1_R <- t(apply(vnir.matrix, MARGIN=1, FUN=function(x){log10(1/x)}))
log1_R <- data.frame(colsMean(log1_R[,ind], 10)); log1_R <- log1_R[,-rmove]; names(log1_R) <- vars
return(log1_R)
}
sg.l.3.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.3.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=3))    #S-G linear, 3 point moving window
sg.l.3.0 <- data.frame(colsMean(sg.l.3.0[,ind], 10)); sg.l.3.0 <- sg.l.3.0[,-rmove]; names(sg.l.3.0) <- vars
return(sg.l.3.0)
}
sg.l.5.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.5.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=5))    #S-G linear, 5 point moving window
sg.l.5.0 <- data.frame(colsMean(sg.l.5.0[,ind], 10)); sg.l.5.0 <- sg.l.5.0[,-rmove]; names(sg.l.5.0) <- vars
return(sg.l.5.0)
}
sg.l.7.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.7.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=7))    #S-G linear, 7 point moving window
sg.l.7.0 <- data.frame(colsMean(sg.l.7.0[,ind], 10)); sg.l.7.0 <- sg.l.7.0[,-rmove]; names(sg.l.7.0) <- vars
return(sg.l.7.0)
}
sg.l.9.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.9.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=9))    #S-G linear, 9 point moving window
sg.l.9.0 <- data.frame(colsMean(sg.l.9.0[,ind], 10)); sg.l.9.0 <- sg.l.9.0[,-rmove]; names(sg.l.9.0) <- vars
return(sg.l.9.0)
}
sg.l.11.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.11.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=11))    #S-G linear, 11 point moving window
sg.l.11.0 <- data.frame(colsMean(sg.l.11.0[,ind], 10)); sg.l.11.0 <- sg.l.11.0[,-rmove]; names(sg.l.11.0) <- vars
return(sg.l.11.0)
}

# SG Quadratic smoother
sg.q.3.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.3.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=3))    #S-G quadratic, 3 point moving window
sg.q.3.0 <- data.frame(colsMean(sg.q.3.0[,ind], 10)); sg.q.3.0 <- sg.q.3.0[,-rmove]; names(sg.q.3.0) <- vars
return(sg.q.3.0)
}
sg.q.5.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.5.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=5))    #S-G quadratic, 5 point moving window
sg.q.5.0 <- data.frame(colsMean(sg.q.5.0[,ind], 10)); sg.q.5.0 <- sg.q.5.0[,-rmove]; names(sg.q.5.0) <- vars
return(sg.q.5.0)
}
sg.q.7.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.7.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=7))    #S-G quadratic, 7 point moving window
sg.q.7.0 <- data.frame(colsMean(sg.q.7.0[,ind], 10)); sg.q.7.0 <- sg.q.7.0[,-rmove]; names(sg.q.7.0) <- vars
return(sg.q.7.0)
}
sg.q.9.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.9.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=9))    #S-G quadratic, 9 point moving window
sg.q.9.0 <- data.frame(colsMean(sg.q.9.0[,ind], 10)); sg.q.9.0 <- sg.q.9.0[,-rmove]; names(sg.q.9.0) <- vars
return(sg.q.7.0)
}
sg.q.11.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.11.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=11))    #S-G quadratic, 11 point moving window
sg.q.11.0 <- data.frame(colsMean(sg.q.11.0[,ind], 10)); sg.q.11.0 <- sg.q.11.0[,-rmove]; names(sg.q.11.0) <- vars
return(sg.q.7.0)
}

# SG Cubic smoother
sg.c.5.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.5.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=5))    #S-G cubic, 5 point moving window
sg.c.5.0 <- data.frame(colsMean(sg.c.5.0[,ind], 10)); sg.c.5.0 <- sg.c.5.0[,-rmove]; names(sg.c.5.0) <- vars
return(sg.c.5.0)
}
sg.c.7.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.7.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=7))    #S-G cubic, 7 point moving window
sg.c.7.0 <- data.frame(colsMean(sg.c.7.0[,ind], 10)); sg.c.7.0 <- sg.c.7.0[,-rmove]; names(sg.c.7.0) <- vars
return(sg.c.7.0)
}
sg.c.9.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.9.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=9))    #S-G cubic, 9 point moving window
sg.c.9.0 <- data.frame(colsMean(sg.c.9.0[,ind], 10)); sg.c.9.0 <- sg.c.9.0[,-rmove]; names(sg.c.9.0) <- vars
return(sg.c.9.0)
}
sg.c.11.0<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.11.0 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=11))    #S-G cubic, 11 point moving window
sg.c.11.0 <- data.frame(colsMean(sg.c.11.0[,ind], 10)); sg.c.11.0 <- sg.c.11.0[,-rmove]; names(sg.c.11.0) <- vars
return(sg.c.11.0)
}


# SG Linear 1st Dervative
sg.l.3.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.3.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=3, m=1))    #S-G linear, 3 point moving window
sg.l.3.1 <- data.frame(colsMean(sg.l.3.1[,ind], 10)); sg.l.3.1 <- sg.l.3.1[,-rmove]; names(sg.l.3.1) <- vars
return(sg.l.3.1)
}
sg.l.5.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.5.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=5, m=1))    #S-G linear, 5 point moving window
sg.l.5.1 <- data.frame(colsMean(sg.l.5.1[,ind], 10)); sg.l.5.1 <- sg.l.5.1[,-rmove]; names(sg.l.5.1) <- vars
return(sg.l.5.1)
}
sg.l.7.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.7.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=7, m=1))    #S-G linear, 7 point moving window
sg.l.7.1 <- data.frame(colsMean(sg.l.7.1[,ind], 10)); sg.l.7.1 <- sg.l.7.1[,-rmove]; names(sg.l.7.1) <- vars
return(sg.l.7.1)
}
sg.l.9.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.9.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=9, m=1))    #S-G linear, 9 point moving window
sg.l.9.1 <- data.frame(colsMean(sg.l.9.1[,ind], 10)); sg.l.9.1 <- sg.l.9.1[,-rmove]; names(sg.l.9.1) <- vars
return(sg.l.9.1)
}
sg.l.11.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.l.11.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=1, n=11, m=1))    #S-G linear, 11 point moving window
sg.l.11.1 <- data.frame(colsMean(sg.l.11.1[,ind], 10)); sg.l.11.1 <- sg.l.11.1[,-rmove]; names(sg.l.11.1) <- vars
return(sg.l.11.1)
}

# SG Quadratic 1st Derivative
sg.q.3.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.3.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=3, m=1))    #S-G quadratic, 3 point moving window
sg.q.3.1 <- data.frame(colsMean(sg.q.3.1[,ind], 10)); sg.q.3.1 <- sg.q.3.1[,-rmove]; names(sg.q.3.1) <- vars
return(sg.q.3.1)
}
sg.q.5.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.5.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=5, m=1))    #S-G quadratic, 5 point moving window
sg.q.5.1 <- data.frame(colsMean(sg.q.5.1[,ind], 10)); sg.q.5.1 <- sg.q.5.1[,-rmove]; names(sg.q.5.1) <- vars
return(sg.q.5.1)
}
sg.q.7.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.7.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=7, m=1))    #S-G quadratic, 7 point moving window
sg.q.7.1 <- data.frame(colsMean(sg.q.7.1[,ind], 10)); sg.q.7.1 <- sg.q.7.1[,-rmove]; names(sg.q.7.1) <- vars
return(sg.q.7.1)
}
sg.q.9.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.9.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=9, m=1))    #S-G quadratic, 9 point moving window
sg.q.9.1 <- data.frame(colsMean(sg.q.9.1[,ind], 10)); sg.q.9.1 <- sg.q.9.1[,-rmove]; names(sg.q.9.1) <- vars
return(sg.q.9.1)
}
sg.q.11.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.11.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=11, m=1))    #S-G quadratic, 11 point moving window
sg.q.11.1 <- data.frame(colsMean(sg.q.11.1[,ind], 10)); sg.q.11.1 <- sg.q.11.1[,-rmove]; names(sg.q.11.1) <- vars
return(sg.q.11.1)
}

# SG Cubic 1st Derivative
sg.c.5.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.5.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=5, m=1))    #S-G cubic, 5 point moving window
sg.c.5.1 <- data.frame(colsMean(sg.c.5.1[,ind], 10)); sg.c.5.1 <- sg.c.5.1[,-rmove]; names(sg.c.5.1) <- vars
return(sg.c.5.1)
}
sg.c.7.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.7.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=7, m=1))    #S-G cubic, 7 point moving window
sg.c.7.1 <- data.frame(colsMean(sg.c.7.1[,ind], 10)); sg.c.7.1 <- sg.c.7.1[,-rmove]; names(sg.c.7.1) <- vars
return(sg.c.7.1)
}
sg.c.9.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.9.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=9, m=1))    #S-G cubic, 9 point moving window
sg.c.9.1 <- data.frame(colsMean(sg.c.9.1[,ind], 10)); sg.c.9.1 <- sg.c.9.1[,-rmove]; names(sg.c.9.1) <- vars
return(sg.c.9.1)
}
sg.c.11.1<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.11.1 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=11, m=1))    #S-G cubic, 11 point moving window
sg.c.11.1 <- data.frame(colsMean(sg.c.11.1[,ind], 10)); sg.c.11.1 <- sg.c.11.1[,-rmove]; names(sg.c.11.1) <- vars
return(sg.c.11.1)
}

# SG Quadratic 2nd Derivative
sg.q.3.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.3.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=3, m=2))    #S-G quadratic, 3 point moving window
sg.q.3.2 <- data.frame(colsMean(sg.q.3.2[,ind], 10)); sg.q.3.2 <- sg.q.3.2[,-rmove]; names(sg.q.3.2) <- vars
return(sg.q.3.2)
}
sg.q.5.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.5.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=5, m=2))    #S-G quadratic, 5 point moving window
sg.q.5.2 <- data.frame(colsMean(sg.q.5.2[,ind], 10)); sg.q.5.2 <- sg.q.5.2[,-rmove]; names(sg.q.5.2) <- vars
return(sg.q.5.2)
}
sg.q.7.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.7.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=7, m=2))    #S-G quadratic, 7 point moving window
sg.q.7.2 <- data.frame(colsMean(sg.q.7.2[,ind], 10)); sg.q.7.2 <- sg.q.7.2[,-rmove]; names(sg.q.7.2) <- vars
return(sg.q.7.2)
}
sg.q.9.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.9.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=9, m=2))    #S-G quadratic, 9 point moving window
sg.q.9.2 <- data.frame(colsMean(sg.q.9.2[,ind], 10)); sg.q.9.2 <- sg.q.9.2[,-rmove]; names(sg.q.9.2) <- vars
return(sg.q.9.2)
}
sg.q.11.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.q.11.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=2, n=11, m=2))    #S-G quadratic, 11 point moving window
sg.q.11.2 <- data.frame(colsMean(sg.q.11.2[,ind], 10)); sg.q.11.2 <- sg.q.11.2[,-rmove]; names(sg.q.11.2) <- vars
return(sg.q.11.2)
}

# SG Cubic 2nd Derivative
sg.c.5.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.5.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=5, m=2))    #S-G cubic, 5 point moving window
sg.c.5.2 <- data.frame(colsMean(sg.c.5.2[,ind], 10)); sg.c.5.2 <- sg.c.5.2[,-rmove]; names(sg.c.5.2) <- vars
return(sg.c.5.2)
}
sg.c.7.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.7.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=7, m=2))    #S-G cubic, 7 point moving window
sg.c.7.2 <- data.frame(colsMean(sg.c.7.2[,ind], 10)); sg.c.7.2 <- sg.c.7.2[,-rmove]; names(sg.c.7.2) <- vars
return(sg.c.7.2)
}
sg.c.9.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.9.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=9, m=2))    #S-G cubic, 9 point moving window
sg.c.9.2 <- data.frame(colsMean(sg.c.9.2[,ind], 10)); sg.c.9.2 <- sg.c.9.2[,-rmove]; names(sg.c.9.2) <- vars
return(sg.c.9.2)
}
sg.c.11.2<-function(vnir.matrix) {
  rm1 <- c(seq.int(1,46,1)); rm2 <- c(seq.int(637,666,1)); rm3 <- c(seq.int(1437,1466,1)); rm4 <- c(seq.int(2127,2151,1))     #for creating formulas with out including detector artifact areas ~1000nm and ~1800nm
  
  ind <- 7:2146  ## indices that will be used in aggregating over columns
  rmove <- c(1:4, 64:66, 144:146, 213:214)  ## remove these bands
  vars <- names(vnir.matrix[,seq(11, 2141, 10)]); vars <- vars[-rmove]  ## spectral band names
sg.c.11.2 <- t(apply(vnir.matrix, MARGIN=1, FUN=sgolayfilt, p=3, n=11, m=2))    #S-G cubic, 11 point moving window
sg.c.11.2 <- data.frame(colsMean(sg.c.11.2[,ind], 10)); sg.c.11.2 <- sg.c.11.2[,-rmove]; names(sg.c.11.2) <- vars
return(sg.c.11.2)
}