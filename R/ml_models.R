#' Soil and Spectral Data Preprocessing for Model Training
#'
#' These functions fit predictive models for soil properties using VNIR spectral data.
#' Each function applies a specific machine learning method:
#' \itemize{
#'   \item \code{pcr_preprocess()} – Principal Component Regression (PCR)
#'   \item \code{plsr_preprocess()} – Partial Least Squares Regression (PLSR)
#'   \item \code{lasso_preprocess()} – LASSO regression
#'   \item \code{rf_preprocess()} – Random Forest regression
#'   \item \code{cubist_preprocess()} – Cubist regression
#' }
#'
#' All functions use the same workflow:
#' \enumerate{
#'   \item Combine the selected soil property with preprocessed spectra.
#'   \item Split data into calibration and validation sets (using sample indices).
#'   \item Fit the chosen model across multiple calibration/validation partitions.
#'   \item Generate predictions and compute performance metrics (MSD-based).
#' }
#'
#' @param soil A data frame of soil properties. Must include the target soil variable.
#' @param vnir.matrix A numeric matrix of VNIR spectral data.
#' @param j A list of index vectors specifying calibration sample sets
#'   (e.g., from \code{\link{merge_of_lab_and_spectrum}}).
#' @param preprocess A preprocessing function to apply to the spectral data
#'   (e.g., smoothing, normalization).
#' @param type_of_soil An integer index selecting which soil property column to model.
#'
#' @return A list of MSD metric objects for calibration and validation sets,
#' specific to the fitted model.
#'
#' @examples
#' \donttest{
#' # Example with PCR
#' results_pcr <- pcr_preprocess(soil, vnir.matrix, j, preprocess = scale, type_of_soil = 2)
#'
#' # Example with Random Forest
#' results_rf <- rf_preprocess(soil, vnir.matrix, j, preprocess = scale, type_of_soil = 2)
#' }
#'
#' @seealso \code{\link{merge_of_lab_and_spectrum}}, \code{\link{ml_f}}
#'
#' @name soil_preprocess

#' @rdname soil_preprocess
#' @export
pcr_preprocess<-function(soil,vnir.matrix,j,preprocess,type_of_soil) {
  i <- type_of_soil     #Starting soil variable column
  soil.var<-names(soil)
  num.trees <- c("501", "1001", "1501", "2001")
  num <- c(501, 1001, 1501, 2001)
  nt <- 4 # for temp storage lists
  
  #for(i in 2:length(soil)){
  #fmla.strg <- paste(soil.var[i],"~", paste(vars, collapse= "+"))  ##regression formula
  #var.fmla <- as.formula(fmla.strg)
  
  # Main Storage Lists (will be saved as workspace files)
  #rf.list <- list(rep(NA, 41))          #stores randomForest models
  metric.list <- list(rep(NA, 41))      #stores cal/val metrics
  #pred.list <- list(rep(NA,41))         #stores observed and predicted values of validation
  #cval.list <- list(rep(NA, 41))        #stores observed and predicted values of calibration
  
  k <- 1   # Storage lists index
  ## 1 Raw Reflectance
  cal<-list()
  val<-list()
  prop.vnir <- cbind(soil[,i], preprocess(vnir.matrix)); names(prop.vnir)[1] <- soil.var[i]    # Combine soil and VNIR data
  for (jj in 1:4) {cal[[jj]] <- prop.vnir[j[[jj]],]; val[[jj]] <- prop.vnir[-j[[jj]],]}                                      # Create calibration and validation datasets 
  # Temp storage lists
  mod.list <- list(rep(NA, nt))  # Models
  msd.list <- list(rep(NA, nt))  # MSD comps
  cal.list <- list(rep(NA, nt))  # Calibration data and predictions
  val.list <- list(rep(NA, nt))  # Validation data and predictions
  # Models
  for(l in 1:nt){
    mod.list[[l]] <- pcr(cal[[l]][,1]~.,data = cal[[l]][,1:203],scale =TRUE, validation = "CV"); names(mod.list)[[l]] <- num.trees[l]
  }
  # Predictions and MSD comps
  for(m in 1:nt){
    cal.pre <- as.matrix(predict(mod.list[[m]],newdata = cal[[m]]))                                                       # Calibration prediction
    val.pre <- as.matrix(predict(mod.list[[m]],newdata = val[[m]]))                                      # Validation prediction
    val.value <- data.frame(obs=val[[m]][,1], pred=val.pre)                                               # Extract obs and predicted values for validation set
    cal.value <- data.frame(obs=cal[[m]][,1], pred=cal.pre)                                                # Extract obs and predicted values for calibration set
    msd.val <- msd.comp(y=val.value$obs, yhat=val.value$pred)                                         # Calculate metrics for validation set
    msd.cal <- msd.comp(y=cal.value$obs, yhat=cal.value$pred)                                         # Calculate metrics for calibration set
    msd.rf <- cbind(cal=msd.cal, val=msd.val)                                                         # Combine metrics
    msd.list[[m]] <- msd.rf; names(msd.list)[[m]] <- num.trees[m]                                     # Store model to list
    val.list[[m]] <- val.value; names(val.list)[[m]] <- num.trees[m]                                  # Store validation values to list
    cal.list[[m]] <- cal.value; names(cal.list)[[m]] <- num.trees[m]                                  # Store calibration values to list
  }
  # Store results in lists
  return(msd.list)
  #}
}

#' @rdname soil_preprocess
#' @export
plsr_preprocess<-function(soil,vnir.matrix,j,preprocess,type_of_soil) {
  i <- type_of_soil     #Starting soil variable column
  soil.var<-names(soil)
  num.trees <- c("501", "1001", "1501", "2001")
  num <- c(501, 1001, 1501, 2001)
  nt <- 4 # for temp storage lists
  
  #for(i in 2:length(soil)){
  #fmla.strg <- paste(soil.var[i],"~", paste(vars, collapse= "+"))  ##regression formula
  #var.fmla <- as.formula(fmla.strg)
  
  # Main Storage Lists (will be saved as workspace files)
  rf.list <- list(rep(NA, 41))          #stores randomForest models
  metric.list <- list(rep(NA, 41))      #stores cal/val metrics
  pred.list <- list(rep(NA,41))         #stores observed and predicted values of validation
  cval.list <- list(rep(NA, 41))        #stores observed and predicted values of calibration
  
  k <- 1   # Storage lists index
  ## 1 Raw Reflectance
  cal<-list()
  val<-list()
  prop.vnir <- cbind(soil[,i], preprocess(vnir.matrix)); names(prop.vnir)[1] <- soil.var[i]    # Combine soil and VNIR data
  for (jj in 1:4) {cal[[jj]] <- prop.vnir[j[[jj]],]; val[[jj]] <- prop.vnir[-j[[jj]],]}                                      # Create calibration and validation datasets 
  # Temp storage lists
  mod.list <- list(rep(NA, nt))  # Models
  msd.list <- list(rep(NA, nt))  # MSD comps
  cal.list <- list(rep(NA, nt))  # Calibration data and predictions
  val.list <- list(rep(NA, nt))  # Validation data and predictions
  # Models
  for(l in 1:nt){
    mod.list[[l]] <- plsr(cal[[l]][,1]~.,data = cal[[l]][,1:203],scale =TRUE, validation = "CV"); names(mod.list)[[l]] <- num.trees[l]
  }
  # Predictions and MSD comps
  for(m in 1:nt){
    cal.pre <- as.matrix(predict(mod.list[[m]],newdata = cal[[m]]))                                                       # Calibration prediction
    val.pre <- as.matrix(predict(mod.list[[m]],newdata = val[[m]]))                                      # Validation prediction
    val.value <- data.frame(obs=val[[m]][,1], pred=val.pre)                                               # Extract obs and predicted values for validation set
    cal.value <- data.frame(obs=cal[[m]][,1], pred=cal.pre)                                                # Extract obs and predicted values for calibration set
    msd.val <- msd.comp(y=val.value$obs, yhat=val.value$pred)                                         # Calculate metrics for validation set
    msd.cal <- msd.comp(y=cal.value$obs, yhat=cal.value$pred)                                         # Calculate metrics for calibration set
    msd.rf <- cbind(cal=msd.cal, val=msd.val)                                                         # Combine metrics
    msd.list[[m]] <- msd.rf; names(msd.list)[[m]] <- num.trees[m]                                     # Store model to list
    val.list[[m]] <- val.value; names(val.list)[[m]] <- num.trees[m]                                  # Store validation values to list
    cal.list[[m]] <- cal.value; names(cal.list)[[m]] <- num.trees[m]                                  # Store calibration values to list
  }
  # Store results in lists
  return(msd.list)
  #}
}
#' @rdname soil_preprocess
#' @export
lasso_preprocess<-function(soil,vnir.matrix,j,preprocess,type_of_soil) {
  i <- type_of_soil     #Starting soil variable column
  soil.var<-names(soil)
  num.trees <- c("501", "1001", "1501", "2001")
  num <- c(501, 1001, 1501, 2001)
  nt <- 4 # for temp storage lists
  
  #for(i in 2:length(soil)){
  #fmla.strg <- paste(soil.var[i],"~", paste(vars, collapse= "+"))  ##regression formula
  #var.fmla <- as.formula(fmla.strg)
  
  # Main Storage Lists (will be saved as workspace files)
  rf.list <- list(rep(NA, 41))          #stores randomForest models
  metric.list <- list(rep(NA, 41))      #stores cal/val metrics
  pred.list <- list(rep(NA,41))         #stores observed and predicted values of validation
  cval.list <- list(rep(NA, 41))        #stores observed and predicted values of calibration
  
  k <- 1   # Storage lists index
  ## 1 Raw Reflectance
  cal<-list()
  val<-list()
  prop.vnir <- cbind(soil[,i], preprocess(vnir.matrix)); names(prop.vnir)[1] <- soil.var[i]    # Combine soil and VNIR data
  for (jj in 1:4) {cal[[jj]] <- prop.vnir[j[[jj]],]; val[[jj]] <- prop.vnir[-j[[jj]],]}                                      # Create calibration and validation datasets 
  # Temp storage lists
  mod.list <- list(rep(NA, nt))  # Models
  msd.list <- list(rep(NA, nt))  # MSD comps
  cal.list <- list(rep(NA, nt))  # Calibration data and predictions
  val.list <- list(rep(NA, nt))  # Validation data and predictions
  # Models
  for(l in 1:nt){
    cv_model<-cv.glmnet(as.matrix(cal[[l]][, 2:203]), cal[[l]][, 1], alpha = 1);best_lambda <- cv_model$lambda.min;mod.list[[l]] <- glmnet(cal[[l]][,2:203], cal[[l]][,1], alpha = 1, lambda = best_lambda)
  }
  # Predictions and MSD comps
  for(m in 1:nt){
    cal.pre <- as.matrix(predict(mod.list[[m]],newx =  as.matrix(cal[[m]][,2:203])))                                                       # Calibration prediction
    val.pre <- as.matrix(predict(mod.list[[m]],newx =  as.matrix(val[[m]][,2:203])))                                      # Validation prediction
    val.value <- data.frame(obs=val[[m]][,1], pred=val.pre)                                               # Extract obs and predicted values for validation set
    cal.value <- data.frame(obs=cal[[m]][,1], pred=cal.pre)                                                # Extract obs and predicted values for calibration set
    msd.val <- msd.comp(y=val.value$obs, yhat=val.value$s0)                                         # Calculate metrics for validation set
    msd.cal <- msd.comp(y=cal.value$obs, yhat=cal.value$s0)                                         # Calculate metrics for calibration set
    msd.rf <- cbind(cal=msd.cal, val=msd.val)                                                         # Combine metrics
    msd.list[[m]] <- msd.rf; names(msd.list)[[m]] <- num.trees[m]                                     # Store model to list
    val.list[[m]] <- val.value; names(val.list)[[m]] <- num.trees[m]                                  # Store validation values to list
    cal.list[[m]] <- cal.value; names(cal.list)[[m]] <- num.trees[m]                                  # Store calibration values to list
  }
  #}
  # Store results in lists
  return(msd.list)
  
}
#' @rdname soil_preprocess
#' @export
rf_preprocess<-function(soil,vnir.matrix,j,preprocess,type_of_soil) {
  i <- type_of_soil     #Starting soil variable column
  soil.var<-names(soil)
  num.trees <- c("501", "1001", "1501", "2001")
  num <- c(501, 1001, 1501, 2001)
  nt <- 4 # for temp storage lists
  
  #for(i in 2:length(soil)){
  #fmla.strg <- paste(soil.var[i],"~", paste(vars, collapse= "+"))  ##regression formula
  #var.fmla <- as.formula(fmla.strg)
  
  # Main Storage Lists (will be saved as workspace files)
  rf.list <- list(rep(NA, 41))          #stores randomForest models
  metric.list <- list(rep(NA, 41))      #stores cal/val metrics
  pred.list <- list(rep(NA,41))         #stores observed and predicted values of validation
  cval.list <- list(rep(NA, 41))        #stores observed and predicted values of calibration
  cal<-list()
  val<-list()
  k <- 1   # Storage lists index
  ## 1 Raw Reflectance
  prop.vnir <- cbind(soil[,i], preprocess(vnir.matrix)); names(prop.vnir)[1] <- soil.var[i]    # Combine soil and VNIR data
  for (jj in 1:4) {cal[[jj]] <- prop.vnir[j[[jj]],]; val[[jj]] <- prop.vnir[-j[[jj]],]}                                       # Create calibration and validation datasets 
  # Temp storage lists
  mod.list <- list(rep(NA, nt*4))  # Models
  msd.list <- list(rep(NA, nt*4))  # MSD comps
  cal.list <- list(rep(NA, nt*4))  # Calibration data and predictions
  val.list <- list(rep(NA, nt*4))  # Validation data and predictions
  # Models
  lll<-1
  lll<-1;for(l in 1:nt){
    for (ll in 1:4) {
      mod.list[[lll]] <- randomForest(x=cal[[ll]][, 2:203],y=cal[[ll]][,1],data = cal[ll], importance = T, proximity = T, ntree = num[l]); names(mod.list)[[lll]] <- num.trees[l]
      lll<-lll+1
    }
    
  }
  # Predictions and MSD comps
  for(m in 1:nt){
    cal.pre <- as.matrix(mod.list[[m]]$predicted)                                                      # Calibration prediction
    val.pre <- as.matrix(predict(mod.list[[m]], newdata = val[[m]]))                                       # Validation prediction
    val.value <- data.frame(obs=val[[m]][,1], pred=val.pre)                                                # Extract obs and predicted values for validation set
    cal.value <- data.frame(obs=cal[[m]][,1], pred=cal.pre)                                                # Extract obs and predicted values for calibration set
    msd.val <- msd.comp(y=val.value$obs, yhat=val.value$pred)                                         # Calculate metrics for validation set
    msd.cal <- msd.comp(y=cal.value$obs, yhat=cal.value$pred)                                         # Calculate metrics for calibration set
    msd.rf <- cbind(cal=msd.cal, val=msd.val)                                                         # Combine metrics
    msd.list[[m]] <- msd.rf; names(msd.list)[[m]] <- num.trees[m]                                     # Store model to list
    val.list[[m]] <- val.value; names(val.list)[[m]] <- num.trees[m]                                  # Store validation values to list
    cal.list[[m]] <- cal.value; names(cal.list)[[m]] <- num.trees[m]                                  # Store calibration values to list
  }
  for(m in (nt+1):(2*nt)){
    cal.pre <- as.matrix(mod.list[[m]]$predicted)                                                      # Calibration prediction
    val.pre <- as.matrix(predict(mod.list[[m]], newdata = val[[m-nt]]))                                       # Validation prediction
    val.value <- data.frame(obs=val[[m-nt]][,1], pred=val.pre)                                                # Extract obs and predicted values for validation set
    cal.value <- data.frame(obs=cal[[m-nt]][,1], pred=cal.pre)                                                # Extract obs and predicted values for calibration set
    msd.val <- msd.comp(y=val.value$obs, yhat=val.value$pred)                                         # Calculate metrics for validation set
    msd.cal <- msd.comp(y=cal.value$obs, yhat=cal.value$pred)                                         # Calculate metrics for calibration set
    msd.rf <- cbind(cal=msd.cal, val=msd.val)                                                         # Combine metrics
    msd.list[[m]] <- msd.rf; names(msd.list)[[m]] <- num.trees[m-nt]                                     # Store model to list
    val.list[[m]] <- val.value; names(val.list)[[m]] <- num.trees[m-nt]                                  # Store validation values to list
    cal.list[[m]] <- cal.value; names(cal.list)[[m]] <- num.trees[m-nt]                                  # Store calibration values to list
  }
  for(m in (2*nt+1):(3*nt)){
    cal.pre <- as.matrix(mod.list[[m]]$predicted)                                                      # Calibration prediction
    val.pre <- as.matrix(predict(mod.list[[m]], newdata = val[[m-2*nt]]))                                       # Validation prediction
    val.value <- data.frame(obs=val[[m-2*nt]][,1], pred=val.pre)                                                # Extract obs and predicted values for validation set
    cal.value <- data.frame(obs=cal[[m-2*nt]][,1], pred=cal.pre)                                                # Extract obs and predicted values for calibration set
    msd.val <- msd.comp(y=val.value$obs, yhat=val.value$pred)                                         # Calculate metrics for validation set
    msd.cal <- msd.comp(y=cal.value$obs, yhat=cal.value$pred)                                         # Calculate metrics for calibration set
    msd.rf <- cbind(cal=msd.cal, val=msd.val)                                                         # Combine metrics
    msd.list[[m]] <- msd.rf; names(msd.list)[[m]] <- num.trees[m-2*nt]                                     # Store model to list
    val.list[[m]] <- val.value; names(val.list)[[m]] <- num.trees[m-2*nt]                                  # Store validation values to list
    cal.list[[m]] <- cal.value; names(cal.list)[[m]] <- num.trees[m-2*nt]                                  # Store calibration values to list
  }
  for(m in (3*nt+1):(4*nt)){
    cal.pre <- as.matrix(mod.list[[m]]$predicted)                                                      # Calibration prediction
    val.pre <- as.matrix(predict(mod.list[[m]], newdata = val[[m-3*nt]]))                                       # Validation prediction
    val.value <- data.frame(obs=val[[m-3*nt]][,1], pred=val.pre)                                                # Extract obs and predicted values for validation set
    cal.value <- data.frame(obs=cal[[m-3*nt]][,1], pred=cal.pre)                                                # Extract obs and predicted values for calibration set
    msd.val <- msd.comp(y=val.value$obs, yhat=val.value$pred)                                         # Calculate metrics for validation set
    msd.cal <- msd.comp(y=cal.value$obs, yhat=cal.value$pred)                                         # Calculate metrics for calibration set
    msd.rf <- cbind(cal=msd.cal, val=msd.val)                                                         # Combine metrics
    msd.list[[m]] <- msd.rf; names(msd.list)[[m]] <- num.trees[m-3*nt]                                     # Store model to list
    val.list[[m]] <- val.value; names(val.list)[[m]] <- num.trees[m-3*nt]                                  # Store validation values to list
    cal.list[[m]] <- cal.value; names(cal.list)[[m]] <- num.trees[m-3*nt]                                  # Store calibration values to list
  }
  # Store results in lists
  #rf.list[[k]] <- mod.list; metric.list[[k]] <- msd.list; pred.list[[k]] <- val.list; cval.list[[k]] <- cal.list    # Save temp lists within storage lists
  #tech <- "RAW"; names(rf.list)[[k]] <- tech; names(metric.list)[[k]] <- tech; names(pred.list)[[k]] <- tech; names(cval.list)[[k]] <- tech    # Names
  #rm(prop.vnir, mod.list, cal.pre, val.pre, val.value, cal.value, msd.val, msd.cal, msd.rf, m, msd.list, cal.list, val.list)    # Remove data structures
  #k <- k+1   # Increment storage list index
  #}
  return(msd.list)
}
#' @rdname soil_preprocess
#' @export
cubist_preprocess<-function(soil,vnir.matrix,j,preprocess,type_of_soil) {
  i <- type_of_soil     #Starting soil variable column
  soil.var<-names(soil)
  num.trees <- c("501", "1001", "1501", "2001")
  num <- c(501, 1001, 1501, 2001)
  nt <- 4 # for temp storage lists
  
  #for(i in 2:length(soil)){
  #fmla.strg <- paste(soil.var[i],"~", paste(vars, collapse= "+"))  ##regression formula
  #var.fmla <- as.formula(fmla.strg)
  
  # Main Storage Lists (will be saved as workspace files)
  rf.list <- list(rep(NA, 41))          #stores randomForest models
  metric.list <- list(rep(NA, 41))      #stores cal/val metrics
  pred.list <- list(rep(NA,41))         #stores observed and predicted values of validation
  cval.list <- list(rep(NA, 41))        #stores observed and predicted values of calibration
  
  k <- 1   # Storage lists index
  ## 1 Raw Reflectance
  cal<-list()
  val<-list()
  prop.vnir <- cbind(soil[,i], preprocess(vnir.matrix)); names(prop.vnir)[1] <- soil.var[i]    # Combine soil and VNIR data
  for (jj in 1:4) {cal[[jj]] <- prop.vnir[j[[jj]],]; val[[jj]] <- prop.vnir[-j[[jj]],]}                                      # Create calibration and validation datasets 
  # Temp storage lists
  mod.list <- list(rep(NA, nt))  # Models
  msd.list <- list(rep(NA, nt))  # MSD comps
  cal.list <- list(rep(NA, nt))  # Calibration data and predictions
  val.list <- list(rep(NA, nt))  # Validation data and predictions
  # Models
  for(l in 1:nt){
    mod.list[[l]]<-cubist(x =cal[[l]][,2:203] , y =cal[[l]][,1])
  }
  # Predictions and MSD comps
  for(m in 1:nt){
    cal.pre <- as.matrix(predict(mod.list[[m]],newdata = cal[[m]]))                                                       # Calibration prediction
    val.pre <- as.matrix(predict(mod.list[[m]],newdata = val[[m]]))                                      # Validation prediction
    val.value <- data.frame(obs=val[[m]][,1], pred=val.pre)                                               # Extract obs and predicted values for validation set
    cal.value <- data.frame(obs=cal[[m]][,1], pred=cal.pre)                                                # Extract obs and predicted values for calibration set
    msd.val <- msd.comp(y=val.value$obs, yhat=val.value$pred)                                         # Calculate metrics for validation set
    msd.cal <- msd.comp(y=cal.value$obs, yhat=cal.value$pred)                                         # Calculate metrics for calibration set
    msd.rf <- cbind(cal=msd.cal, val=msd.val)                                                         # Combine metrics
    msd.list[[m]] <- msd.rf; names(msd.list)[[m]] <- num.trees[m]                                     # Store model to list
    val.list[[m]] <- val.value; names(val.list)[[m]] <- num.trees[m]                                  # Store validation values to list
    cal.list[[m]] <- cal.value; names(cal.list)[[m]] <- num.trees[m]                                  # Store calibration values to list
  }
  # Store results in lists
  return(msd.list)
  #}
}