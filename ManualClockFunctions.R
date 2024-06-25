### Re-Writing Manual Clock Computation Codes ### ----

### Reading in Libraries ### ----
library(methylCIPHER)
library(DunedinPACE)
GrimAge2_data <- readRDS("DNAmGrimAge2_final.Rds")

#######################################################
### Horvath1 Clock Function ### ----
#######################################################

calcHorvath1Edit <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){
  
  #######################
  ### Read in the Data###
  #######################
  if(!exists("anti.trafo") || !exists("trafo")){
    stop("Calculation of Horvath1 requires that the functions trafo and anti.trafo is loaded. \n Ensure that you loaded the entire calcAllCpGClocks library")
  }
  
  #data("Horvath1_CpGs")
  
  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################
  
  CpGCheck <- length(Horvath1_CpGs$CpGmarker) == sum(Horvath1_CpGs$CpGmarker %in% colnames(DNAm))
  
  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################
  
  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){
    
    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")
    
  } else if(CpGCheck == T | imputation == F){
    
    present <- Horvath1_CpGs$CpGmarker %in% colnames(DNAm)
    
    betas <- DNAm[,na.omit(match(Horvath1_CpGs$CpGmarker,colnames(DNAm)))]
    tt <- sweep(betas, MARGIN = 2, Horvath1_CpGs$CoefficientTraining[present], `*`)
    
    Horvath1 <- as.numeric(anti.trafo(rowSums(tt,na.rm=T)+0.696))
    if(is.null(pheno)){
      Horvath1
    } else{
      pheno$Horvath1 <- Horvath1
      pheno
    }
    
  } else {
    message("Imputation of mean CpG Values occured for Horvath1")
    missingCpGs <- Horvath1_CpGs$CpGmarker[!(Horvath1_CpGs$CpGmarker %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))
    
    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)
    
    betas <- DNAm[,match(Horvath1_CpGs$CpGmarker,colnames(DNAm))]
    tt <- sweep(betas, MARGIN = 2, Horvath1_CpGs$CoefficientTraining, `*`)
    
    Horvath1 <- as.numeric(anti.trafo(rowSums(tt,na.rm=T)+0.696))
    if(is.null(pheno)){
      Horvath1
    } else{
      pheno$Horvath1 <- Horvath1
      pheno
    }
    
  }
  
}




#######################################################
### Hannum Clock Function ### ----
#######################################################

calcHannumEdit <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){
  
  #######################
  ### Read in the Data###
  #######################
  
  #data("Hannum_CpGs")
  
  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################
  
  CpGCheck <- length(Hannum_CpGs$Marker) == sum(Hannum_CpGs$Marker %in% colnames(DNAm))
  
  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################
  
  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){
    
    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")
    
  } else if(CpGCheck == T | imputation == F){
    
    present <- Hannum_CpGs$Marker %in% colnames(DNAm)
    
    betas <- DNAm[,na.omit(match(Hannum_CpGs$Marker,colnames(DNAm)))]
    tt <- rowSums(sweep(as.matrix(betas) ,MARGIN = 2, Hannum_CpGs$Coefficient[present], `*`), na.rm = T)
    
    if(is.null(pheno)){
      tt
    } else{
      pheno$Hannum <- tt
      pheno
    }
    
  } else {
    message("Imputation of mean CpG Values occured for Hannum")
    missingCpGs <- Hannum_CpGs$Marker[!(Hannum_CpGs$Marker %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))
    
    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)
    
    betas <- DNAm[,match(Hannum_CpGs$Marker,colnames(DNAm))]
    tt <- rowSums(sweep(betas ,MARGIN = 2, Hannum_CpGs$Coefficient, `*`), na.rm = TRUE)
    
    if(is.null(pheno)){
      tt
    } else{
      pheno$Hannum <- tt
      pheno
    }
    
  }
  
}




#######################################################
### Horvath2 Clock Function ### ----
#######################################################
calcHorvath2Edit <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){
  
  #######################
  ### Read in the Data###
  #######################
  
  if(!exists("anti.trafo") || !exists("trafo")){
    stop("Calculation of Horvath2 requires that the functions trafo and anti.trafo is loaded. \n Ensure that you loaded the entire calcAllCpGClocks library")
  }
  
  #data("Horvath2_CpGs")
  
  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################
  
  CpGCheck <- length(Horvath2_CpGs$ID) == sum(Horvath2_CpGs$ID %in% colnames(DNAm))
  
  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################
  
  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){
    
    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")
    
  } else if(CpGCheck == T | imputation == F){
    
    present <- Horvath2_CpGs$ID[-1] %in% colnames(DNAm)
    
    betas <- DNAm[,na.omit(match(Horvath2_CpGs$ID[-1],colnames(DNAm)))]
    tt <- sweep(as.matrix(betas), MARGIN = 2, Horvath2_CpGs$Coef[-1][present], `*`)
    
    Horvath2 <- as.numeric(anti.trafo(rowSums(tt,na.rm=T)-0.447119319))
    if(is.null(pheno)){
      Horvath2
    } else{
      pheno$Horvath2 <- Horvath2
      pheno
    }
    
  } else {
    message("Imputation of mean CpG Values occured for Horvath2")
    missingCpGs <- Horvath2_CpGs$ID[-1][!(Horvath2_CpGs$ID[-1] %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))
    
    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)
    
    betas <- DNAm[,na.omit(match(Horvath2_CpGs$ID[-1],colnames(DNAm)))]
    tt <- sweep(betas, MARGIN = 2, Horvath2_CpGs$Coef[-1], `*`)
    
    Horvath2 <- as.numeric(anti.trafo(rowSums(tt,na.rm=T)-0.447119319))
    if(is.null(pheno)){
      Horvath2
    } else{
      pheno$Horvath2 <- Horvath2
      pheno
    }
    
  }
  
}






#######################################################
### PhenoAge Clock Function ### ----
#######################################################
calcPhenoAgeEdit <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){
  
  #######################
  ### Read in the Data###
  #######################
  
  #data("PhenoAge_CpGs")
  
  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################
  CpGCheck <- length(PhenoAge_CpGs$CpG) == sum(PhenoAge_CpGs$CpG %in% colnames(DNAm))
  
  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################
  
  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){
    
    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")
    
  } else if(CpGCheck == T | imputation == F){
    
    present <- PhenoAge_CpGs$CpG %in% colnames(DNAm)
    
    betas <- DNAm[,na.omit(match(PhenoAge_CpGs$CpG,colnames(DNAm)))]
    tt <- rep(0, dim(DNAm)[1])
    
    tt <- rowSums(sweep(as.matrix(betas), MARGIN = 2, PhenoAge_CpGs$Weight[present],`*`), na.rm = T) + 60.664
    
    if(is.null(pheno)){
      tt
    } else{
      pheno$PhenoAge <- tt
      pheno
    }
    
  } else {
    message("Imputation of mean CpG Values occured for PhenoAge")
    missingCpGs <- PhenoAge_CpGs$CpG[!(PhenoAge_CpGs$CpG %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))
    
    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)
    
    betas <- DNAm[,match(PhenoAge_CpGs$CpG,colnames(DNAm))]
    betas <- betas %>% mutate_all(~ifelse(is.nan(.), NA, .))
    betas <- betas %>% mutate_all(~ifelse(is.na(.), 0, .))
    tt <- rep(0, dim(DNAm)[1])
    tt <- rowSums(sweep(betas, MARGIN = 2, PhenoAge_CpGs$Weight,`*`)) + 60.664
    
    if(is.null(pheno)){
      tt
    } else{
      pheno$PhenoAge <- tt
      pheno
    }
    
  }
  
}




#######################################################
### PedBE Clock Function ### ----
#######################################################
calcPEDBEEdit <- function(DNAm, pheno = NULL, CpGImputation = NULL, imputation = F){
  
  #######################
  ### Read in the Data###
  #######################
  if(!exists("anti.trafo") || !exists("trafo")){
    stop("Calculation of PEDBE requires that the functions trafo and anti.trafo is loaded. \n Ensure that you loaded the entire calcAllCpGClocks library")
  }
  
  #data("PEDBE_CpGs")
  
  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################
  
  CpGCheck <- length(PEDBE_CpGs$ID) == sum(PEDBE_CpGs$ID %in% colnames(DNAm))
  
  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################
  
  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){
    
    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")
    
  } else if(CpGCheck == T | imputation == F){
    
    present <- PEDBE_CpGs$ID %in% colnames(DNAm)
    
    betas <- DNAm[,na.omit(match(PEDBE_CpGs$ID,colnames(DNAm)))]
    tt <- sweep(betas, MARGIN = 2, PEDBE_CpGs$Coef[present], `*`)
    
    PEDBE <- as.numeric(anti.trafo(rowSums(tt,na.rm=T)-2.10))
    if(is.null(pheno)){
      PEDBE
    } else{
      pheno$PEDBE <- PEDBE
      pheno
    }
    
  } else {
    message("Imputation of mean CpG Values occured for PEDBE")
    missingCpGs <- PEDBE_CpGs$ID[!(PEDBE_CpGs$ID %in% colnames(DNAm))]
    tempDNAm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs))
    
    for(j in 1:length(missingCpGs)){
      meanVals <- CpGImputation[match(missingCpGs[j],names(CpGImputation))]
      tempDNAm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    colnames(tempDNAm) <- missingCpGs
    DNAm <- cbind(DNAm,tempDNAm)
    
    betas <- DNAm[,match(PEDBE_CpGs$ID,colnames(DNAm))]
    tt <- sweep(betas, MARGIN = 2, PEDBE_CpGs$Coef, `*`)
    
    PEDBE <- as.numeric(anti.trafo(rowSums(tt,na.rm=T)-2.10))
    if(is.null(pheno)){
      PEDBE
    } else{
      pheno$PEDBE <- PEDBE
      pheno
    }
    
  }
  
}



#######################################################
### DunedinPACE Function ### ----
#######################################################

calcPACEEdit = function( betas, proportionOfProbesRequired=0.8, GoldenStandard) {
  requireNamespace("preprocessCore")
  
  if( any(grepl("TC|BC", rownames(betas))) )
  {
    # Print message to user
    print("This looks like EPICv2 array data. If EPICv2, DunedinPACE will lower the proportionOfProbesRequired to 0.7 and proceed with missing probes present on 450k and EPICv1. Averaging the opposite strand replicates may take a bit more time.")
    
    # Proceed with averaging: Code taken from ENmix package.
    cgid = sapply(strsplit(rownames(betas), split = "_"),  unlist)[1, ]
    dupcg = unique(cgid[duplicated(cgid)])
    betas2 = betas[cgid %in% dupcg, ]
    cid = sapply(strsplit(rownames(betas2), split = "_"),   unlist)[1, ]
    betas2 = aggregate(betas2, by = list(cid), FUN = function(x) mean(x, na.rm = TRUE))
    rownames(betas2) = betas2[, 1]
    betas2 = as.matrix(betas2[, -1])
    betas = betas[!(cgid %in% dupcg), ]
    rownames(betas) = sapply(strsplit(rownames(betas), split = "_"), unlist)[1, ]
    betas <-rbind(betas, betas2)
    
    # Set proportionOfProbesRequired to 0.7 if condition is TRUE
    proportionOfProbesRequired <- 0.7
  }
  
  # loop through models
  model_results <- lapply(mPACE_Models$model_names, function(model_name) {
    # make sure it has been converted to a matrix
    if( !is.numeric(as.matrix(betas)) ) { stop("betas matrix/data.frame is not numeric!") }
    probeOverlap <- length(which(rownames(betas) %in% mPACE_Models$model_probes[[model_name]])) / length(mPACE_Models$model_probes[[model_name]])
    probeOverlap_background <- length(which(rownames(betas) %in% mPACE_Models$gold_standard_probes[[model_name]])) / length(mPACE_Models$gold_standard_probes[[model_name]])
    # make sure enough of the probes are present in the data file
    if( probeOverlap < proportionOfProbesRequired | probeOverlap_background < proportionOfProbesRequired ) {
      result <- rep(NA, ncol(betas))
      names(result) <- colnames(betas)
      result
    } else {
      # Work with a numeric matrix of betas
      betas.mat <- as.matrix(betas[which(rownames(betas) %in% mPACE_Models$gold_standard_probes[[model_name]]),])
      # If probes don't exist, we'll add them as rows of values based on their mean in the gold standard dataset
      probesNotInMatrix <- mPACE_Models$gold_standard_probes[[model_name]][which(mPACE_Models$gold_standard_probes[[model_name]] %in% rownames(betas.mat) == F)]
      if( length(probesNotInMatrix) > 0 ) {
        for( probe in probesNotInMatrix ) {
          tmp.mat <- matrix(0, nrow=1, ncol=ncol(betas.mat))
          rownames(tmp.mat) <- probe
          colnames(tmp.mat) <- colnames(betas.mat)
          tmp.mat[probe,] <- rep(GoldenStandard[probe], ncol(tmp.mat))
          betas.mat <- rbind(betas.mat, tmp.mat)
        }
      }
      
      # Identify samples with too many missing probes and remove them from the matrix
      samplesToRemove <- colnames(betas.mat)[which(apply(betas.mat, 2, function(x) { 1 - ( length(which(is.na(x))) / length(x) ) < proportionOfProbesRequired}))]
      if( length(samplesToRemove) > 0 ) {
        betas.mat <- betas.mat[,-which(colnames(betas.mat) %in% samplesToRemove)]
      }
      if(ncol(betas.mat) > 0) {
        # Identify missingness on a probe level
        pctValuesPresent <- apply( betas.mat, 1, function(x) { 1 - (length(which(is.na(x))) / length(x)) } )
        # If they're missing values, but less than the proportion required, we impute to the cohort mean
        probesToAdjust <- which(pctValuesPresent < 1 & pctValuesPresent >= proportionOfProbesRequired)
        if( length(probesToAdjust) > 0 ) {
          if( length(probesToAdjust) > 1 ) {
            betas.mat[probesToAdjust,] <- t(apply( betas.mat[probesToAdjust,], 1 , function(x) {
              x[is.na(x)] = mean( x, na.rm = TRUE )
              x
            }))
          } else {
            betas.mat[probesToAdjust,which(is.na(betas.mat[probesToAdjust,]))] <- mean(betas.mat[probesToAdjust,], na.rm=T)
          }
        }
        # If they're missing too many values, everyones value gets replaced with the mean from the Dunedin cohort
        if( length(which(pctValuesPresent < proportionOfProbesRequired)) > 0 ) {
          probesToReplaceWithMean <- rownames(betas.mat)[which(pctValuesPresent < proportionOfProbesRequired)]
          for( probe in probesToReplaceWithMean ) {
            betas.mat[probe,] <- rep(GoldenStandard[probe], ncol(betas.mat))
          }
        }
        
        # Normalize the matrix to the gold standard dataset
        betas.norm <- preprocessCore::normalize.quantiles.use.target(betas.mat, target=mPACE_Models$gold_standard_means[[model_name]])
        rownames(betas.norm) <- rownames(betas.mat)
        colnames(betas.norm) <- colnames(betas.mat)
        # Calculate score:
        score = mPACE_Models$model_intercept[[model_name]] + rowSums(t(betas.norm[mPACE_Models$model_probes[[model_name]],]) %*% diag(mPACE_Models$model_weights[[model_name]]))
        names(score) <- colnames(betas.norm)
        if( length(samplesToRemove) > 0 ) {
          score.tmp <- rep(NA, length(samplesToRemove))
          names(score.tmp) <- samplesToRemove
          score <- c(score, score.tmp)
        }
        score <- score[colnames(betas)]
        score
      } else {
        result <- rep(NA, ncol(betas.mat))
        names(result) <- colnames(betas.mat)
        result
      }
    }
  })
  names(model_results) <- mPACE_Models$model_names
  model_results
}


#######################################################
### GrimAge2 Function ### ----
#######################################################

calcGrimAge2Edit <- function(DNAm, Ages, Sex, pheno = NULL, CpGImputation = NULL, imputation = F){
  
  #######################
  ### Read in the Data###
  #######################
  
  GrimAge2_data <- readRDS("DNAmGrimAge2_final.Rds")
  
  #############################################
  ### Generate List of Unique Probes Needed ###
  #############################################
  
  # All GrimAge2 Probes #
  GrimAge2_CpGs_unique <- unique(GrimAge2_data$dnam.all$var)[-c(1:2)]
  
  # Probes for Specific Subcomponents of GrimAge2 #
  DNAmadm_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmadm"][-c(1:2)]
  DNAmB2M_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmB2M"][-c(1:2)]
  DNAmCystatin_C_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmCystatin_C"][-c(1:2)]
  DNAmGDF_15_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmGDF_15"][-c(1:2)]
  DNAmleptin_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmleptin"][-c(1)]
  DNAmlog.A1C_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.A1C"][-c(1:2)]
  DNAmlog.CRP_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.CRP"][-c(1)]
  DNAmPACKYRS_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmPACKYRS"][-c(1:2)]
  DNAmpai_1_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmpai_1"][-c(1)]
  DNAmTIMP_1_CpGs <- GrimAge2_data$dnam.all$var[GrimAge2_data$dnam.all$Y.pred == "DNAmTIMP_1"][-c(1:2)]
  
  # Weights for Specific Subcomponents of GrimAge2 #
  DNAmadm_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmadm"][-c(1:2)]
  DNAmB2M_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmB2M"][-c(1:2)]
  DNAmCystatin_C_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmCystatin_C"][-c(1:2)]
  DNAmGDF_15_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmGDF_15"][-c(1:2)]
  DNAmleptin_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmleptin"][-c(1)]
  DNAmlog.A1C_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.A1C"][-c(1:2)]
  DNAmlog.CRP_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.CRP"][-c(1)]
  DNAmPACKYRS_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmPACKYRS"][-c(1:2)]
  DNAmpai_1_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmpai_1"][-c(1)]
  DNAmTIMP_1_Weights <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmTIMP_1"][-c(1:2)]
  
  ###################################################
  ### Check if all necessary CpGs are in the data ###
  ###################################################
  
  CpGCheck <- length(GrimAge2_CpGs_unique) == sum(GrimAge2_CpGs_unique %in% colnames(DNAm))
  
  ###################################################################################
  ### The calculation will be performed or an error will be thrown as appropriate ###
  ###################################################################################
  
  if(CpGCheck == F && is.null(CpGImputation) && imputation == T){
    
    stop("Need to provide of named vector of CpG Imputations; Necessary CpGs are missing!")
    
  } else if(CpGCheck == T | imputation == F){
    
    # Initialize empty dataframe for results #
    SubcomponentResults <- data.frame(DNAmGDF_15_Total = rep(NA,nrow(DNAm)),
                                      DNAmB2M_Total = rep(NA,nrow(DNAm)),
                                      DNAmCystatin_C_Total = rep(NA,nrow(DNAm)),
                                      DNAmTIMP_1_Total = rep(NA,nrow(DNAm)),
                                      DNAmadm_Total = rep(NA,nrow(DNAm)),
                                      DNAmpai_1_Total = rep(NA,nrow(DNAm)),
                                      DNAmleptin_Total = rep(NA,nrow(DNAm)),
                                      DNAmPACKYRS_Total = rep(NA,nrow(DNAm)),
                                      DNAmlog.CRP_Total = rep(NA,nrow(DNAm)),
                                      DNAmlog.A1C_Total = rep(NA,nrow(DNAm)),
                                      Age_Total = rep(NA,nrow(DNAm)),
                                      Female_Total = rep(NA,nrow(DNAm)))
    
    # DNAmadm calculations #
    present_DNAmadm <- DNAmadm_CpGs %in% colnames(DNAm)
    betas_DNAmadm <- DNAm[,na.omit(match(DNAmadm_CpGs,colnames(DNAm)))]
    tt_DNAmadm <- rowSums(sweep(as.matrix(betas_DNAmadm), MARGIN = 2, DNAmadm_Weights[present_DNAmadm], `*`), na.rm = T)
    DNAmadm_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmadm" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmadm_Total <- tt_DNAmadm + DNAmadm_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmadm"][1]
    
    # DNAmB2M calculations #
    present_DNAmB2M <- DNAmB2M_CpGs %in% colnames(DNAm)
    betas_DNAmB2M <- DNAm[,na.omit(match(DNAmB2M_CpGs,colnames(DNAm)))]
    tt_DNAmB2M <- rowSums(sweep(as.matrix(betas_DNAmB2M), MARGIN = 2, DNAmB2M_Weights[present_DNAmB2M], `*`), na.rm = T)
    DNAmB2M_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmB2M" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmB2M_Total <- tt_DNAmB2M + DNAmB2M_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmB2M"][1]
    
    # DNAmCystatin_C calculations #
    present_DNAmCystatin_C <- DNAmCystatin_C_CpGs %in% colnames(DNAm)
    betas_DNAmCystatin_C <- DNAm[,na.omit(match(DNAmCystatin_C_CpGs,colnames(DNAm)))]
    tt_DNAmCystatin_C <- rowSums(sweep(as.matrix(betas_DNAmCystatin_C), MARGIN = 2, DNAmCystatin_C_Weights[present_DNAmCystatin_C], `*`), na.rm = T)
    DNAmCystatin_C_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmCystatin_C" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmCystatin_C_Total <- tt_DNAmCystatin_C + DNAmCystatin_C_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmCystatin_C"][1]
    
    # DNAmGDF_15 calculations #
    present_DNAmGDF_15 <- DNAmGDF_15_CpGs %in% colnames(DNAm)
    betas_DNAmGDF_15 <- DNAm[,na.omit(match(DNAmGDF_15_CpGs,colnames(DNAm)))]
    tt_DNAmGDF_15 <- rowSums(sweep(as.matrix(betas_DNAmGDF_15), MARGIN = 2, DNAmGDF_15_Weights[present_DNAmGDF_15], `*`), na.rm = T)
    DNAmGDF_15_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmGDF_15" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmGDF_15_Total <- tt_DNAmGDF_15 + DNAmGDF_15_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmGDF_15"][1]
    
    # DNAmleptin calculations #
    present_DNAmleptin <- DNAmleptin_CpGs %in% colnames(DNAm)
    betas_DNAmleptin <- DNAm[,na.omit(match(DNAmleptin_CpGs,colnames(DNAm)))]
    tt_DNAmleptin <- rowSums(sweep(as.matrix(betas_DNAmleptin), MARGIN = 2, DNAmleptin_Weights[present_DNAmleptin], `*`), na.rm = T)
    SubcomponentResults$DNAmleptin_Total <- tt_DNAmleptin + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmleptin"][1]
    
    # DNAmlog.A1C calculations #
    present_DNAmlog.A1C <- DNAmlog.A1C_CpGs %in% colnames(DNAm)
    betas_DNAmlog.A1C <- DNAm[,na.omit(match(DNAmlog.A1C_CpGs,colnames(DNAm)))]
    tt_DNAmlog.A1C <- rowSums(sweep(as.matrix(betas_DNAmlog.A1C), MARGIN = 2, DNAmlog.A1C_Weights[present_DNAmlog.A1C], `*`), na.rm = T)
    DNAmlog.A1C_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.A1C" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmlog.A1C_Total <- tt_DNAmlog.A1C + DNAmlog.A1C_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.A1C"][1]
    
    # DNAmlog.CRP calculations #
    present_DNAmlog.CRP <- DNAmlog.CRP_CpGs %in% colnames(DNAm)
    betas_DNAmlog.CRP <- DNAm[,na.omit(match(DNAmlog.CRP_CpGs,colnames(DNAm)))]
    tt_DNAmlog.CRP <- rowSums(sweep(as.matrix(betas_DNAmlog.CRP), MARGIN = 2, DNAmlog.CRP_Weights[present_DNAmlog.CRP], `*`), na.rm = T)
    DNAmlog.CRP_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.CRP" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmlog.CRP_Total <- tt_DNAmlog.CRP + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.CRP"][1]
    
    # DNAmPACKYRS calculations #
    present_DNAmPACKYRS <- DNAmPACKYRS_CpGs %in% colnames(DNAm)
    betas_DNAmPACKYRS <- DNAm[,na.omit(match(DNAmPACKYRS_CpGs,colnames(DNAm)))]
    tt_DNAmPACKYRS <- rowSums(sweep(as.matrix(betas_DNAmPACKYRS), MARGIN = 2, DNAmPACKYRS_Weights[present_DNAmPACKYRS], `*`), na.rm = T)
    DNAmPACKYRS_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmPACKYRS" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmPACKYRS_Total <- tt_DNAmPACKYRS + DNAmPACKYRS_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmPACKYRS"][1]
    
    # DNAmpai_1 calculations #
    present_DNAmpai_1 <- DNAmpai_1_CpGs %in% colnames(DNAm)
    betas_DNAmpai_1 <- DNAm[,na.omit(match(DNAmpai_1_CpGs,colnames(DNAm)))]
    tt_DNAmpai_1 <- rowSums(sweep(as.matrix(betas_DNAmpai_1), MARGIN = 2, DNAmpai_1_Weights[present_DNAmpai_1], `*`), na.rm = T)
    DNAmpai_1_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmpai_1" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmpai_1_Total <- tt_DNAmpai_1 + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmpai_1"][1]
    
    # DNAmTIMP_1 calculations #
    present_DNAmTIMP_1 <- DNAmTIMP_1_CpGs %in% colnames(DNAm)
    betas_DNAmTIMP_1 <- DNAm[,na.omit(match(DNAmTIMP_1_CpGs,colnames(DNAm)))]
    tt_DNAmTIMP_1 <- rowSums(sweep(as.matrix(betas_DNAmTIMP_1), MARGIN = 2, DNAmTIMP_1_Weights[present_DNAmTIMP_1], `*`), na.rm = T)
    DNAmTIMP_1_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmTIMP_1" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmTIMP_1_Total <- tt_DNAmTIMP_1 + DNAmTIMP_1_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmTIMP_1"][1]
    
    # Add age and sex information #
    SubcomponentResults$Age_Total <- Ages
    SubcomponentResults$Female_Total <- Sex
    
    # Compute GrimAge2 #
    tt <- rowSums(sweep(as.matrix(SubcomponentResults), MARGIN = 2, GrimAge2_data$glmnet.final1$beta, `*`), na.rm = T)
    tt <- 8.271105 * tt - 61.03936
    
    if(is.null(pheno)){
      tt
    } else{
      pheno$GrimAge2 <- tt
      pheno
    }
    
  } else {
    
    # Print message for imputation #
    message("Imputation of mean CpG Values occured for GrimAge2")
    
    # Initialize empty dataframe for results #
    SubcomponentResults <- data.frame(DNAmGDF_15_Total = rep(NA,nrow(DNAm)),
                                      DNAmB2M_Total = rep(NA,nrow(DNAm)),
                                      DNAmCystatin_C_Total = rep(NA,nrow(DNAm)),
                                      DNAmTIMP_1_Total = rep(NA,nrow(DNAm)),
                                      DNAmadm_Total = rep(NA,nrow(DNAm)),
                                      DNAmpai_1_Total = rep(NA,nrow(DNAm)),
                                      DNAmleptin_Total = rep(NA,nrow(DNAm)),
                                      DNAmPACKYRS_Total = rep(NA,nrow(DNAm)),
                                      DNAmlog.CRP_Total = rep(NA,nrow(DNAm)),
                                      DNAmlog.A1C_Total = rep(NA,nrow(DNAm)),
                                      Age_Total = rep(NA,nrow(DNAm)),
                                      Female_Total = rep(NA,nrow(DNAm)))
    
    ### Impute missing DNAmadm CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmadm <- DNAmadm_CpGs[!(DNAmadm_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmadm <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmadm))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmadm)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmadm[j],names(CpGImputation))]
      tempDNAm_DNAmadm[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmadm) <- missingCpGs_DNAmadm
    # Add missing CpGs to DNAm #
    DNAm_DNAmadm <- cbind(DNAm,tempDNAm_DNAmadm)
    # Extract out needed CpGs from DNAm #
    betas_DNAmadm <- DNAm_DNAmadm[,match(DNAmadm_CpGs,colnames(DNAm_DNAmadm))]
    # Compute subcomponent scores #
    tt_DNAmadm <- rowSums(sweep(betas_DNAmadm, MARGIN = 2, DNAmadm_Weights, `*`), na.rm = TRUE)
    DNAmadm_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmadm" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmadm_Total <- tt_DNAmadm + DNAmadm_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmadm"][1]
    
    ### Impute missing DNAmB2M CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmB2M <- DNAmB2M_CpGs[!(DNAmB2M_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmB2M <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmB2M))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmB2M)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmB2M[j],names(CpGImputation))]
      tempDNAm_DNAmB2M[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmB2M) <- missingCpGs_DNAmB2M
    # Add missing CpGs to DNAm #
    DNAm_DNAmB2M <- cbind(DNAm,tempDNAm_DNAmB2M)
    # Extract out needed CpGs from DNAm #
    betas_DNAmB2M <- DNAm_DNAmB2M[,match(DNAmB2M_CpGs,colnames(DNAm_DNAmB2M))]
    # Compute subcomponent scores #
    tt_DNAmB2M <- rowSums(sweep(betas_DNAmB2M, MARGIN = 2, DNAmB2M_Weights, `*`), na.rm = TRUE)
    DNAmB2M_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmB2M" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmB2M_Total <- tt_DNAmB2M + DNAmB2M_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmB2M"][1]
    
    ### Impute missing DNAmCystatin_C CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmCystatin_C <- DNAmCystatin_C_CpGs[!(DNAmCystatin_C_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmCystatin_C <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmCystatin_C))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmCystatin_C)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmCystatin_C[j],names(CpGImputation))]
      tempDNAm_DNAmCystatin_C[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmCystatin_C) <- missingCpGs_DNAmCystatin_C
    # Add missing CpGs to DNAm #
    DNAm_DNAmCystatin_C <- cbind(DNAm,tempDNAm_DNAmCystatin_C)
    # Extract out needed CpGs from DNAm #
    betas_DNAmCystatin_C <- DNAm_DNAmCystatin_C[,match(DNAmCystatin_C_CpGs,colnames(DNAm_DNAmCystatin_C))]
    # Compute subcomponent scores #
    tt_DNAmCystatin_C <- rowSums(sweep(betas_DNAmCystatin_C, MARGIN = 2, DNAmCystatin_C_Weights, `*`), na.rm = TRUE)
    DNAmCystatin_C_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmCystatin_C" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmCystatin_C_Total <- tt_DNAmCystatin_C + DNAmCystatin_C_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmCystatin_C"][1]
    
    ### Impute missing DNAmGDF_15 CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmGDF_15 <- DNAmGDF_15_CpGs[!(DNAmGDF_15_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmGDF_15 <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmGDF_15))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmGDF_15)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmGDF_15[j],names(CpGImputation))]
      tempDNAm_DNAmGDF_15[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmGDF_15) <- missingCpGs_DNAmGDF_15
    # Add missing CpGs to DNAm #
    DNAm_DNAmGDF_15 <- cbind(DNAm,tempDNAm_DNAmGDF_15)
    # Extract out needed CpGs from DNAm #
    betas_DNAmGDF_15 <- DNAm_DNAmGDF_15[,match(DNAmGDF_15_CpGs,colnames(DNAm_DNAmGDF_15))]
    # Compute subcomponent scores #
    tt_DNAmGDF_15 <- rowSums(sweep(betas_DNAmGDF_15, MARGIN = 2, DNAmGDF_15_Weights, `*`), na.rm = TRUE)
    DNAmGDF_15_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmGDF_15" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmGDF_15_Total <- tt_DNAmGDF_15 + DNAmGDF_15_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmGDF_15"][1]
    
    ### Impute missing DNAmleptin CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmleptin <- DNAmleptin_CpGs[!(DNAmleptin_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmleptin <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmleptin))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmleptin)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmleptin[j],names(CpGImputation))]
      tempDNAm_DNAmleptin[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmleptin) <- missingCpGs_DNAmleptin
    # Add missing CpGs to DNAm #
    DNAm_DNAmleptin <- cbind(DNAm,tempDNAm_DNAmleptin)
    # Extract out needed CpGs from DNAm #
    betas_DNAmleptin <- DNAm_DNAmleptin[,match(DNAmleptin_CpGs,colnames(DNAm_DNAmleptin))]
    # Compute subcomponent scores #
    tt_DNAmleptin <- rowSums(sweep(betas_DNAmleptin, MARGIN = 2, DNAmleptin_Weights, `*`), na.rm = TRUE)
    SubcomponentResults$DNAmleptin_Total <- tt_DNAmleptin + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmleptin"][1]
    
    ### Impute missing DNAmlog.A1C CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmlog.A1C <- DNAmlog.A1C_CpGs[!(DNAmlog.A1C_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmlog.A1C <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmlog.A1C))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmlog.A1C)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmlog.A1C[j],names(CpGImputation))]
      tempDNAm_DNAmlog.A1C[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmlog.A1C) <- missingCpGs_DNAmlog.A1C
    # Add missing CpGs to DNAm #
    DNAm_DNAmlog.A1C <- cbind(DNAm,tempDNAm_DNAmlog.A1C)
    # Extract out needed CpGs from DNAm #
    betas_DNAmlog.A1C <- DNAm_DNAmlog.A1C[,match(DNAmlog.A1C_CpGs,colnames(DNAm_DNAmlog.A1C))]
    # Compute subcomponent scores #
    tt_DNAmlog.A1C <- rowSums(sweep(betas_DNAmlog.A1C, MARGIN = 2, DNAmlog.A1C_Weights, `*`), na.rm = TRUE)
    DNAmlog.A1C_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.A1C" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmlog.A1C_Total <- tt_DNAmlog.A1C + DNAmlog.A1C_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.A1C"][1]
    
    ### Impute missing DNAmlog.CRP CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmlog.CRP <- DNAmlog.CRP_CpGs[!(DNAmlog.CRP_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmlog.CRP <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmlog.CRP))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmlog.CRP)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmlog.CRP[j],names(CpGImputation))]
      tempDNAm_DNAmlog.CRP[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmlog.CRP) <- missingCpGs_DNAmlog.CRP
    # Add missing CpGs to DNAm #
    DNAm_DNAmlog.CRP <- cbind(DNAm,tempDNAm_DNAmlog.CRP)
    # Extract out needed CpGs from DNAm #
    betas_DNAmlog.CRP <- DNAm_DNAmlog.CRP[,match(DNAmlog.CRP_CpGs,colnames(DNAm_DNAmlog.CRP))]
    # Compute subcomponent scores #
    tt_DNAmlog.CRP <- rowSums(sweep(betas_DNAmlog.CRP, MARGIN = 2, DNAmlog.CRP_Weights, `*`), na.rm = TRUE)
    SubcomponentResults$DNAmlog.CRP_Total <- tt_DNAmlog.CRP + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmlog.CRP"][1]
    
    ### Impute missing DNAmPACKYRS CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmPACKYRS <- DNAmPACKYRS_CpGs[!(DNAmPACKYRS_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmPACKYRS <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmPACKYRS))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmPACKYRS)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmPACKYRS[j],names(CpGImputation))]
      tempDNAm_DNAmPACKYRS[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmPACKYRS) <- missingCpGs_DNAmPACKYRS
    # Add missing CpGs to DNAm #
    DNAm_DNAmPACKYRS <- cbind(DNAm,tempDNAm_DNAmPACKYRS)
    # Extract out needed CpGs from DNAm #
    betas_DNAmPACKYRS <- DNAm_DNAmPACKYRS[,match(DNAmPACKYRS_CpGs,colnames(DNAm_DNAmPACKYRS))]
    # Compute subcomponent scores #
    tt_DNAmPACKYRS <- rowSums(sweep(betas_DNAmPACKYRS, MARGIN = 2, DNAmPACKYRS_Weights, `*`), na.rm = TRUE)
    DNAmPACKYRS_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmPACKYRS" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmPACKYRS_Total <- tt_DNAmPACKYRS + DNAmPACKYRS_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmPACKYRS"][1]
    
    ### Impute missing DNAmpai_1 CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmpai_1 <- DNAmpai_1_CpGs[!(DNAmpai_1_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmpai_1 <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmpai_1))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmpai_1)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmpai_1[j],names(CpGImputation))]
      tempDNAm_DNAmpai_1[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmpai_1) <- missingCpGs_DNAmpai_1
    # Add missing CpGs to DNAm #
    DNAm_DNAmpai_1 <- cbind(DNAm,tempDNAm_DNAmpai_1)
    # Extract out needed CpGs from DNAm #
    betas_DNAmpai_1 <- DNAm_DNAmpai_1[,match(DNAmpai_1_CpGs,colnames(DNAm_DNAmpai_1))]
    # Compute subcomponent scores #
    tt_DNAmpai_1 <- rowSums(sweep(betas_DNAmpai_1, MARGIN = 2, DNAmpai_1_Weights, `*`), na.rm = TRUE)
    SubcomponentResults$DNAmpai_1_Total <- tt_DNAmpai_1 + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmpai_1"][1]
    
    ### Impute missing DNAmTIMP_1 CpGs ###
    # Define missing CpGs #
    missingCpGs_DNAmTIMP_1 <- DNAmTIMP_1_CpGs[!(DNAmTIMP_1_CpGs %in% colnames(DNAm))]
    # Create temporary matrix to house missing CpGs #
    tempDNAm_DNAmTIMP_1 <- matrix(nrow = dim(DNAm)[1], ncol = length(missingCpGs_DNAmTIMP_1))
    # Input GoldenStandard mean for each missing CpG #
    for(j in 1:length(missingCpGs_DNAmTIMP_1)){
      meanVals <- CpGImputation[match(missingCpGs_DNAmTIMP_1[j],names(CpGImputation))]
      tempDNAm_DNAmTIMP_1[,j] <- rep(meanVals,dim(DNAm)[1])
    }
    # Change temp CpG column names #
    colnames(tempDNAm_DNAmTIMP_1) <- missingCpGs_DNAmTIMP_1
    # Add missing CpGs to DNAm #
    DNAm_DNAmTIMP_1 <- cbind(DNAm,tempDNAm_DNAmTIMP_1)
    # Extract out needed CpGs from DNAm #
    betas_DNAmTIMP_1 <- DNAm_DNAmTIMP_1[,match(DNAmTIMP_1_CpGs,colnames(DNAm_DNAmTIMP_1))]
    # Compute subcomponent scores #
    tt_DNAmTIMP_1 <- rowSums(sweep(betas_DNAmTIMP_1, MARGIN = 2, DNAmTIMP_1_Weights, `*`), na.rm = TRUE)
    DNAmTIMP_1_ages <- GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmTIMP_1" & GrimAge2_data$dnam.all$var == "Age"] * Ages
    SubcomponentResults$DNAmTIMP_1_Total <- tt_DNAmTIMP_1 + DNAmTIMP_1_ages + GrimAge2_data$dnam.all$beta[GrimAge2_data$dnam.all$Y.pred == "DNAmTIMP_1"][1]
    
    # Add age and sex information #
    SubcomponentResults$Age_Total <- Ages
    SubcomponentResults$Female_Total <- Sex
    
    # Compute GrimAge2 #
    tt <- rowSums(sweep(as.matrix(SubcomponentResults), MARGIN = 2, GrimAge2_data$glmnet.final1$beta, `*`), na.rm = T)
    tt <- 8.271105 * tt - 61.03936
    
    if(is.null(pheno)){
      tt
    } else{
      pheno$GrimAge2 <- tt
      pheno
    }
    
  }
  
}


#######################################################
### Principal Component Clocks Function ### ----
#######################################################

calcPCClocksEdit <- function(path_to_PCClocks_directory, datMeth, datPheno, GoldenStandard){
  #path_to_PCClocks should end with a "/"
  #datMeth is a matrix of methylation Beta values, where row names are samples, and 
  #   column names are CpGs
  #datPheno has rows as samples and columns as phenotype variables. This can also
  #   include the original clocks if you used the Horvath online calculator as well.
  #   It MUST include a column named "Age" and a column named "Female"
  
  if(!("Age" %in% variable.names(datPheno))){
    stop("Error: datPheno must have a column named Age")
  }
  if(!("Female" %in% variable.names(datPheno))){
    stop("Error: datPheno must have a column named Female")
  }
  if(sum(startsWith(colnames(datMeth),"cg")) == 0){
    warning("Warning: It looks like you may need to format datMeth using t(datMeth) to get samples as rows!")
  }
  
  
  #Note: this code assumes all your files are in one working directory. Alter the code as needed based on file locations.
  #Load packages
  pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }
  pkgTest("dplyr")
  library(dplyr)
  pkgTest("tibble")
  library(tibble)
  pkgTest("tidyr")
  library(tidyr)
  
  
  #In datPheno, rows are samples and columns are phenotypic variables. 
  #One of the phenotypic variables must be "Age", and another one "Female" (coded as Female = 1, Male = 0; should be a numeric variable as this will be included in PCGrimAge calculation)
  #Also ensure that the order of datMeth sample IDs matches your phenotype data sample IDs, otherwise your data will be scrambled
  
  load(file = paste(path_to_PCClocks_directory,"CalcAllPCClocks.RData", sep = ""))
  
  message("PCClocks Data successfully loaded")
  
  #If needed: Fill in missing CpGs needed for calculation of PCs; use mean values from GSE40279 (Hannum 2013; blood)- note that for other tissues you might prefer to use a different one
  datMeth <- as.data.frame(datMeth)
  if(length(c(CpGs[!(CpGs %in% colnames(datMeth))],CpGs[apply(datMeth[,colnames(datMeth) %in% CpGs], 2, function(x)all(is.na(x)))])) == 0){
    message("No CpGs were NA for all samples")
  } else{
    missingCpGs <- c(CpGs[!(CpGs %in% colnames(datMeth))])
    datMeth[,missingCpGs] <- NA
    datMeth = datMeth[,CpGs]
    missingCpGs <- CpGs[apply(datMeth[,CpGs], 2, function(x)all(is.na(x)))]
    for(i in 1:length(missingCpGs)){
      datMeth[,missingCpGs[i]] <- GoldenStandard[missingCpGs[i]] # This needs to be edited to new GS means #
    }
    message("Any missing CpGs successfully filled in (see function for more details)")
  }
  
  # If any probes are not in the golden standard dataset supplied, set them all equal to 0 #
  
  
  # Prepare methylation data for calculation of PC Clocks (subset to 78,464 CpGs and perform imputation if needed)
  datMeth <- datMeth[,CpGs]
  meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
  datMeth <- apply(datMeth,2,meanimpute)
  # Note: you may substitute another imputation method of your choice (e.g. KNN), but we have not found the method makes a significant difference.
  # If any probes are not in the golden standard dataset supplied, set them all equal to 0 #
  datMeth[is.na(datMeth)] <- 0
  message("Mean imputation successfully completed for any missing CpG values")
  
  #Initialize a data frame for PC clocks
  DNAmAge <- datPheno
  
  #var = readline(prompt = "To check whether datMeth and datPheno match up, type the column name in datPheno with sample names (or type skip):")
  var = "skip"
  if(var != "skip"){
    if(sum(DNAmAge[,var] == rownames(datMeth)) != dim(DNAmAge[,var])[1]){
      warning("Warning: It would appear that datPheno and datMeth do not have matching sample order! Check your inputs!")
    } else message("datPheno and datMeth sample order verified to match!")
  }
  
  message("Calculating PC Clocks now")
  
  #Calculate PC Clocks
  DNAmAge$PCHorvath1 <- as.numeric(anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath1$center) %*% CalcPCHorvath1$rotation %*% CalcPCHorvath1$model + CalcPCHorvath1$intercept))
  DNAmAge$PCHorvath2 <- as.numeric(anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath2$center) %*% CalcPCHorvath2$rotation %*% CalcPCHorvath2$model + CalcPCHorvath2$intercept))
  DNAmAge$PCHannum <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCHannum$center) %*% CalcPCHannum$rotation %*% CalcPCHannum$model + CalcPCHannum$intercept)
  DNAmAge$PCPhenoAge <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCPhenoAge$center) %*% CalcPCPhenoAge$rotation %*% CalcPCPhenoAge$model + CalcPCPhenoAge$intercept)
  DNAmAge$PCDNAmTL <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCDNAmTL$center) %*% CalcPCDNAmTL$rotation %*% CalcPCDNAmTL$model + CalcPCDNAmTL$intercept)
  temp <- cbind(sweep(as.matrix(datMeth),2,CalcPCGrimAge$center) %*% CalcPCGrimAge$rotation,Female = DNAmAge$Female,Age = DNAmAge$Age)
  DNAmAge$PCPACKYRS <- as.numeric(temp[,names(CalcPCGrimAge$PCPACKYRS.model)] %*% CalcPCGrimAge$PCPACKYRS.model + CalcPCGrimAge$PCPACKYRS.intercept)
  DNAmAge$PCADM <- as.numeric(temp[,names(CalcPCGrimAge$PCADM.model)] %*% CalcPCGrimAge$PCADM.model + CalcPCGrimAge$PCADM.intercept)
  DNAmAge$PCB2M <- as.numeric(temp[,names(CalcPCGrimAge$PCB2M.model)] %*% CalcPCGrimAge$PCB2M.model + CalcPCGrimAge$PCB2M.intercept)
  DNAmAge$PCCystatinC <- as.numeric(temp[,names(CalcPCGrimAge$PCCystatinC.model)] %*% CalcPCGrimAge$PCCystatinC.model + CalcPCGrimAge$PCCystatinC.intercept)
  DNAmAge$PCGDF15 <- as.numeric(temp[,names(CalcPCGrimAge$PCGDF15.model)] %*% CalcPCGrimAge$PCGDF15.model + CalcPCGrimAge$PCGDF15.intercept)
  DNAmAge$PCLeptin <- as.numeric(temp[,names(CalcPCGrimAge$PCLeptin.model)] %*% CalcPCGrimAge$PCLeptin.model + CalcPCGrimAge$PCLeptin.intercept)
  DNAmAge$PCPAI1 <- as.numeric(temp[,names(CalcPCGrimAge$PCPAI1.model)] %*% CalcPCGrimAge$PCPAI1.model + CalcPCGrimAge$PCPAI1.intercept)
  DNAmAge$PCTIMP1 <- as.numeric(temp[,names(CalcPCGrimAge$PCTIMP1.model)] %*% CalcPCGrimAge$PCTIMP1.model + CalcPCGrimAge$PCTIMP1.intercept)
  DNAmAge$PCGrimAge <- as.numeric(as.matrix(DNAmAge[,CalcPCGrimAge$components]) %*% CalcPCGrimAge$PCGrimAge.model + CalcPCGrimAge$PCGrimAge.intercept)
  rm(CalcPCHorvath1,CalcPCHorvath2,CalcPCHannum,CalcPCPhenoAge,CalcPCDNAmTL,CalcPCGrimAge,temp,imputeMissingCpGs)
  
  message("PC Clocks successfully calculated!")
  
  return(DNAmAge)
}