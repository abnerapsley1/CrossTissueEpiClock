### Cross-Tissue EpiClock Statistical Analyses ### ----

### Loading Packages ### ----
library(knitr)
library(dplyr)
library(datapasta)
library(ggplot2)
library(psych)
library(corrplot)
library(tibble)
library(ggvenn)
library('org.Hs.eg.db')
library(pathview)
library(clusterProfiler)
library(glmnet)
library(caret)
library(glmnetUtils)
library(umx)
library(TOAST) 
library(matlib)
library(ggpubr)
library(lme4)
library(lmerTest)
library(impute)
options(digits=5)
options(device = "quartz")

### Read in Data ### ----
# Online Horvath Calculator Results #
data_horv_b1 <- read.csv("Horvath_ClockProjResults_Batch1.csv")
data_horv_b2 <- read.csv("Horvath_ClockProjResults_Batch2.csv")
data_horv_b3 <- read.csv("Horvath_ClockProjResults_Batch3.csv")
data_horv_b4 <- read.csv("Horvath_ClockProjResults_Batch4.csv")
# MethyAge (PACE) Results #
data_PACE_b1 <- read.csv("MethyAge_EpiClocks_Batch1.csv")[,-1]
data_PACE_b2 <- read.csv("MethyAge_EpiClocks_Batch2.csv")[,-1]
data_PACE_b3 <- read.csv("MethyAge_EpiClocks_Batch3.csv")[,-1]
data_PACE_b4 <- read.csv("MethyAge_EpiClocks_Batch4.csv")[,-1]
# PC Clocks Results #
data_PC_b1 <- read.csv("PCClocks_Batch1.csv")[,-1]
data_PC_b2 <- read.csv("PCClocks_Batch2.csv")[,-1]
data_PC_b3 <- read.csv("PCClocks_Batch3.csv")[,-1]
data_PC_b4 <- read.csv("PCClocks_Batch4.csv")[,-1]
# PedBE Clock Results #
data_PBE_b1 <- read.csv("PedBe_age_Batch1.csv")[,-1]
data_PBE_b2 <- read.csv("PedBe_age_Batch2.csv")[,-1]
data_PBE_b3 <- read.csv("PedBe_age_Batch3.csv")[,-1]
data_PBE_b4 <- read.csv("PedBe_age_Batch4.csv")[,-1]
# DunedinPoAm Clock Results #
data_PoAm_b1 <- read.csv("PoAm_Clock_Batch1.csv")[,]
data_PoAm_b2 <- read.csv("PoAm_Clock_Batch2.csv")[,]
data_PoAm_b3 <- read.csv("PoAm_Clock_Batch3.csv")[,]
data_PoAm_b4 <- read.csv("PoAm_Clock_Batch4.csv")[,]
# U01 Phenotypes #
data_pheno_TRN <- read.csv("TL_Study_Data_2.28.2020.csv")
# CHS Phenotypes #
data_pheno_CHS  <- read.csv("CHS_Phenotypes.csv")
# Manual Clock Computation Results #
data_mcc <- read.csv("ManualClockComputationData.csv")[,-c(1,3:5)]


### Clean up Data ### ----
## data_horv ##
# Add Batch variable #
data_horv_b1$Batch <- "1"
data_horv_b2$Batch <- "2"
data_horv_b3$Batch <- "3"
data_horv_b4$Batch <- "4"
# Ensure columns are in same order #
identical(names(data_horv_b1), names(data_horv_b2))
identical(names(data_horv_b1), names(data_horv_b3))
identical(names(data_horv_b2), names(data_horv_b3))
identical(names(data_horv_b4), names(data_horv_b1))
# Combine data frames #
data_horv <- rbind(data_horv_b1, data_horv_b2, data_horv_b3, data_horv_b4)
# Remove X from ID in data_horv #
data_horv$SID <- sub('.', '', data_horv$SID)
# Rename ID column #
colnames(data_horv)[1] <- "SampleID"
# Add Study_Cohort variable #
data_horv$Study_Cohort <- NA
for (i in 1:nrow(data_horv)){
  if (nchar(data_horv[i,1]) > 3) {data_horv$Study_Cohort[i] = "CHS"} else {data_horv$Study_Cohort[i] = "Telomere"}
}
# Remove unnecessary columns #
data_horv <- data_horv[,-c(2:5,9:35,42:91,101,131,146,151:152,162,164)]
# Ensure variables are correct format #
data_horv$Tissue <- as.factor(data_horv$Tissue)
## data_PACE ##
# Ensure columns are in same order #
identical(names(data_PACE_b1), names(data_PACE_b2))
identical(names(data_PACE_b1), names(data_PACE_b3))
identical(names(data_PACE_b2), names(data_PACE_b3))
identical(names(data_PACE_b1), names(data_PACE_b4))
# Combine data frames #
data_PACE <- rbind(data_PACE_b1, data_PACE_b2, data_PACE_b3, data_PACE_b4)
# Rename variables #
colnames(data_PACE)[2:5] <- c("mAge_Horvath_MA","mAge_Hunnum_MA","PhenoAge_MA","PACE_MA")
## data_PC ##
# Ensure columns are in same order #
identical(names(data_PC_b1), names(data_PC_b2))
identical(names(data_PC_b1), names(data_PC_b3))
identical(names(data_PC_b2), names(data_PC_b3))
identical(names(data_PC_b1), names(data_PC_b4))
data_PC_b1 <- data_PC_b1[,-c(12:15)]
data_PC_b2 <- data_PC_b2[,-c(12:15)]
data_PC_b3 <- data_PC_b3[,-c(12:15)]
# Combine data frames #
data_PC <- rbind(data_PC_b1, data_PC_b2, data_PC_b3, data_PC_b4)
# Remove unnecessary columns #
data_PC <- data_PC[,-c(2:10,12:14)]
# Rename variables #
colnames(data_PC)[2] <- "Individual"
## data_PBE ##
# Combine data frames #
data_PBE <- rbind(data_PBE_b1, data_PBE_b2, data_PBE_b3, data_PBE_b4)
# Rename variables #
colnames(data_PBE)[1] <- "SampleID"
colnames(data_PBE)[2] <- "DNAmAgePedBE"
## data_PoAm ##
# Combine data frames #
data_PoAm <- rbind(data_PoAm_b1, data_PoAm_b2, data_PoAm_b3, data_PoAm_b4)
# Rename variables #
colnames(data_PoAm)[1] <- "SampleID"
colnames(data_PoAm)[2] <- "PoAm"
## Combine all dataframes ##
data_horv_pace <- merge(data_horv, data_PACE, by = "SampleID")
data_horv_pace_PC <- merge(data_horv_pace, data_PC, by = "SampleID")
data_horv_pace_PC_PBE <- merge(data_horv_pace_PC, data_PBE, by = "SampleID")
data <- merge(data_horv_pace_PC_PBE, data_PoAm, by = "SampleID")

## Remove Problem Samples ##
removeSamples <- data$SampleID[c(57,77,213,228)]
data <- data[-c(57,77,213,228),]


## Add Manual Clock Computations ##
# Ensure same SampleIDs #
table(data$SampleID %in% data_mcc$SampleID)
# Change data SampleID name for one individual #
data$SampleID[116] <- "2029215"
# Combine data and manually computed clocks #
data <- merge(data, data_mcc, by = "SampleID")

## Compute Age Accelerations for Manually Computed Clocks ##
## Tissue-Agnostic Intercept Included ##
# Horvath1 #
AccelModelHorvath1 <- lm(mcc_Horvath1 ~ 1 + Age, data)
data$mcc_Horvath1_TissAg_Int <- AccelModelHorvath1$residuals
# Hannum #
AccelModelHannum <- lm(mcc_Hannum ~ 1 + Age, data)
data$mcc_Hannum_TissAg_Int <- AccelModelHannum$residuals
# PhenoAge #
AccelModelPhenoAge <- lm(mcc_PhenoAge ~ 1 + Age, data)
data$mcc_PhenoAge_TissAg_Int <- AccelModelPhenoAge$residuals
# GrimAge2 #
AccelModelGrimAge2 <- lm(mcc_GrimAge2 ~ 1 + Age, data)
data$mcc_GrimAge2_TissAg_Int <- AccelModelGrimAge2$residuals
# Horvath2 #
AccelModelHorvath2 <- lm(mcc_Horvath2 ~ 1 + Age, data)
data$mcc_Horvath2_TissAg_Int <- AccelModelHorvath2$residuals
# PedBE #
AccelModelPedBE <- lm(mcc_PedBE ~ 1 + Age, data)
data$mcc_PedBE_TissAg_Int <- AccelModelPedBE$residuals

## Tissue-Specific Intercept Included ##
# Horvath1 #
AccelModelHorvath1_TS_Buccal <- lm(mcc_Horvath1 ~ 1 + Age, data[data$Tissue == "Buccal",])
AccelModelHorvath1_TS_Saliva <- lm(mcc_Horvath1 ~ 1 + Age, data[data$Tissue == "Saliva",])
AccelModelHorvath1_TS_DBS <- lm(mcc_Horvath1 ~ 1 + Age, data[data$Tissue == "DBS",])
AccelModelHorvath1_TS_BuffyCoat <- lm(mcc_Horvath1 ~ 1 + Age, data[data$Tissue == "Buffy Coat",])
AccelModelHorvath1_TS_PBMC <- lm(mcc_Horvath1 ~ 1 + Age, data[data$Tissue == "PBMC",])

data_TS <- data.frame(SampleID = c(data$SampleID[data$Tissue == "Buccal"],
                                            data$SampleID[data$Tissue == "Saliva"],
                                            data$SampleID[data$Tissue == "DBS"],
                                            data$SampleID[data$Tissue == "Buffy Coat"],
                                            data$SampleID[data$Tissue == "PBMC"]),
                               mcc_Horvath1_TissS_Int = c(AccelModelHorvath1_TS_Buccal$residuals,
                                                          AccelModelHorvath1_TS_Saliva$residuals,
                                                          AccelModelHorvath1_TS_DBS$residuals,
                                                          AccelModelHorvath1_TS_BuffyCoat$residuals,
                                                          AccelModelHorvath1_TS_PBMC$residuals))

# Hannum #
AccelModelHannum_TS_Buccal <- lm(mcc_Hannum ~ 1 + Age, data[data$Tissue == "Buccal",])
AccelModelHannum_TS_Saliva <- lm(mcc_Hannum ~ 1 + Age, data[data$Tissue == "Saliva",])
AccelModelHannum_TS_DBS <- lm(mcc_Hannum ~ 1 + Age, data[data$Tissue == "DBS",])
AccelModelHannum_TS_BuffyCoat <- lm(mcc_Hannum ~ 1 + Age, data[data$Tissue == "Buffy Coat",])
AccelModelHannum_TS_PBMC <- lm(mcc_Hannum ~ 1 + Age, data[data$Tissue == "PBMC",])

data_TS$mccHannum_TissS_Int <- c(AccelModelHannum_TS_Buccal$residuals,
                                 AccelModelHannum_TS_Saliva$residuals,
                                 AccelModelHannum_TS_DBS$residuals,
                                 AccelModelHannum_TS_BuffyCoat$residuals,
                                 AccelModelHannum_TS_PBMC$residuals)

# PhenoAge #
AccelModelPhenoAge_TS_Buccal <- lm(mcc_PhenoAge ~ 1 + Age, data[data$Tissue == "Buccal",])
AccelModelPhenoAge_TS_Saliva <- lm(mcc_PhenoAge ~ 1 + Age, data[data$Tissue == "Saliva",])
AccelModelPhenoAge_TS_DBS <- lm(mcc_PhenoAge ~ 1 + Age, data[data$Tissue == "DBS",])
AccelModelPhenoAge_TS_BuffyCoat <- lm(mcc_PhenoAge ~ 1 + Age, data[data$Tissue == "Buffy Coat",])
AccelModelPhenoAge_TS_PBMC <- lm(mcc_PhenoAge ~ 1 + Age, data[data$Tissue == "PBMC",])

data_TS$mccPhenoAge_TissS_Int <- c(AccelModelPhenoAge_TS_Buccal$residuals,
                                 AccelModelPhenoAge_TS_Saliva$residuals,
                                 AccelModelPhenoAge_TS_DBS$residuals,
                                 AccelModelPhenoAge_TS_BuffyCoat$residuals,
                                 AccelModelPhenoAge_TS_PBMC$residuals)

# GrimAge2 #
AccelModelGrimAge2_TS_Buccal <- lm(mcc_GrimAge2 ~ 1 + Age, data[data$Tissue == "Buccal",])
AccelModelGrimAge2_TS_Saliva <- lm(mcc_GrimAge2 ~ 1 + Age, data[data$Tissue == "Saliva",])
AccelModelGrimAge2_TS_DBS <- lm(mcc_GrimAge2 ~ 1 + Age, data[data$Tissue == "DBS",])
AccelModelGrimAge2_TS_BuffyCoat <- lm(mcc_GrimAge2 ~ 1 + Age, data[data$Tissue == "Buffy Coat",])
AccelModelGrimAge2_TS_PBMC <- lm(mcc_GrimAge2 ~ 1 + Age, data[data$Tissue == "PBMC",])

data_TS$mccGrimAge2_TissS_Int <- c(AccelModelGrimAge2_TS_Buccal$residuals,
                                 AccelModelGrimAge2_TS_Saliva$residuals,
                                 AccelModelGrimAge2_TS_DBS$residuals,
                                 AccelModelGrimAge2_TS_BuffyCoat$residuals,
                                 AccelModelGrimAge2_TS_PBMC$residuals)

# Horvath2 #
AccelModelHorvath2_TS_Buccal <- lm(mcc_Horvath2 ~ 1 + Age, data[data$Tissue == "Buccal",])
AccelModelHorvath2_TS_Saliva <- lm(mcc_Horvath2 ~ 1 + Age, data[data$Tissue == "Saliva",])
AccelModelHorvath2_TS_DBS <- lm(mcc_Horvath2 ~ 1 + Age, data[data$Tissue == "DBS",])
AccelModelHorvath2_TS_BuffyCoat <- lm(mcc_Horvath2 ~ 1 + Age, data[data$Tissue == "Buffy Coat",])
AccelModelHorvath2_TS_PBMC <- lm(mcc_Horvath2 ~ 1 + Age, data[data$Tissue == "PBMC",])

data_TS$mccHorvath2_TissS_Int <- c(AccelModelHorvath2_TS_Buccal$residuals,
                                 AccelModelHorvath2_TS_Saliva$residuals,
                                 AccelModelHorvath2_TS_DBS$residuals,
                                 AccelModelHorvath2_TS_BuffyCoat$residuals,
                                 AccelModelHorvath2_TS_PBMC$residuals)

# PedBE #
AccelModelPedBE_TS_Buccal <- lm(mcc_PedBE ~ 1 + Age, data[data$Tissue == "Buccal",])
AccelModelPedBE_TS_Saliva <- lm(mcc_PedBE ~ 1 + Age, data[data$Tissue == "Saliva",])
AccelModelPedBE_TS_DBS <- lm(mcc_PedBE ~ 1 + Age, data[data$Tissue == "DBS",])
AccelModelPedBE_TS_BuffyCoat <- lm(mcc_PedBE ~ 1 + Age, data[data$Tissue == "Buffy Coat",])
AccelModelPedBE_TS_PBMC <- lm(mcc_PedBE ~ 1 + Age, data[data$Tissue == "PBMC",])

data_TS$mccPedBE_TissS_Int <- c(AccelModelPedBE_TS_Buccal$residuals,
                                 AccelModelPedBE_TS_Saliva$residuals,
                                 AccelModelPedBE_TS_DBS$residuals,
                                 AccelModelPedBE_TS_BuffyCoat$residuals,
                                 AccelModelPedBE_TS_PBMC$residuals)

data <- merge(data, data_TS, by = "SampleID")

## Compute PCs for Control Probes ##
ControlProbes_Batch1 <- readRDS("Batch1_ControlBeta.RDS")
ControlProbes_Batch2 <- readRDS("Batch2_ControlBeta.RDS")
ControlProbes_Batch3 <- readRDS("Batch3_ControlBeta.RDS")
ControlProbes_Batch4 <- readRDS("Batch4_ControlBeta.RDS")
ControlProbes_all <- cbind(ControlProbes_Batch1,
                           ControlProbes_Batch2,
                           ControlProbes_Batch3,
                           ControlProbes_Batch4)
# Remove bad samples #
ControlProbes_all <- t(ControlProbes_all[,!(colnames(ControlProbes_all) %in% removeSamples)])
# Perform PCA #
PCA_Results <- prcomp(ControlProbes_all)
PC_ControlProbes <- PCA_Results[["x"]]
Summary_PCA <- summary(PCA_Results)
Summary_PCA
# Add first 30 PCs to data #
PC_ControlProbes <- tibble::rownames_to_column(as.data.frame(PC_ControlProbes), "SampleID")
data <- merge(data, PC_ControlProbes[,1:31])

## Add Batch and Cell Composition Residualized Clocks ##
## Buccal (K=3) ##
#Bval_b1 <- readRDS("bVals_SampFilter_BMIQNorm_Batch1.RDS")
#Bval_b1 <- as.data.frame(Bval_b1)
#Bval_b1 <- tibble::rownames_to_column(Bval_b1, "ProbeID")
#Bval_b2 <- readRDS("bVals_SampFilter_BMIQNorm_Batch2.RDS")
#Bval_b2 <- as.data.frame(Bval_b2)
#Bval_b2 <- tibble::rownames_to_column(Bval_b2, "ProbeID")
#Bval_b3 <- readRDS("bVals_SampFilter_BMIQNorm_Batch3.RDS")
#Bval_b3 <- as.data.frame(Bval_b3)
#Bval_b3 <- tibble::rownames_to_column(Bval_b3, "ProbeID")
#Bval_b4 <- readRDS("bVals_SampFilter_BMIQNorm_Batch4.RDS")
#Bval_b4 <- as.data.frame(Bval_b4)
#Bval_b4 <- tibble::rownames_to_column(Bval_b4, "ProbeID")
# Combine all matrices #
#Bval_b1b2 <- merge(Bval_b1, Bval_b2, by = "ProbeID")
#Bval_b1b2b3 <- merge(Bval_b1b2, Bval_b3, by = "ProbeID")
#Bval <- merge(Bval_b1b2b3, Bval_b4, by = "ProbeID")
#Bval_b1 <- NULL
#Bval_b2 <- NULL
#Bval_b3 <- NULL
#Bval_b4 <- NULL
#Bval_b1b2 <- NULL
#Bval_b1b2b3 <- NULL
# Split Bval by tissue type for Buccal cells #
#keep <- colnames(Bval) %in% dplyr::filter(data, data$Tissue == "Buccal")$SampleID
#keep[1] <- TRUE
#Bval_Buccal <- Bval[,keep]
# Estimate cell proportions with reference-free method #
#RefFree_Deconvolution <- myRefFreeCellMix(as.matrix(Bval_Buccal[,-1]), K = 3)
#estProp_RF <- RefFree_Deconvolution$Omega
# Move row names to first column #
#estProp_RF <- tibble::rownames_to_column(as.data.frame(estProp_RF), "SampleID")
# Rename columns #
#colnames(estProp_RF) <- c("SampleID", "RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal")
# Write estimated cell proportions #
#write.csv(estProp_RF, "ReffFree_CellProp_Estimations_Buccal_K3.csv")
# Read in previously estimated cell props #
estProp_RF <- read.csv("ReffFree_CellProp_Estimations_Buccal_K3.csv")[,-1]
# Merge with all data #
data <- merge(data, estProp_RF, by = "SampleID", all.x = TRUE)
# Horvath #
data_temp <- umx_residualize("mcc_Horvath1", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_Horvath1_BatchCell_Res <- data_temp$mcc_Horvath1
# Horvath age acceleration #
data_temp <- umx_residualize("mcc_Horvath1_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_Horvath1_AccMinus_BatchCell_Res <- data_temp$mcc_Horvath1_AccMinus
# Hannum age #
data_temp <- umx_residualize("mcc_Hannum", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_Hannum_BatchCell_Res <- data_temp$mcc_Hannum
# Hannum age acceleration #
data_temp <- umx_residualize("mcc_Hannum_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_Hannum_AccMinus_BatchCell_Res <- data_temp$mcc_Hannum_AccMinus
# PedBE age #
data_temp <- umx_residualize("mcc_PedBE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_PedBE_BatchCell_Res <- data_temp$mcc_PedBE
# PedBE age acceleration #
data_temp <- umx_residualize("mcc_PedBE_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_PedBE_AccMinus_BatchCell_Res <- data_temp$mcc_PedBE_AccMinus
# SkinBloodClock age #
data_temp <- umx_residualize("mcc_Horvath2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_Horvath2_BatchCell_Res <- data_temp$mcc_Horvath2
# SkinBloodClock age acceleration #
data_temp <- umx_residualize("mcc_Horvath2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_Horvath2_AccMinus_BatchCell_Res <- data_temp$mcc_Horvath2_AccMinus
# PhenoAge #
data_temp <- umx_residualize("mcc_PhenoAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_PhenoAge_BatchCell_Res <- data_temp$mcc_PhenoAge
# PhenoAge acceleration #
data_temp <- umx_residualize("mcc_PhenoAge_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_PhenoAge_AccMinus_BatchCell_Res <- data_temp$mcc_PhenoAge_AccMinus
# GrimAge2 #
data_temp <- umx_residualize("mcc_GrimAge2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_GrimAge2_BatchCell_Res <- data_temp$mcc_GrimAge2
# GrimAge2 acceleration #
data_temp <- umx_residualize("mcc_GrimAge2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_GrimAge2_AccMinus_BatchCell_Res <- data_temp$mcc_GrimAge2_AccMinus
# DunedinPoAm #
data_temp <- umx_residualize("PoAm", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$PoAm_BatchCell_Res <- data_temp$PoAm
# DunedinPACE #
data_temp <- umx_residualize("mcc_DunedinPACE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$mcc_DunedinPACE_BatchCell_Res <- data_temp$mcc_DunedinPACE
# DNAmTL #
data_temp <- umx_residualize("DNAmTL", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$DNAmTL_BatchCell_Res <- data_temp$DNAmTL
# DNAmTL age adjusted #
data_temp <- umx_residualize("DNAmTLAdjAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Buccal", "RF_Cell2_Buccal", "RF_Cell3_Buccal"), data = data)
data$DNAmTLAdjAge_BatchCell_Res <- data_temp$DNAmTLAdjAge

## Add Batch and Cell Composition Residualized Clocks ##
## Saliva (K=5) ##
#Bval_b1 <- readRDS("bVals_SampFilter_BMIQNorm_Batch1.RDS")
#Bval_b1 <- as.data.frame(Bval_b1)
#Bval_b1 <- tibble::rownames_to_column(Bval_b1, "ProbeID")
#Bval_b2 <- readRDS("bVals_SampFilter_BMIQNorm_Batch2.RDS")
#Bval_b2 <- as.data.frame(Bval_b2)
#Bval_b2 <- tibble::rownames_to_column(Bval_b2, "ProbeID")
#Bval_b3 <- readRDS("bVals_SampFilter_BMIQNorm_Batch3.RDS")
#Bval_b3 <- as.data.frame(Bval_b3)
#Bval_b3 <- tibble::rownames_to_column(Bval_b3, "ProbeID")
#Bval_b4 <- readRDS("bVals_SampFilter_BMIQNorm_Batch4.RDS")
#Bval_b4 <- as.data.frame(Bval_b4)
#Bval_b4 <- tibble::rownames_to_column(Bval_b4, "ProbeID")
# Combine all matrices #
#Bval_b1b2 <- merge(Bval_b1, Bval_b2, by = "ProbeID")
#Bval_b1b2b3 <- merge(Bval_b1b2, Bval_b3, by = "ProbeID")
#Bval <- merge(Bval_b1b2b3, Bval_b4, by = "ProbeID")
#Bval_b1 <- NULL
#Bval_b2 <- NULL
#Bval_b3 <- NULL
#Bval_b4 <- NULL
#Bval_b1b2 <- NULL
#Bval_b1b2b3 <- NULL
# Split Bval by tissue type for Saliva cells #
#keep <- colnames(Bval) %in% dplyr::filter(data, data$Tissue == "Saliva")$SampleID
#keep[1] <- TRUE
#Bval_Saliva <- Bval[,keep]
# Estimate cell proportions with reference-free method #
#RefFree_Deconvolution <- myRefFreeCellMix(as.matrix(Bval_Saliva[,-1]), K = 5)
#estProp_RF <- RefFree_Deconvolution$Omega
# Move row names to first column #
#estProp_RF <- tibble::rownames_to_column(as.data.frame(estProp_RF), "SampleID")
# Rename columns #
#colnames(estProp_RF) <- c("SampleID", "RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva")
# Write estimated cell proportions #
#write.csv(estProp_RF, "ReffFree_CellProp_Estimations_Saliva_K5.csv")
# Read in previously estimated cell props #
estProp_RF <- read.csv("ReffFree_CellProp_Estimations_Saliva_K5.csv")[,-1]
# Merge with all data #
data <- merge(data, estProp_RF, by = "SampleID", all.x = TRUE)
# Horvath #
data_temp <- umx_residualize("mcc_Horvath1", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva","RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_Horvath1_BatchCell_Res[i] <- data_temp$mcc_Horvath1[i]}
}
# Horvath age acceleration #
data_temp <- umx_residualize("mcc_Horvath1_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_Horvath1_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Horvath1_AccMinus[i]}
}
# Hannum age #
data_temp <- umx_residualize("mcc_Hannum", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_Hannum_BatchCell_Res[i] <- data_temp$mcc_Hannum[i]}
}
# Hannum age acceleration #
data_temp <- umx_residualize("mcc_Hannum_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_Hannum_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Hannum_AccMinus[i]}
}
# PedBE age #
data_temp <- umx_residualize("mcc_PedBE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_PedBE_BatchCell_Res[i] <- data_temp$mcc_PedBE[i]}
}
# PedBE age acceleration #
data_temp <- umx_residualize("mcc_PedBE_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_PedBE_AccMinus_BatchCell_Res[i] <- data_temp$mcc_PedBE_AccMinus[i]}
}
# SkinBloodClock age #
data_temp <- umx_residualize("mcc_Horvath2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_Horvath2_BatchCell_Res[i] <- data_temp$mcc_Horvath2[i]}
}
# SkinBloodClock age acceleration #
data_temp <- umx_residualize("mcc_Horvath2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_Horvath2_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Horvath2_AccMinus[i]}
}
# PhenoAge #
data_temp <- umx_residualize("mcc_PhenoAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_PhenoAge_BatchCell_Res[i] <- data_temp$mcc_PhenoAge[i]}
}
# PhenoAge acceleration #
data_temp <- umx_residualize("mcc_PhenoAge_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_PhenoAge_AccMinus_BatchCell_Res[i] <- data_temp$mcc_PhenoAge_AccMinus[i]}
}
# GrimAge2 #
data_temp <- umx_residualize("mcc_GrimAge2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_GrimAge2_BatchCell_Res[i] <- data_temp$mcc_GrimAge2[i]}
}
# GrimAge2 acceleration #
data_temp <- umx_residualize("mcc_GrimAge2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_GrimAge2_AccMinus_BatchCell_Res[i] <- data_temp$mcc_GrimAge2_AccMinus[i]}
}
# DunedinPoAm #
data_temp <- umx_residualize("PoAm", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$PoAm_BatchCell_Res[i] <- data_temp$PoAm[i]}
}
# DunedinPACE #
data_temp <- umx_residualize("mcc_DunedinPACE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$mcc_DunedinPACE_BatchCell_Res[i] <- data_temp$mcc_DunedinPACE[i]}
}
# DNAmTL #
data_temp <- umx_residualize("DNAmTL", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$DNAmTL_BatchCell_Res[i] <- data_temp$DNAmTL[i]}
}
# DNAmTL age adjusted #
data_temp <- umx_residualize("DNAmTLAdjAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","RF_Cell1_Saliva", "RF_Cell2_Saliva", "RF_Cell3_Saliva", "RF_Cell4_Saliva", "RF_Cell5_Saliva"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Saliva") {data$DNAmTLAdjAge_BatchCell_Res[i] <- data_temp$DNAmTLAdjAge[i]}
}

## Add Batch and Cell Composition Residualized Clocks ##
## DBS (DNAm Clock Foundation Cell Compositions) ##
# Horvath #
data_temp <- umx_residualize("mcc_Horvath1", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_Horvath1_BatchCell_Res[i] <- data_temp$mcc_Horvath1[i]}
}
# Horvath age acceleration #
data_temp <- umx_residualize("mcc_Horvath1_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_Horvath1_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Horvath1_AccMinus[i]}
}
# Hannum age #
data_temp <- umx_residualize("mcc_Hannum", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_Hannum_BatchCell_Res[i] <- data_temp$mcc_Hannum[i]}
}
# Hannum age acceleration #
data_temp <- umx_residualize("mcc_Hannum_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_Hannum_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Hannum_AccMinus[i]}
}
# PedBE age #
data_temp <- umx_residualize("mcc_PedBE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_PedBE_BatchCell_Res[i] <- data_temp$mcc_PedBE[i]}
}
# PedBE age acceleration #
data_temp <- umx_residualize("mcc_PedBE_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_PedBE_AccMinus_BatchCell_Res[i] <- data_temp$mcc_PedBE_AccMinus[i]}
}
# SkinBloodClock age #
data_temp <- umx_residualize("mcc_Horvath2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_Horvath2_BatchCell_Res[i] <- data_temp$mcc_Horvath2[i]}
}
# SkinBloodClock age acceleration #
data_temp <- umx_residualize("mcc_Horvath2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_Horvath2_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Horvath2_AccMinus[i]}
}
# PhenoAge #
data_temp <- umx_residualize("mcc_PhenoAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_PhenoAge_BatchCell_Res[i] <- data_temp$mcc_PhenoAge[i]}
}
# PhenoAge acceleration #
data_temp <- umx_residualize("mcc_PhenoAge_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_PhenoAge_AccMinus_BatchCell_Res[i] <- data_temp$mcc_PhenoAge_AccMinus[i]}
}
# GrimAge2 #
data_temp <- umx_residualize("mcc_GrimAge2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_GrimAge2_BatchCell_Res[i] <- data_temp$mcc_GrimAge2[i]}
}
# GrimAge2 acceleration #
data_temp <- umx_residualize("mcc_GrimAge2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_GrimAge2_AccMinus_BatchCell_Res[i] <- data_temp$mcc_GrimAge2_AccMinus[i]}
}
# DunedinPoAm #
data_temp <- umx_residualize("PoAm", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$PoAm_BatchCell_Res[i] <- data_temp$PoAm[i]}
}
# DunedinPACE #
data_temp <- umx_residualize("mcc_DunedinPACE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$mcc_DunedinPACE_BatchCell_Res[i] <- data_temp$mcc_DunedinPACE[i]}
}
# DNAmTL #
data_temp <- umx_residualize("DNAmTL", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$DNAmTL_BatchCell_Res[i] <- data_temp$DNAmTL[i]}
}
# DNAmTL age adjusted #
data_temp <- umx_residualize("DNAmTLAdjAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "DBS") {data$DNAmTLAdjAge_BatchCell_Res[i] <- data_temp$DNAmTLAdjAge[i]}
}

## Add Batch and Cell Composition Residualized Clocks ##
## Buffy Coat (DNAm Clock Foundation Cell Compositions) ##
# Horvath #
data_temp <- umx_residualize("mcc_Horvath1", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_Horvath1_BatchCell_Res[i] <- data_temp$mcc_Horvath1[i]}
}
# Horvath age acceleration #
data_temp <- umx_residualize("mcc_Horvath1_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_Horvath1_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Horvath1_AccMinus[i]}
}
# Hannum age #
data_temp <- umx_residualize("mcc_Hannum", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_Hannum_BatchCell_Res[i] <- data_temp$mcc_Hannum[i]}
}
# Hannum age acceleration #
data_temp <- umx_residualize("mcc_Hannum_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_Hannum_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Hannum_AccMinus[i]}
}
# PedBE age #
data_temp <- umx_residualize("mcc_PedBE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_PedBE_BatchCell_Res[i] <- data_temp$mcc_PedBE[i]}
}
# PedBE age acceleration #
data_temp <- umx_residualize("mcc_PedBE_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_PedBE_AccMinus_BatchCell_Res[i] <- data_temp$mcc_PedBE_AccMinus[i]}
}
# SkinBloodClock age #
data_temp <- umx_residualize("mcc_Horvath2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_Horvath2_BatchCell_Res[i] <- data_temp$mcc_Horvath2[i]}
}
# SkinBloodClock age acceleration #
data_temp <- umx_residualize("mcc_Horvath2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_Horvath2_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Horvath2_AccMinus[i]}
}
# PhenoAge #
data_temp <- umx_residualize("mcc_PhenoAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_PhenoAge_BatchCell_Res[i] <- data_temp$mcc_PhenoAge[i]}
}
# PhenoAge acceleration #
data_temp <- umx_residualize("mcc_PhenoAge_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_PhenoAge_AccMinus_BatchCell_Res[i] <- data_temp$mcc_PhenoAge_AccMinus[i]}
}
# GrimAge2 #
data_temp <- umx_residualize("mcc_GrimAge2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_GrimAge2_BatchCell_Res[i] <- data_temp$mcc_GrimAge2[i]}
}
# GrimAge2 acceleration #
data_temp <- umx_residualize("mcc_GrimAge2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_GrimAge2_AccMinus_BatchCell_Res[i] <- data_temp$mcc_GrimAge2_AccMinus[i]}
}
# DunedinPoAm #
data_temp <- umx_residualize("PoAm", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$PoAm_BatchCell_Res[i] <- data_temp$PoAm[i]}
}
# DunedinPACE #
data_temp <- umx_residualize("mcc_DunedinPACE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$mcc_DunedinPACE_BatchCell_Res[i] <- data_temp$mcc_DunedinPACE[i]}
}
# DNAmTL #
data_temp <- umx_residualize("DNAmTL", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$DNAmTL_BatchCell_Res[i] <- data_temp$DNAmTL[i]}
}
# DNAmTL age adjusted #
data_temp <- umx_residualize("DNAmTLAdjAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono","Gran"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "Buffy Coat") {data$DNAmTLAdjAge_BatchCell_Res[i] <- data_temp$DNAmTLAdjAge[i]}
}

## Add Batch and Cell Composition Residualized Clocks ##
## PBMC (DNAm Clock Foundation Cell Compositions) ##
# Horvath #
# Horvath #
data_temp <- umx_residualize("mcc_Horvath1", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_Horvath1_BatchCell_Res[i] <- data_temp$mcc_Horvath1[i]}
}
# Horvath age acceleration #
data_temp <- umx_residualize("mcc_Horvath1_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_Horvath1_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Horvath1_AccMinus[i]}
}
# Hannum age #
data_temp <- umx_residualize("mcc_Hannum", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_Hannum_BatchCell_Res[i] <- data_temp$mcc_Hannum[i]}
}
# Hannum age acceleration #
data_temp <- umx_residualize("mcc_Hannum_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_Hannum_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Hannum_AccMinus[i]}
}
# PedBE age #
data_temp <- umx_residualize("mcc_PedBE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_PedBE_BatchCell_Res[i] <- data_temp$mcc_PedBE[i]}
}
# PedBE age acceleration #
data_temp <- umx_residualize("mcc_PedBE_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_PedBE_AccMinus_BatchCell_Res[i] <- data_temp$mcc_PedBE_AccMinus[i]}
}
# SkinBloodClock age #
data_temp <- umx_residualize("mcc_Horvath2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_Horvath2_BatchCell_Res[i] <- data_temp$mcc_Horvath2[i]}
}
# SkinBloodClock age acceleration #
data_temp <- umx_residualize("mcc_Horvath2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_Horvath2_AccMinus_BatchCell_Res[i] <- data_temp$mcc_Horvath2_AccMinus[i]}
}
# PhenoAge #
data_temp <- umx_residualize("mcc_PhenoAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_PhenoAge_BatchCell_Res[i] <- data_temp$mcc_PhenoAge[i]}
}
# PhenoAge acceleration #
data_temp <- umx_residualize("mcc_PhenoAge_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_PhenoAge_AccMinus_BatchCell_Res[i] <- data_temp$mcc_PhenoAge_AccMinus[i]}
}
# GrimAge2 #
data_temp <- umx_residualize("mcc_GrimAge2", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_GrimAge2_BatchCell_Res[i] <- data_temp$mcc_GrimAge2[i]}
}
# GrimAge2 acceleration #
data_temp <- umx_residualize("mcc_GrimAge2_AccMinus", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_GrimAge2_AccMinus_BatchCell_Res[i] <- data_temp$mcc_GrimAge2_AccMinus[i]}
}
# DunedinPoAm #
data_temp <- umx_residualize("PoAm", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$PoAm_BatchCell_Res[i] <- data_temp$PoAm[i]}
}
# DunedinPACE #
data_temp <- umx_residualize("mcc_DunedinPACE", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$mcc_DunedinPACE_BatchCell_Res[i] <- data_temp$mcc_DunedinPACE[i]}
}
# DNAmTL #
data_temp <- umx_residualize("DNAmTL", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$DNAmTL_BatchCell_Res[i] <- data_temp$DNAmTL[i]}
}
# DNAmTL age adjusted #
data_temp <- umx_residualize("DNAmTLAdjAge", covs = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30","CD8T","CD4T","NK","Bcell","Mono"), data = data)
for (i in 1:nrow(data)){
  if (data$Tissue[i] == "PBMC") {data$DNAmTLAdjAge_BatchCell_Res[i] <- data_temp$DNAmTLAdjAge[i]}
}

## Generate long dataframe ##
data_long <- reshape(data, idvar = "Individual", timevar = "Tissue", direction = "wide")
data_CHS <- dplyr::filter(data, data$Study_Cohort == "CHS")
data_TRN <- dplyr::filter(data, data$Study_Cohort == "Telomere")
data_long_CHS <- reshape(data_CHS, idvar = "Individual", timevar = "Tissue", direction = "wide")
data_long_TRN <- reshape(data_TRN, idvar = "Individual", timevar = "Tissue", direction = "wide")

## Add Dataframe with Bio and Phenotype Data to TRN ##
# Keep only phenotype variables wanted #
data_pheno_TRN <- data_pheno_TRN[,c(2,6,22,88,109,202)]
# Rename ID column #
colnames(data_pheno_TRN)[1] <- "SampleID"
# Add new variables to data_TRN #
data_TRN$BMI <- NA
data_TRN$PSS_Sum <- NA
data_TRN$BAI_Sum <- NA
data_TRN$BDI_Sum <- NA
data_TRN$LESS_Sum <- NA
# Merge dataframes #
for (i in 1:nrow(data_TRN)){
  for (j in 1:nrow(data_pheno_TRN)){
    if (data_TRN$Individual[i] == data_pheno_TRN$SampleID[j]){data_TRN[i,142:146] <- data_pheno_TRN[j,2:6]}
  }
}

## Add Dataframe with Bio and Phenotype Data to CHS ##
# Keep only phenotype variables wanted #
data_pheno_CHS <- data_pheno_CHS[,c(2,15)]
# Rename ID column #
colnames(data_pheno_CHS) <- c("SampleID","BMI")
# Add new variables to data_CHS #
data_CHS$BMI <- NA
# Merge dataframes #
for (i in 1:nrow(data_CHS)){
  for (j in 1:nrow(data_pheno_CHS)){
    if (data_CHS$Individual[i] == data_pheno_CHS$SampleID[j]){data_CHS[i,142] <- data_pheno_CHS[j,2]}
  }
}

## Add Dataframe with Bio and Phenotype Data to Main Data ##
# Add new variables to data #
data$BMI <- NA
# Merge dataframes - Step 1 #
for (i in 1:nrow(data)){
  for (j in 1:nrow(data_pheno_TRN)){
    if (data$Individual[i] == data_pheno_TRN$SampleID[j]){data[i,142] <- data_pheno_TRN[j,2]}
  }
}
# Merge dataframes - Step 2 #
for (i in 1:nrow(data)){
  for (j in 1:nrow(data_pheno_CHS)){
    if (data$Individual[i] == data_pheno_CHS$SampleID[j]){data[i,142] <- data_pheno_CHS[j,2]}
  }
}

## Remove unnecessary data frames ##
data_horv <- NULL
data_horv_b1 <- NULL
data_horv_b2 <- NULL
data_horv_b3 <- NULL
data_horv_pace <- NULL
data_horv_pace_PC <- NULL
data_horv_pace_PC_PBE <- NULL
data_PACE <- NULL
data_PACE_b1 <- NULL
data_PACE_b2 <- NULL
data_PACE_b3 <- NULL
data_PBE <- NULL
data_PBE_b1 <- NULL
data_PBE_b2 <- NULL
data_PBE_b3 <- NULL
data_PC <- NULL
data_PC_b1 <- NULL
data_PC_b2 <- NULL
data_PC_b3 <- NULL
data_PoAm <- NULL
data_PoAm_b1 <- NULL
data_PoAm_b2 <- NULL
data_PoAm_b3 <- NULL
data_temp <- NULL
data_TRN_pheno <- NULL
data_CHS_pheno <- NULL

### Variable Engineering/Creation ### ----
## Add plotting data frame ##
data_plot <- data
data_plot$Tissue <- as.character(data_plot$Tissue)
data_plot$Tissue[data_plot$Tissue == "Buffy Coat"] <- c("Buffy Coat\n(children only)")
data_plot$Tissue[data_plot$Tissue == "PBMC"] <- c("PBMC\n(adult only)")

## Reorder Tissue Factors ##
data$Tissue <- factor(data$Tissue, levels=c('Buccal', 'Saliva', 'DBS', 'Buffy Coat', 'PBMC'))
data_plot$Tissue <- factor(data_plot$Tissue, levels=c('Buccal', 'Saliva', 'DBS', 'Buffy Coat\n(children only)', 'PBMC\n(adult only)'))





### Testing Clock Foundation vs. Manual Calculation of Clocks ### ----
## Actual Age Estimates ##
# Horvath Pan-Tissue Clock #
cor(data$DNAmAge, data$mcc_Horvath1)
t.test(data$DNAmAge, data$mcc_Horvath1, paired = TRUE)
plot(data$DNAmAge, data$mcc_Horvath1, 
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Horvath Pan-Tissue Age")
abline(0,1)
# Hannum Clock #
cor(data$DNAmAgeHannum, data$mcc_Hannum)
t.test(data$DNAmAgeHannum, data$mcc_Hannum, paired = TRUE)
plot(data$DNAmAgeHannum, data$mcc_Hannum,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Hannum Age")
abline(0,1)
# PhenoAge #
cor(data$DNAmPhenoAge, data$mcc_PhenoAge)
t.test(data$DNAmPhenoAge, data$mcc_PhenoAge, paired = TRUE)
plot(data$DNAmPhenoAge, data$mcc_PhenoAge,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "PhenoAge")
abline(0,1)
# GrimAge2 #
cor(data$mcc_GrimAge2, data$mcc_GrimAge2)
t.test(data$mcc_GrimAge2, data$mcc_GrimAge2, paired = TRUE)
plot(data$mcc_GrimAge2, data$mcc_GrimAge2,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "GrimAge2")
abline(0,1)
# Horvath2 #
cor(data$DNAmAgeSkinBloodClock, data$mcc_Horvath2)
t.test(data$DNAmAgeSkinBloodClock, data$mcc_Horvath2, paired = TRUE)
plot(data$DNAmAgeSkinBloodClock, data$mcc_Horvath2,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Skin and Blood Age")
abline(0,1)
# PedBE #
cor(data$DNAmAgePedBE, data$mcc_PedBE)
t.test(data$DNAmAgePedBE, data$mcc_PedBE, paired = TRUE)
plot(data$DNAmAgePedBE, data$mcc_PedBE,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "PedBE Age")
abline(0,1)

## Acceleration Estimates (Difference Method) ##
# Horvath Pan-Tissue Clock #
cor(data$AgeAccelerationResidual, data$mcc_Horvath1_AccMinus)
t.test(data$AgeAccelerationResidual, data$mcc_Horvath1_AccMinus, paired = TRUE)
plot(data$AgeAccelerationResidual, data$mcc_Horvath1_AccMinus, 
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Horvath Pan-Tissue Age Accel. (Difference)")
abline(0,1)
# Hannum Clock #
cor(data$AgeAccelerationResidualHannum, data$mcc_Hannum_AccMinus)
t.test(data$AgeAccelerationResidualHannum, data$mcc_Hannum_AccMinus, paired = TRUE)
plot(data$AgeAccelerationResidualHannum, data$mcc_Hannum_AccMinus,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Hannum Accel. (Difference)")
abline(0,1)
# PhenoAge #
cor(data$AgeAccelPheno, data$mcc_PhenoAge_AccMinus)
t.test(data$AgeAccelPheno, data$mcc_PhenoAge_AccMinus, paired = TRUE)
plot(data$AgeAccelPheno, data$mcc_PhenoAge_AccMinus,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "PhenoAge Accel. (Difference)")
abline(0,1)
# GrimAge2 #
cor(data$AgeAccelGrim2, data$mcc_GrimAge2_AccMinus)
t.test(data$AgeAccelGrim2, data$mcc_GrimAge2_AccMinus, paired = TRUE)
plot(data$AgeAccelGrim2, data$mcc_GrimAge2_AccMinus,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "GrimAge2 Accel. (Difference)")
abline(0,1)
# Horvath2 #
cor(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_AccMinus)
t.test(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_AccMinus, paired = TRUE)
plot(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_AccMinus,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Skin and Blood Accel. (Difference)")
abline(0,1)
# PedBE #
cor(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_AccMinus)
t.test(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_AccMinus, paired = TRUE)
plot(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_AccMinus,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "PedBE Accel. (Difference)")
abline(0,1)

## Acceleration Estimates (Tissue-Agnostic Intercept-Included Method) ##
# Horvath Pan-Tissue Clock #
cor(data$AgeAccelerationResidual, data$mcc_Horvath1_TissAg_Int)
t.test(data$AgeAccelerationResidual, data$mcc_Horvath1_TissAg_Int, paired = TRUE)
plot(data$AgeAccelerationResidual, data$mcc_Horvath1_TissAg_Int, 
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Horvath Pan-Tissue Age Accel. (TisAg-Int)")
abline(0,1)
# Hannum Clock #
cor(data$AgeAccelerationResidualHannum, data$mcc_Hannum_TissAg_Int)
t.test(data$AgeAccelerationResidualHannum, data$mcc_Hannum_TissAg_Int, paired = TRUE)
plot(data$AgeAccelerationResidualHannum, data$mcc_Hannum_TissAg_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Hannum Accel. (TisAg-Int)")
abline(0,1)
# PhenoAge #
cor(data$AgeAccelPheno, data$mcc_PhenoAge_TissAg_Int)
t.test(data$AgeAccelPheno, data$mcc_PhenoAge_TissAg_Int, paired = TRUE)
plot(data$AgeAccelPheno, data$mcc_PhenoAge_TissAg_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "PhenoAge Accel. (TisAg-Int)")
abline(0,1)
# GrimAge2 #
cor(data$AgeAccelGrim2, data$mcc_GrimAge2_TissAg_Int)
t.test(data$AgeAccelGrim2, data$mcc_GrimAge2_TissAg_Int, paired = TRUE)
plot(data$AgeAccelGrim2, data$mcc_GrimAge2_TissAg_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "GrimAge2 Accel. (Difference)")
abline(0,1)
# Horvath2 #
cor(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_TissAg_Int)
t.test(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_TissAg_Int, paired = TRUE)
plot(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_TissAg_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Skin and Blood Accel. (Difference)")
abline(0,1)
# PedBE #
cor(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_TissAg_Int)
t.test(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_TissAg_Int, paired = TRUE)
plot(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_TissAg_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "PedBE Accel. (Difference)")
abline(0,1)

## Acceleration Estimates (Tissue-Specific Intercept-Included Method) ##
# Horvath Pan-Tissue Clock #
cor(data$AgeAccelerationResidual, data$mcc_Horvath1_TissS_Int)
t.test(data$AgeAccelerationResidual, data$mcc_Horvath1_TissS_Int, paired = TRUE)
plot(data$AgeAccelerationResidual, data$mcc_Horvath1_TissS_Int, 
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Horvath Pan-Tissue Age Accel. (TisS-Int)")
abline(0,1)
# Hannum Clock #
cor(data$AgeAccelerationResidualHannum, data$mccHannum_TissS_Int)
t.test(data$AgeAccelerationResidualHannum, data$mccHannum_TissS_Int, paired = TRUE)
plot(data$AgeAccelerationResidualHannum, data$mccHannum_TissS_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Hannum Accel. (TisS-Int)")
abline(0,1)
# PhenoAge #
cor(data$AgeAccelPheno, data$mccPhenoAge_TissS_Int)
t.test(data$AgeAccelPheno, data$mccPhenoAge_TissS_Int, paired = TRUE)
plot(data$AgeAccelPheno, data$mccPhenoAge_TissS_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "PhenoAge Accel. (TisS-Int)")
abline(0,1)
# GrimAge2 #
cor(data$AgeAccelGrim2, data$mcc_GrimAge2_TissS_Int)
t.test(data$AgeAccelGrim2, data$mcc_GrimAge2_TissS_Int, paired = TRUE)
plot(data$AgeAccelGrim2, data$mcc_GrimAge2_TissS_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "GrimAge2 Accel. (Difference)")
abline(0,1)
# Horvath2 #
cor(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_TissS_Int)
t.test(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_TissS_Int, paired = TRUE)
plot(data$AgeAccelerationResidual_SkinBloodClock, data$mcc_Horvath2_TissS_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "Skin and Blood Accel. (Difference)")
abline(0,1)
# PedBE #
cor(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_TissS_Int)
t.test(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_TissS_Int, paired = TRUE)
plot(data$AgeAccelerationResidualPedBE, data$mcc_PedBE_TissS_Int,
     xlab = "Clock Foundation", ylab = "Manual Computations",
     main = "PedBE Accel. (Difference)")
abline(0,1)

###############################################################
### Overall Descriptive Statistics ### ----
## Overall ##
# Total number of samples #
nrow(data)
# Number of individuals #
length(unique(data$Individual))
# Number of individuals from each cohort #
dplyr::filter(data, Study_Cohort == "Telomere") %>% group_by(Individual)
dplyr::filter(data, Study_Cohort == "CHS") %>% group_by(Individual)
# Number of samples from each tissue #
table(data$Tissue)
# Number of tissues per sample #
table(table(data$Individual))
# Average age of each cohort #
# Create data_baseline with only one line per person #
data_baseline <- data %>% 
  distinct(Individual, .keep_all = T)
# Telomere U01 #
dplyr::filter(data_baseline, Study_Cohort == "Telomere") %>% pull(Age) %>% range()
# CHS P50 #
dplyr::filter(data_baseline, Study_Cohort == "CHS") %>% pull(Age) %>% range()
# Average age of both cohorts combined #
data_baseline %>% pull(Age) %>% sd()
# Number of males and females in each cohort #
# Telomere U01 #
dplyr::filter(data_baseline, Study_Cohort == "Telomere") %>% pull(Female) %>% table()
# CHS P50 #
dplyr::filter(data_baseline, Study_Cohort == "CHS") %>% pull(Female) %>% table()
# Number of CHS samples that have EPIC v1 and EPIC v2 measurements for Buffy Coat #
data_CHS_All <- read.csv("Manual_SampleID_Bookkeeping2.csv")
data_CHS_All <- dplyr::filter(data_CHS_All, data_CHS_All$Visit == "V1")
data_CHS_buffy <- dplyr::filter(data_CHS, data_CHS$Tissue == "Buffy Coat")
data_CHS_buffy$SampleID <- substr(data_CHS_buffy$SampleID, 1, nchar(data_CHS_buffy$SampleID)-1)
keep <- data_CHS_All$SampleID %in% data_CHS_buffy$SampleID
data_CHS_All <- data_CHS_All[keep,]
data_CHS_All <- data_CHS_All[order(data_CHS_All$Sample_Combined),]
write.csv(data_CHS_All, "SampleSheet_CHS_EPICv1_Overlap.csv")
data_CHS_BC <- dplyr::filter(data_CHS, data_CHS$Tissue == "Buffy Coat")
write.csv(data_CHS_BC, "Clocks_CHS_Overlap_MultiTissue_Results.csv")


###############################################################
### Statistical Analyses - Horvath Clock Raw ### ----
## Descriptive Statistics - Horvath Age ##
# Full cohort #
describe(data_long$mcc_Horvath1.Buccal)
describe(data_long$mcc_Horvath1.Saliva)
describe(data_long$mcc_Horvath1.DBS)
describe(data_long$`mcc_Horvath1.Buffy Coat`)
describe(data_long$mcc_Horvath1.PBMC)
# TRN #
describe(data_long$mcc_Horvath1.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_Horvath1.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_Horvath1.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_Horvath1.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_Horvath1.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_Horvath1.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_Horvath1.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_Horvath1.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - Horvath Acceleration ##
# Full cohort #
describe(data_long$mcc_Horvath1_AccMinus.Buccal)
describe(data_long$mcc_Horvath1_AccMinus.Saliva)
describe(data_long$mcc_Horvath1_AccMinus.DBS)
describe(data_long$`mcc_Horvath1_AccMinus.Buffy Coat`)
describe(data_long$mcc_Horvath1_AccMinus.PBMC)
# TRN #
describe(data_long$mcc_Horvath1_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_Horvath1_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_Horvath1_AccMinus.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_Horvath1_AccMinus.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_Horvath1_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_Horvath1_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_Horvath1_AccMinus.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_Horvath1_AccMinus.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Violin Plots ##
# Violin plot Horvath Age stratified by tissue #
Horvath_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_Horvath1)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath1), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath1, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "Horvath Pan-Tissue Age", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age stratified by tissue and cohort #
# TRN #
Horvath_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_Horvath1)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath1), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath1, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "Horvath Pan-Tissue Age (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
Horvath_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_Horvath1)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath1), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath1, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "Horvath Pan-Tissue Age (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age Accel stratified by tissue #
Horvath_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_Horvath1_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath1_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath1_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "Horvath Pan-Tissue Age Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age Acceleration stratified by tissue and cohort #
# TRN #
Horvath_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_Horvath1_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath1_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath1_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "Horvath Pan-Tissue Age Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
Horvath_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_Horvath1_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath1_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath1_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "Horvath Pan-Tissue Age Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()


### Statistical Analyses - Hannum Clock Raw ### ----
## Descriptive Statistics - Hannum Age ##
# Full cohort #
describe(data_long$mcc_Hannum.Buccal)
describe(data_long$mcc_Hannum.Saliva)
describe(data_long$mcc_Hannum.DBS)
describe(data_long$`mcc_Hannum.Buffy Coat`)
describe(data_long$mcc_Hannum.PBMC)
# TRN #
describe(data_long$mcc_Hannum.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_Hannum.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_Hannum.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_Hannum.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_Hannum.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_Hannum.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_Hannum.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_Hannum.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - Hannum Acceleration ##
# Full cohort #
describe(data_long$mcc_Hannum_AccMinus.Buccal)
describe(data_long$mcc_Hannum_AccMinus.Saliva)
describe(data_long$mcc_Hannum_AccMinus.DBS)
describe(data_long$`mcc_Hannum_AccMinus.Buffy Coat`)
describe(data_long$mcc_Hannum_AccMinus.PBMC)
# TRN #
describe(data_long$mcc_Hannum_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_Hannum_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_Hannum_AccMinus.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_Hannum_AccMinus.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_Hannum_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_Hannum_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_Hannum_AccMinus.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_Hannum_AccMinus.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Violin Plots ##
# Violin plot Hannum Age stratified by tissue #
Hannum_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_Hannum)) +
  geom_violin(aes(x = Tissue, y = mcc_Hannum), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Hannum, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "Hannum Age", x = element_blank()) + 
  theme_bw()
# Violin plot Hannum Age stratified by tissue and cohort #
# TRN #
Hannum_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_Hannum)) +
  geom_violin(aes(x = Tissue, y = mcc_Hannum), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Hannum, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "Hannum Age (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
Hannum_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_Hannum)) +
  geom_violin(aes(x = Tissue, y = mcc_Hannum), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Hannum, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "Hannum Age (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot Hannum Age Accel stratified by tissue #
Hannum_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_Hannum_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Hannum_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Hannum_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "Hannum Age Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot Hannum Age Acceleration stratified by tissue and cohort #
# TRN #
Hannum_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_Hannum_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Hannum_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Hannum_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "Hannum Age Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
Hannum_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_Hannum_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Hannum_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Hannum_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "Hannum Age Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()

## Generate Combined First-Generation Plots ##
# All ages #
ggarrange(Horvath_All_RawAge, Horvath_All_AgeAccel,
          Hannum_All_RawAge, Hannum_All_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))
# TRN #
ggarrange(Horvath_TRN_RawAge, Horvath_TRN_AgeAccel,
          Hannum_TRN_RawAge, Hannum_TRN_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))
# CHS #
ggarrange(Horvath_CHS_RawAge, Horvath_CHS_AgeAccel,
          Hannum_CHS_RawAge, Hannum_CHS_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))

### Statistical Analyses - Skin and Blood Clock Raw ### ----
## Descriptive Statistics - Horvath Age ##
# Full cohort #
describe(data_long$mcc_Horvath2.Buccal)
describe(data_long$mcc_Horvath2.Saliva)
describe(data_long$mcc_Horvath2.DBS)
describe(data_long$`mcc_Horvath2.Buffy Coat`)
describe(data_long$mcc_Horvath2.PBMC)
# TRN #
describe(data_long$mcc_Horvath2.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_Horvath2.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_Horvath2.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_Horvath2.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_Horvath2.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_Horvath2.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_Horvath2.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_Horvath2.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - Horvath Acceleration ##
# Full cohort #
describe(data_long$mcc_Horvath2_AccMinus.Buccal)
describe(data_long$mcc_Horvath2_AccMinus.Saliva)
describe(data_long$mcc_Horvath2_AccMinus.DBS)
describe(data_long$`mcc_Horvath2_AccMinus.Buffy Coat`)
describe(data_long$mcc_Horvath2_AccMinus.PBMC)
# TRN #
describe(data_long$mcc_Horvath2_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_Horvath2_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_Horvath2_AccMinus.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_Horvath2_AccMinus.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_Horvath2_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_Horvath2_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_Horvath2_AccMinus.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_Horvath2_AccMinus.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Violin Plots ##
# Violin plot Horvath Age stratified by tissue #
SkinBlood_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_Horvath2)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "Skin and Blood Age", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age stratified by tissue and cohort #
# TRN #
SkinBlood_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_Horvath2)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "Skin and Blood Age (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
SkinBlood_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_Horvath2)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "Skin and Blood Age (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age Accel stratified by tissue #
SkinBlood_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_Horvath2_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath2_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath2_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "Skin and Blood Age Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age Acceleration stratified by tissue and cohort #
# TRN #
SkinBlood_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_Horvath2_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath2_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath2_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "Skin and Blood Age Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
SkinBlood_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_Horvath2_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_Horvath2_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_Horvath2_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "Skin and Blood Age Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()

## Generate Combined Plots ##
# TRN #
ggarrange(SkinBlood_TRN_RawAge, SkinBlood_TRN_AgeAccel,
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom",
          labels = c("A","B"))
# CHS #
ggarrange(SkinBlood_CHS_RawAge, SkinBlood_CHS_AgeAccel,
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom",
          labels = c("A","B"))

### Statistical Analyses - PedBE Clock Raw ### ----
## Descriptive Statistics - PedBE Age ##
# Full cohort #
describe(data_long$mcc_PedBE.Buccal)
describe(data_long$mcc_PedBE.Saliva)
describe(data_long$mcc_PedBE.DBS)
describe(data_long$`mcc_PedBE.Buffy Coat`)
describe(data_long$mcc_PedBE.PBMC)
# TRN #
describe(data_long$mcc_PedBE.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_PedBE.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_PedBE.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_PedBE.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_PedBE.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_PedBE.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_PedBE.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_PedBE.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - PedBE Acceleration ##
# Full cohort #
describe(data_long$mcc_PedBE_AccMinus.Buccal)
describe(data_long$mcc_PedBE_AccMinus.Saliva)
describe(data_long$mcc_PedBE_AccMinus.DBS)
describe(data_long$`mcc_PedBE_AccMinus.Buffy Coat`)
describe(data_long$mcc_PedBE_AccMinus.PBMC)
# TRN #
describe(data_long$mcc_PedBE_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_PedBE_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_PedBE_AccMinus.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_PedBE_AccMinus.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_PedBE_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_PedBE_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_PedBE_AccMinus.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_PedBE_AccMinus.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Violin Plots ##
# CHS #
PedBE_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_PedBE)) +
  geom_violin(aes(x = Tissue, y = mcc_PedBE), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_PedBE, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PedBE Age", x = element_blank()) + 
  theme_bw()
# CHS #
PedBE_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_PedBE_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_PedBE_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_PedBE_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PedBE Age Acceleration", x = element_blank()) + 
  theme_bw()

## Generate Combined Other Clocks Plots ##
# All ages #
ggarrange(SkinBlood_All_RawAge, SkinBlood_All_AgeAccel,
          PedBE_CHS_RawAge, PedBE_CHS_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))

### Statistical Analyses - PhenoAge Raw ### ----
## Descriptive Statistics - PhenoAge ##
# Full cohort #
describe(data_long$mcc_PhenoAge.Buccal)
describe(data_long$mcc_PhenoAge.Saliva)
describe(data_long$mcc_PhenoAge.DBS)
describe(data_long$`mcc_PhenoAge.Buffy Coat`)
describe(data_long$mcc_PhenoAge.PBMC)
# TRN #
describe(data_long$mcc_PhenoAge.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_PhenoAge.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_PhenoAge.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_PhenoAge.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_PhenoAge.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_PhenoAge.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_PhenoAge.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_PhenoAge.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - PhenoAge Acceleration ##
# Full cohort #
describe(data_long$mcc_PhenoAge_AccMinus.Buccal)
describe(data_long$mcc_PhenoAge_AccMinus.Saliva)
describe(data_long$mcc_PhenoAge_AccMinus.DBS)
describe(data_long$`mcc_PhenoAge_AccMinus.Buffy Coat`)
describe(data_long$mcc_PhenoAge_AccMinus.PBMC)
# TRN #
describe(data_long$mcc_PhenoAge_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_PhenoAge_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_PhenoAge_AccMinus.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_PhenoAge_AccMinus.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_PhenoAge_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_PhenoAge_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_PhenoAge_AccMinus.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_PhenoAge_AccMinus.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Violin Plots ##
# Violin plot PhenoAge stratified by tissue #
PhenoAge_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_PhenoAge)) +
  geom_violin(aes(x = Tissue, y = mcc_PhenoAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_PhenoAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PhenoAge", x = element_blank()) + 
  theme_bw()
# Violin plot PhenoAge stratified by tissue and cohort #
# TRN #
PhenoAge_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_PhenoAge)) +
  geom_violin(aes(x = Tissue, y = mcc_PhenoAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_PhenoAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PhenoAge (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
PhenoAge_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_PhenoAge)) +
  geom_violin(aes(x = Tissue, y = mcc_PhenoAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_PhenoAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PhenoAge (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot PhenoAge Accel stratified by tissue #
PhenoAge_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_PhenoAge_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_PhenoAge_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_PhenoAge_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PhenoAge Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot PhenoAge Acceleration stratified by tissue and cohort #
# TRN #
PhenoAge_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_PhenoAge_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_PhenoAge_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_PhenoAge_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PhenoAge Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
PhenoAge_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_PhenoAge_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_PhenoAge_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_PhenoAge_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PhenoAge Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()

### Statistical Analyses - GrimAge2 Raw ### ----
## Descriptive Statistics - GrimAge2 ##
# Full cohort #
describe(data_long$mcc_GrimAge2.Buccal)
describe(data_long$mcc_GrimAge2.Saliva)
describe(data_long$mcc_GrimAge2.DBS)
describe(data_long$`mcc_GrimAge2.Buffy Coat`)
describe(data_long$mcc_GrimAge2.PBMC)
# TRN #
describe(data_long$mcc_GrimAge2.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_GrimAge2.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_GrimAge2.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_GrimAge2.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_GrimAge2.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_GrimAge2.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_GrimAge2.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_GrimAge2.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - GrimAge2 Acceleration ##
# Full cohort #
describe(data_long$mcc_GrimAge2_AccMinus.Buccal)
describe(data_long$mcc_GrimAge2_AccMinus.Saliva)
describe(data_long$mcc_GrimAge2_AccMinus.DBS)
describe(data_long$`mcc_GrimAge2_AccMinus.Buffy Coat`)
describe(data_long$mcc_GrimAge2_AccMinus.PBMC)
# TRN #
describe(data_long$mcc_GrimAge2_AccMinus.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_GrimAge2_AccMinus.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_GrimAge2_AccMinus.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_GrimAge2_AccMinus.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$AgeAccelGrim.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$AgeAccelGrim.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$AgeAccelGrim.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`AgeAccelGrim.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Violin Plots ##
# Violin plot GrimAge2 stratified by tissue #
GrimAge2_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_GrimAge2)) +
  geom_violin(aes(x = Tissue, y = mcc_GrimAge2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_GrimAge2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "GrimAge2", x = element_blank()) + 
  theme_bw()
# Violin plot GrimAge2 stratified by tissue and cohort #
# TRN #
GrimAge2_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_GrimAge2)) +
  geom_violin(aes(x = Tissue, y = mcc_GrimAge2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_GrimAge2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "GrimAge2 (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
GrimAge2_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_GrimAge2)) +
  geom_violin(aes(x = Tissue, y = mcc_GrimAge2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_GrimAge2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "GrimAge2 (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot GrimAge2 Accel stratified by tissue #
GrimAge2_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = mcc_GrimAge2_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_GrimAge2_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_GrimAge2_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "GrimAge2 Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot GrimAge2 Acceleration stratified by tissue and cohort #
# TRN #
GrimAge2_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_GrimAge2_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_GrimAge2_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_GrimAge2_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "GrimAge2 Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
GrimAge2_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_GrimAge2_AccMinus)) +
  geom_violin(aes(x = Tissue, y = mcc_GrimAge2_AccMinus), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_GrimAge2_AccMinus, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "GrimAge2 Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()

## Generate Combined First-Generation Plots ##
# All ages #
ggarrange(PhenoAge_All_RawAge, PhenoAge_All_AgeAccel,
          GrimAge2_All_RawAge, GrimAge2_All_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))
# TRN #
ggarrange(PhenoAge_TRN_RawAge, PhenoAge_TRN_AgeAccel,
          GrimAge2_TRN_RawAge, GrimAge2_TRN_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))
# CHS #
ggarrange(PhenoAge_CHS_RawAge, PhenoAge_CHS_AgeAccel,
          GrimAge2_CHS_RawAge, GrimAge2_CHS_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))

### Statistical Analyses - DunedinPoAm Raw ### ----
## Descriptive Statistics - PoAm ##
# Full cohort #
describe(data_long$PoAm.Buccal)
describe(data_long$PoAm.Saliva)
describe(data_long$PoAm.DBS)
describe(data_long$`PoAm.Buffy Coat`)
describe(data_long$PoAm.PBMC)
# TRN #
describe(data_long$PoAm.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PoAm.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PoAm.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PoAm.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PoAm.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PoAm.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PoAm.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PoAm.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Histograms ##
# Histograms of PoAm stratified by tissue #
ggplot(data, aes(x=PoAm, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="PoAm Measurement", y="Count")
## Violin Plots ##
# Violin plot PoAm stratified by tissue #
ggplot(data = data, aes(x = Tissue, y = PoAm)) +
  geom_violin(aes(x = Tissue, y = PoAm), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PoAm, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = "PoAm") + 
  theme_bw()
# Violin plot PoAm stratified by tissue and cohort #
# TRN #
ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PoAm)) +
  geom_violin(aes(x = Tissue, y = PoAm), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PoAm, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = "PoAm") + 
  theme_bw()
# CHS #
ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PoAm)) +
  geom_violin(aes(x = Tissue, y = PoAm), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PoAm, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = "PoAm") + 
  theme_bw()
### Statistical Analyses - DunedinPACE Raw ### ----
## Descriptive Statistics - PACE ##
# Full cohort #
describe(data_long$mcc_DunedinPACE.Buccal)
describe(data_long$mcc_DunedinPACE.Saliva)
describe(data_long$mcc_DunedinPACE.DBS)
describe(data_long$`mcc_DunedinPACE.Buffy Coat`)
describe(data_long$mcc_DunedinPACE.PBMC)
# TRN #
describe(data_long$mcc_DunedinPACE.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$mcc_DunedinPACE.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$mcc_DunedinPACE.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$mcc_DunedinPACE.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$mcc_DunedinPACE.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$mcc_DunedinPACE.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$mcc_DunedinPACE.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`mcc_DunedinPACE.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Histograms ##
# Histograms of PACE stratified by tissue #
ggplot(data, aes(x=mcc_DunedinPACE, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="PACE Measurement", y="Count")
## Violin Plots ##
# Violin plot PACE stratified by tissue #
ggplot(data = data_plot, aes(x = Tissue, y = mcc_DunedinPACE)) +
  geom_violin(aes(x = Tissue, y = mcc_DunedinPACE), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_DunedinPACE, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "DunedinPACE", x = element_blank()) + 
  theme_bw() +
  theme(legend.position="bottom")
# Violin plot PACE stratified by tissue and cohort #
# TRN #
ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = mcc_DunedinPACE)) +
  geom_violin(aes(x = Tissue, y = mcc_DunedinPACE), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_DunedinPACE, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "DunedinPACE (Adults Only)", x = element_blank()) + 
  theme_bw() +
  theme(legend.position="bottom")
# CHS #
ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = mcc_DunedinPACE)) +
  geom_violin(aes(x = Tissue, y = mcc_DunedinPACE), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = mcc_DunedinPACE, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "DunedinPACE (Children Only)", x = element_blank()) + 
  theme_bw() +
  theme(legend.position="bottom")

### Statistical Analyses - DNAmTL Raw ### ----
## Descriptive Statistics - DNAmTL ##
# Full cohort #
describe(data_long$DNAmTL.Buccal)
describe(data_long$DNAmTL.Saliva)
describe(data_long$DNAmTL.DBS)
describe(data_long$`DNAmTL.Buffy Coat`)
describe(data_long$DNAmTL.PBMC)
# TRN #
describe(data_long$DNAmTL.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$DNAmTL.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$DNAmTL.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$DNAmTL.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$DNAmTL.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$DNAmTL.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$DNAmTL.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`DNAmTL.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - DNAmTL Age Adjusted ##
# Full cohort #
describe(data_long$DNAmTLAdjAge.Buccal)
describe(data_long$DNAmTLAdjAge.Saliva)
describe(data_long$DNAmTLAdjAge.DBS)
describe(data_long$`DNAmTLAdjAge.Buffy Coat`)
describe(data_long$DNAmTLAdjAge.PBMC)
# TRN #
describe(data_long$DNAmTLAdjAge.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$DNAmTLAdjAge.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$DNAmTLAdjAge.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$DNAmTLAdjAge.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$DNAmTLAdjAge.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$DNAmTLAdjAge.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$DNAmTLAdjAge.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`DNAmTLAdjAge.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Histograms ##
# Histograms of DNAmTL stratified by tissue #
ggplot(data, aes(x=DNAmTL, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="DNAmTL", y="Count")
# Histograms of DNAmTL Age Adjusted stratified by tissue #
ggplot(data, aes(x=DNAmTLAdjAge, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="DNAmTL Age-Adjusted", y="Count")
## Violin Plots ##
# Violin plot DNAmTL stratified by tissue #
ggplot(data = data, aes(x = Tissue, y = DNAmTL)) +
  geom_violin(aes(x = Tissue, y = DNAmTL), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = DNAmTL, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = "DNAmTL (kb)") + 
  theme_bw()
# Violin plot DNAmTL stratified by tissue and cohort #
# TRN #
ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = DNAmTL)) +
  geom_violin(aes(x = Tissue, y = DNAmTL), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = DNAmTL, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = "TRN DNAmTL") + 
  theme_bw()
# CHS #
ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = DNAmTL)) +
  geom_violin(aes(x = Tissue, y = DNAmTL), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = DNAmTL, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = "CHS DNAmTL") + 
  theme_bw()
# Violin plot DNAmTL Age Adjusted stratified by tissue #
ggplot(data = data, aes(x = Tissue, y = DNAmTLAdjAge)) +
  geom_violin(aes(x = Tissue, y = DNAmTLAdjAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = DNAmTLAdjAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = "DNAmTL Age Adjusted") + 
  theme_bw()
# Violin plot DNAmTL Age Adjusted stratified by tissue and cohort #
# TRN #
ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = DNAmTLAdjAge)) +
  geom_violin(aes(x = Tissue, y = DNAmTLAdjAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = DNAmTLAdjAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = "TRN DNAmTL Age Adjusted") + 
  theme_bw()
# CHS #
ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = DNAmTLAdjAge)) +
  geom_violin(aes(x = Tissue, y = DNAmTLAdjAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = DNAmTLAdjAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = "CHS DNAmTL Age Adjusted") + 
  theme_bw()

### Statistical Analyses - PC Horvath Clock Raw ### ----
## Descriptive Statistics - Horvath Age ##
# Full cohort #
describe(data_long$PCHorvath1.Buccal)
describe(data_long$PCHorvath1.Saliva)
describe(data_long$PCHorvath1.DBS)
describe(data_long$`PCHorvath1.Buffy Coat`)
describe(data_long$PCHorvath1.PBMC)
# TRN #
describe(data_long$PCHorvath1.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCHorvath1.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCHorvath1.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCHorvath1.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCHorvath1.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCHorvath1.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCHorvath1.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCHorvath1.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - Horvath Acceleration ##
# Full cohort #
describe(data_long$PCHorvath1Resid.Buccal)
describe(data_long$PCHorvath1Resid.Saliva)
describe(data_long$PCHorvath1Resid.DBS)
describe(data_long$`PCHorvath1Resid.Buffy Coat`)
describe(data_long$PCHorvath1Resid.PBMC)
# TRN #
describe(data_long$PCHorvath1Resid.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCHorvath1Resid.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCHorvath1Resid.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCHorvath1Resid.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCHorvath1Resid.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCHorvath1Resid.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCHorvath1Resid.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCHorvath1Resid.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Histograms ##
# Histograms of Horvath Age stratified by tissue #
ggplot(data, aes(x=PCHorvath1, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="Horvath Clock", y="Count")
# Histograms of Horvath Age Accel stratified by tissue #
ggplot(data, aes(x=PCHorvath1Resid, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="Horvath Clock Accel", y="Count")
## Violin Plots ##
# Violin plot Horvath Age stratified by tissue #
HorvathPC_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = PCHorvath1)) +
  geom_violin(aes(x = Tissue, y = PCHorvath1), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath1, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC Horvath Pan-Tissue Age", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age stratified by tissue and cohort #
# TRN #
HorvathPC_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCHorvath1)) +
  geom_violin(aes(x = Tissue, y = PCHorvath1), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath1, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC Horvath Pan-Tissue Age (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
HorvathPC_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCHorvath1)) +
  geom_violin(aes(x = Tissue, y = PCHorvath1), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath1, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC Horvath Pan-Tissue Age (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age Accel stratified by tissue #
HorvathPC_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = PCHorvath1Resid)) +
  geom_violin(aes(x = Tissue, y = PCHorvath1Resid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath1Resid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC Horvath Pan-Tissue Age Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot Horvath Age Acceleration stratified by tissue and cohort #
# TRN #
HorvathPC_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCHorvath1Resid)) +
  geom_violin(aes(x = Tissue, y = PCHorvath1Resid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath1Resid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC Horvath Pan-Tissue Age Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
HorvathPC_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCHorvath1Resid)) +
  geom_violin(aes(x = Tissue, y = PCHorvath1Resid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath1Resid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC Horvath Pan-Tissue Age Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()


### Statistical Analyses - PC Hannum Clock Raw ### ----
## Descriptive Statistics - Hannum Age ##
# Full cohort #
describe(data_long$PCHannum.Buccal)
describe(data_long$PCHannum.Saliva)
describe(data_long$PCHannum.DBS)
describe(data_long$`PCHannum.Buffy Coat`)
describe(data_long$PCHannum.PBMC)
# TRN #
describe(data_long$PCHannum.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCHannum.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCHannum.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCHannum.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCHannum.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCHannum.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCHannum.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCHannum.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - Hannum Acceleration ##
# Full cohort #
describe(data_long$PCHannumResid.Buccal)
describe(data_long$PCHannumResid.Saliva)
describe(data_long$PCHannumResid.DBS)
describe(data_long$`PCHannumResid.Buffy Coat`)
describe(data_long$PCHannumResid.PBMC)
# TRN #
describe(data_long$PCHannumResid.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCHannumResid.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCHannumResid.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCHannumResid.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCHannumResid.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCHannumResid.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCHannumResid.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCHannumResid.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Histograms ##
# Histograms of Hannum Age stratified by tissue #
ggplot(data, aes(x=PCHannum, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="Hannum Clock", y="Count")
# Histograms of Hannum Age Accel stratified by tissue #
ggplot(data, aes(x=PCHannumResid, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="Hannum Clock Accel", y="Count")
## Violin Plots ##
# Violin plot Hannum Age stratified by tissue #
HannumPC_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = PCHannum)) +
  geom_violin(aes(x = Tissue, y = PCHannum), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHannum, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC Hannum Age", x = element_blank()) + 
  theme_bw()
# Violin plot Hannum Age stratified by tissue and cohort #
# TRN #
HannumPC_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCHannum)) +
  geom_violin(aes(x = Tissue, y = PCHannum), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHannum, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC Hannum Age (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
HannumPC_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCHannum)) +
  geom_violin(aes(x = Tissue, y = PCHannum), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHannum, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC Hannum Age (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot Hannum Age Accel stratified by tissue #
HannumPC_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = PCHannumResid)) +
  geom_violin(aes(x = Tissue, y = PCHannumResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHannumResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC Hannum Age Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot Hannum Age Acceleration stratified by tissue and cohort #
# TRN #
HannumPC_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCHannumResid)) +
  geom_violin(aes(x = Tissue, y = PCHannumResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHannumResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC Hannum Age Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
HannumPC_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCHannumResid)) +
  geom_violin(aes(x = Tissue, y = PCHannumResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHannumResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC Hannum Age Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()

## Generate Combined First-Generation Plots ##
# All ages #
ggarrange(HorvathPC_All_RawAge, HorvathPC_All_AgeAccel,
          HannumPC_All_RawAge, HannumPC_All_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))
# TRN #
ggarrange(HorvathPC_TRN_RawAge, HorvathPC_TRN_AgeAccel,
          HannumPC_TRN_RawAge, HannumPC_TRN_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))
# CHS #
ggarrange(HorvathPC_CHS_RawAge, HorvathPC_CHS_AgeAccel,
          HannumPC_CHS_RawAge, HannumPC_CHS_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))


### Statistical Analyses - PC PhenoAge Raw ### ----
## Descriptive Statistics - PhenoAge ##
# Full cohort #
describe(data_long$PCPhenoAge.Buccal)
describe(data_long$PCPhenoAge.Saliva)
describe(data_long$PCPhenoAge.DBS)
describe(data_long$`PCPhenoAge.Buffy Coat`)
describe(data_long$PCPhenoAge.PBMC)
# TRN #
describe(data_long$PCPhenoAge.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCPhenoAge.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCPhenoAge.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCPhenoAge.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCPhenoAge.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCPhenoAge.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCPhenoAge.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCPhenoAge.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - PhenoAge Acceleration ##
# Full cohort #
describe(data_long$PCPhenoAgeResid.Buccal)
describe(data_long$PCPhenoAgeResid.Saliva)
describe(data_long$PCPhenoAgeResid.DBS)
describe(data_long$`PCPhenoAgeResid.Buffy Coat`)
describe(data_long$PCPhenoAgeResid.PBMC)
# TRN #
describe(data_long$PCPhenoAgeResid.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCPhenoAgeResid.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCPhenoAgeResid.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCPhenoAgeResid.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCPhenoAgeResid.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCPhenoAgeResid.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCPhenoAgeResid.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCPhenoAgeResid.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Histograms ##
# Histograms of PhenoAge stratified by tissue #
ggplot(data, aes(x=PCPhenoAge, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="PhenoAge", y="Count")
# Histograms of PhenoAge Accel stratified by tissue #
ggplot(data, aes(x=PCPhenoAgeResid, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="PhenoAge Accel", y="Count")
## Violin Plots ##
# Violin plot PhenoAge stratified by tissue #
PhenoAgePC_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = PCPhenoAge)) +
  geom_violin(aes(x = Tissue, y = PCPhenoAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCPhenoAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC PhenoAge", x = element_blank()) + 
  theme_bw()
# Violin plot PhenoAge stratified by tissue and cohort #
# TRN #
PhenoAgePC_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCPhenoAge)) +
  geom_violin(aes(x = Tissue, y = PCPhenoAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCPhenoAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC PhenoAge (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
PhenoAgePC_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCPhenoAge)) +
  geom_violin(aes(x = Tissue, y = PCPhenoAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCPhenoAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC PhenoAge (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot PhenoAge Accel stratified by tissue #
PhenoAgePC_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = PCPhenoAgeResid)) +
  geom_violin(aes(x = Tissue, y = PCPhenoAgeResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCPhenoAgeResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC PhenoAge Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot PhenoAge Acceleration stratified by tissue and cohort #
# TRN #
PhenoAgePC_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCPhenoAgeResid)) +
  geom_violin(aes(x = Tissue, y = PCPhenoAgeResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCPhenoAgeResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC PhenoAge Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
PhenoAgePC_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCPhenoAgeResid)) +
  geom_violin(aes(x = Tissue, y = PCPhenoAgeResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCPhenoAgeResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC PhenoAge Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()

### Statistical Analyses - PC GrimAge Raw ### ----
## Descriptive Statistics - GrimAge2 ##
# Full cohort #
describe(data_long$PCGrimAge.Buccal)
describe(data_long$PCGrimAge.Saliva)
describe(data_long$PCGrimAge.DBS)
describe(data_long$`PCGrimAge.Buffy Coat`)
describe(data_long$PCGrimAge.PBMC)
# TRN #
describe(data_long$PCGrimAge.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCGrimAge.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCGrimAge.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCGrimAge.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCGrimAge.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCGrimAge.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCGrimAge.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCGrimAge.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - GrimAge2 Acceleration ##
# Full cohort #
describe(data_long$PCGrimAgeResid.Buccal)
describe(data_long$PCGrimAgeResid.Saliva)
describe(data_long$PCGrimAgeResid.DBS)
describe(data_long$`PCGrimAgeResid.Buffy Coat`)
describe(data_long$PCGrimAgeResid.PBMC)
# TRN #
describe(data_long$PCGrimAgeResid.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCGrimAgeResid.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCGrimAgeResid.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCGrimAgeResid.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCGrimAgeResid.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCGrimAgeResid.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCGrimAgeResid.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCGrimAgeResid.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Histograms ##
# Histograms of GrimAge2 stratified by tissue #
ggplot(data, aes(x=PCGrimAge, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="GrimAge2", y="Count")
# Histograms of GrimAge2 Accel stratified by tissue #
ggplot(data, aes(x=PCGrimAgeResid, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="GrimAge2 Accel", y="Count")
## Violin Plots ##
# Violin plot GrimAge2 stratified by tissue #
GrimAgePC_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = PCGrimAge)) +
  geom_violin(aes(x = Tissue, y = PCGrimAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCGrimAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC GrimAge", x = element_blank()) + 
  theme_bw()
# Violin plot GrimAge2 stratified by tissue and cohort #
# TRN #
GrimAgePC_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCGrimAge)) +
  geom_violin(aes(x = Tissue, y = PCGrimAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCGrimAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC GrimAge (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
GrimAgePC_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCGrimAge)) +
  geom_violin(aes(x = Tissue, y = PCGrimAge), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCGrimAge, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC GrimAge (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot GrimAge2 Accel stratified by tissue #
GrimAgePC_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = PCGrimAgeResid)) +
  geom_violin(aes(x = Tissue, y = PCGrimAgeResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCGrimAgeResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC GrimAge Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot GrimAge2 Acceleration stratified by tissue and cohort #
# TRN #
GrimAgePC_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCGrimAgeResid)) +
  geom_violin(aes(x = Tissue, y = PCGrimAgeResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCGrimAgeResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC GrimAge Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
GrimAgePC_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCGrimAgeResid)) +
  geom_violin(aes(x = Tissue, y = PCGrimAgeResid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCGrimAgeResid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC GrimAge Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()

## Generate Combined First-Generation Plots ##
# All ages #
ggarrange(PhenoAgePC_All_RawAge, PhenoAgePC_All_AgeAccel,
          GrimAgePC_All_RawAge, GrimAgePC_All_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))
# TRN #
ggarrange(PhenoAgePC_TRN_RawAge, PhenoAgePC_TRN_AgeAccel,
          GrimAgePC_TRN_RawAge, GrimAgePC_TRN_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))
# CHS #
ggarrange(PhenoAgePC_CHS_RawAge, PhenoAgePC_CHS_AgeAccel,
          GrimAgePC_CHS_RawAge, GrimAgePC_CHS_AgeAccel,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
          labels = c("A","B","D","E"))

### Statistical Analyses - PC Skin and Blood Raw ### ----
## Descriptive Statistics - Skin and Blood ##
# Full cohort #
describe(data_long$PCHorvath2.Buccal)
describe(data_long$PCHorvath2.Saliva)
describe(data_long$PCHorvath2.DBS)
describe(data_long$`PCHorvath2.Buffy Coat`)
describe(data_long$PCHorvath2.PBMC)
# TRN #
describe(data_long$PCHorvath2.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCHorvath2.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCHorvath2.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCHorvath2.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCHorvath2.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCHorvath2.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCHorvath2.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCHorvath2.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Descriptive Statistics - Skin and Blood Acceleration ##
# Full cohort #
describe(data_long$PCHorvath2Resid.Buccal)
describe(data_long$PCHorvath2Resid.Saliva)
describe(data_long$PCHorvath2Resid.DBS)
describe(data_long$`PCHorvath2Resid.Buffy Coat`)
describe(data_long$PCHorvath2Resid.PBMC)
# TRN #
describe(data_long$PCHorvath2Resid.Buccal[data_long$Study_Cohort.Buccal == "Telomere"])
describe(data_long$PCHorvath2Resid.Saliva[data_long$Study_Cohort.Saliva == "Telomere"])
describe(data_long$PCHorvath2Resid.DBS[data_long$Study_Cohort.DBS == "Telomere"])
describe(data_long$PCHorvath2Resid.PBMC[data_long$Study_Cohort.PBMC == "Telomere"])
# CHS #
describe(data_long$PCHorvath2Resid.Buccal[data_long$Study_Cohort.Buccal == "CHS"])
describe(data_long$PCHorvath2Resid.Saliva[data_long$Study_Cohort.Saliva == "CHS"])
describe(data_long$PCHorvath2Resid.DBS[data_long$Study_Cohort.DBS == "CHS"])
describe(data_long$`PCHorvath2Resid.Buffy Coat`[data_long$`Study_Cohort.Buffy Coat` == "CHS"])
## Histograms ##
# Histograms of Skin and Blood stratified by tissue #
ggplot(data, aes(x=PCHorvath2, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="Skin and Blood", y="Count")
# Histograms of Skin and Blood Accel stratified by tissue #
ggplot(data, aes(x=PCHorvath2Resid, color=Tissue, fill=Tissue)) +
  geom_histogram(position="identity", alpha=0.5)+
  theme(legend.position="top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Tissue ~ .) +
  labs(x="Skin and Blood Accel", y="Count")
## Violin Plots ##
# Violin plot Skin and Blood stratified by tissue #
SkinBloodPC_All_RawAge <- ggplot(data = data_plot, aes(x = Tissue, y = PCHorvath2)) +
  geom_violin(aes(x = Tissue, y = PCHorvath2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC Skin and Blood Age", x = element_blank()) + 
  theme_bw()
# Violin plot Skin and Blood stratified by tissue and cohort #
# TRN #
SkinBloodPC_TRN_RawAge <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCHorvath2)) +
  geom_violin(aes(x = Tissue, y = PCHorvath2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC Skin and Blood Age (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
SkinBloodPC_CHS_RawAge <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCHorvath2)) +
  geom_violin(aes(x = Tissue, y = PCHorvath2), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath2, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC Skin and Blood Age (Children Only)", x = element_blank()) + 
  theme_bw()
# Violin plot Skin and Blood Accel stratified by tissue #
SkinBloodPC_All_AgeAccel <- ggplot(data = data_plot, aes(x = Tissue, y = PCHorvath2Resid)) +
  geom_violin(aes(x = Tissue, y = PCHorvath2Resid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath2Resid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = element_blank(), title = "PC Skin and Blood Age Acceleration", x = element_blank()) + 
  theme_bw()
# Violin plot Skin and Blood Acceleration stratified by tissue and cohort #
# TRN #
SkinBloodPC_TRN_AgeAccel <- ggplot(data = data[data$Study_Cohort == "Telomere",], aes(x = Tissue, y = PCHorvath2Resid)) +
  geom_violin(aes(x = Tissue, y = PCHorvath2Resid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath2Resid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","orange")) + 
  scale_color_manual(values=c("red","blue","green4","orange")) +
  labs(y = element_blank(), title = "PC Skin and Blood Age Acceleration (Adults Only)", x = element_blank()) + 
  theme_bw()
# CHS #
SkinBloodPC_CHS_AgeAccel <- ggplot(data = data[data$Study_Cohort == "CHS",], aes(x = Tissue, y = PCHorvath2Resid)) +
  geom_violin(aes(x = Tissue, y = PCHorvath2Resid), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = PCHorvath2Resid, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple")) + 
  scale_color_manual(values=c("red","blue","green4","purple")) +
  labs(y = element_blank(), title = "PC Skin and Blood Age Acceleration (Children Only)", x = element_blank()) + 
  theme_bw()

## Generate Combined First-Generation Plots ##
# All ages #
ggarrange(SkinBloodPC_All_RawAge, SkinBloodPC_All_AgeAccel,
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom",
          labels = c("A","B"))
# TRN #
ggarrange(SkinBloodPC_TRN_RawAge, SkinBloodPC_TRN_AgeAccel,
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom",
          labels = c("A","B"))
# CHS #
ggarrange(SkinBloodPC_CHS_RawAge, SkinBloodPC_CHS_AgeAccel,
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom",
          labels = c("A","B"))

### Statistical Analyses - Tissue EpiClock Correlations Raw ### ----
## Define color palette ##
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
## Horvath ##
# Horvath age acceleration #
Horvath_Age_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_Horvath1_AccMinus.Buccal","mcc_Horvath1_AccMinus.Saliva","mcc_Horvath1_AccMinus.DBS",
                                                 "mcc_Horvath1_AccMinus.Buffy Coat","mcc_Horvath1_AccMinus.PBMC")])
Horvath_Age_Acc_Tissue_cor$r[4,5] <- 0
Horvath_Age_Acc_Tissue_cor$r[5,4] <- 0
Horvath_Age_Acc_Tissue_cor$p[4,5] <- 1
Horvath_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(Horvath_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(Horvath_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(Horvath_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(Horvath_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(Horvath_Age_Acc_Tissue_cor$r, method="ellipse",  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## Hannum ##
# Hannum age acceleration #
Hannum_Age_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_Hannum_AccMinus.Buccal","mcc_Hannum_AccMinus.Saliva","mcc_Hannum_AccMinus.DBS",
                                                     "mcc_Hannum_AccMinus.Buffy Coat","mcc_Hannum_AccMinus.PBMC")])
Hannum_Age_Acc_Tissue_cor$r[4,5] <- 0
Hannum_Age_Acc_Tissue_cor$r[5,4] <- 0
Hannum_Age_Acc_Tissue_cor$p[4,5] <- 1
Hannum_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(Hannum_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(Hannum_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(Hannum_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(Hannum_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(Hannum_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## Skin and Blood ##
# SkinBlood age acceleration #
SkinBlood_Age_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_Horvath2_AccMinus.Buccal","mcc_Horvath2_AccMinus.Saliva","mcc_Horvath2_AccMinus.DBS",
                                                     "mcc_Horvath2_AccMinus.Buffy Coat","mcc_Horvath2_AccMinus.PBMC")])
SkinBlood_Age_Acc_Tissue_cor$r[4,5] <- 0
SkinBlood_Age_Acc_Tissue_cor$r[5,4] <- 0
SkinBlood_Age_Acc_Tissue_cor$p[4,5] <- 1
SkinBlood_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(SkinBlood_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(SkinBlood_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(SkinBlood_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(SkinBlood_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(SkinBlood_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = SkinBlood_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PedBE ##
# PedBE age acceleration #
PedBE_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_PedBE_AccMinus.Buccal","mcc_PedBE_AccMinus.Saliva","mcc_PedBE_AccMinus.DBS",
                                                    "mcc_PedBE_AccMinus.Buffy Coat")])
rownames(PedBE_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PedBE_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(PedBE_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PedBE_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(PedBE_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PhenoAge ##
# PhenoAge acceleration #
PhenoAge_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_PhenoAge_AccMinus.Buccal","mcc_PhenoAge_AccMinus.Saliva","mcc_PhenoAge_AccMinus.DBS",
                                                    "mcc_PhenoAge_AccMinus.Buffy Coat","mcc_PhenoAge_AccMinus.PBMC")])
PhenoAge_Acc_Tissue_cor$r[4,5] <- 0
PhenoAge_Acc_Tissue_cor$r[5,4] <- 0
PhenoAge_Acc_Tissue_cor$p[4,5] <- 1
PhenoAge_Acc_Tissue_cor$p[5,4] <- 1
rownames(PhenoAge_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PhenoAge_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(PhenoAge_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PhenoAge_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(PhenoAge_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## GrimAge2 ##
# GrimAge2 acceleration #
GrimAge2_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_GrimAge2_AccMinus.Buccal","mcc_GrimAge2_AccMinus.Saliva","mcc_GrimAge2_AccMinus.DBS",
                                                 "mcc_GrimAge2_AccMinus.Buffy Coat","mcc_GrimAge2_AccMinus.PBMC")])
GrimAge2_Acc_Tissue_cor$r[4,5] <- 0
GrimAge2_Acc_Tissue_cor$r[5,4] <- 0
GrimAge2_Acc_Tissue_cor$p[4,5] <- 1
GrimAge2_Acc_Tissue_cor$p[5,4] <- 1
rownames(GrimAge2_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(GrimAge2_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(GrimAge2_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(GrimAge2_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(GrimAge2_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PoAm ##
# PoAm #
PoAm_Tissue_cor <- corr.test(data_long[,c("PoAm.Buccal","PoAm.Saliva","PoAm.DBS",
                                          "PoAm.Buffy Coat","PoAm.PBMC")])
PoAm_Tissue_cor$r[4,5] <- 0
PoAm_Tissue_cor$r[5,4] <- 0
PoAm_Tissue_cor$p[4,5] <- 1
PoAm_Tissue_cor$p[5,4] <- 1
rownames(PoAm_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PoAm_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(PoAm_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PoAm_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(PoAm_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)

## PACE ##
# PACE #
PACE_Tissue_cor <- corr.test(data_long[,c("mcc_DunedinPACE.Buccal","mcc_DunedinPACE.Saliva","mcc_DunedinPACE.DBS",
                                              "mcc_DunedinPACE.Buffy Coat","mcc_DunedinPACE.PBMC")])
PACE_Tissue_cor$r[4,5] <- 0
PACE_Tissue_cor$r[5,4] <- 0
PACE_Tissue_cor$p[4,5] <- 1
PACE_Tissue_cor$p[5,4] <- 1
rownames(PACE_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PACE_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(PACE_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PACE_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(PACE_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PACE_Tissue_cor$p, insig = "label_sig", sig.level = 0.01, pch = "*", pch.cex = 2,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## DNAmTL ##
# DNAmTL #
DNAmTL_Tissue_cor <- corr.test(data_long[,c("DNAmTL.Buccal","DNAmTL.Saliva","DNAmTL.DBS",
                                              "DNAmTL.Buffy Coat","DNAmTL.PBMC")])
DNAmTL_Tissue_cor$r[4,5] <- 0
DNAmTL_Tissue_cor$r[5,4] <- 0
DNAmTL_Tissue_cor$p[4,5] <- 1
DNAmTL_Tissue_cor$p[5,4] <- 1
rownames(DNAmTL_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(DNAmTL_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(DNAmTL_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(DNAmTL_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(DNAmTL_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)
# DNAmTL age adjusted #
DNAmTL_Adj_Tissue_cor <- corr.test(data_long[,c("DNAmTLAdjAge.Buccal","DNAmTLAdjAge.Saliva","DNAmTLAdjAge.DBS",
                                                  "DNAmTLAdjAge.Buffy Coat","DNAmTLAdjAge.PBMC")])
DNAmTL_Adj_Tissue_cor$r[4,5] <- 0
DNAmTL_Adj_Tissue_cor$r[5,4] <- 0
DNAmTL_Adj_Tissue_cor$p[4,5] <- 1
DNAmTL_Adj_Tissue_cor$p[5,4] <- 1
rownames(DNAmTL_Adj_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(DNAmTL_Adj_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(DNAmTL_Adj_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(DNAmTL_Adj_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(DNAmTL_Adj_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)

## PC Horvath ##
# Horvath age acceleration #
HorvathPC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCHorvath1Resid.Buccal","PCHorvath1Resid.Saliva","PCHorvath1Resid.DBS",
                                                     "PCHorvath1Resid.Buffy Coat","PCHorvath1Resid.PBMC")])
HorvathPC_Age_Acc_Tissue_cor$r[4,5] <- 0
HorvathPC_Age_Acc_Tissue_cor$r[5,4] <- 0
HorvathPC_Age_Acc_Tissue_cor$p[4,5] <- 1
HorvathPC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(HorvathPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(HorvathPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(HorvathPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(HorvathPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(HorvathPC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PC Hannum ##
# Hannum age acceleration #
HannumPC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCHannumResid.Buccal","PCHannumResid.Saliva","PCHannumResid.DBS",
                                                       "PCHannumResid.Buffy Coat","PCHannumResid.PBMC")])
HannumPC_Age_Acc_Tissue_cor$r[4,5] <- 0
HannumPC_Age_Acc_Tissue_cor$r[5,4] <- 0
HannumPC_Age_Acc_Tissue_cor$p[4,5] <- 1
HannumPC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(HannumPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(HannumPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(HannumPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(HannumPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(HannumPC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PC PhenoAge ##
# PhenoAge age acceleration #
PhenoAgePC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCPhenoAgeResid.Buccal","PCPhenoAgeResid.Saliva","PCPhenoAgeResid.DBS",
                                                      "PCPhenoAgeResid.Buffy Coat","PCPhenoAgeResid.PBMC")])
PhenoAgePC_Age_Acc_Tissue_cor$r[4,5] <- 0
PhenoAgePC_Age_Acc_Tissue_cor$r[5,4] <- 0
PhenoAgePC_Age_Acc_Tissue_cor$p[4,5] <- 1
PhenoAgePC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(PhenoAgePC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PhenoAgePC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(PhenoAgePC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PhenoAgePC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(PhenoAgePC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PC GrimAge ##
# GrimAge age acceleration #
GrimAgePC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCGrimAgeResid.Buccal","PCGrimAgeResid.Saliva","PCGrimAgeResid.DBS",
                                                        "PCGrimAgeResid.Buffy Coat","PCGrimAgeResid.PBMC")])
GrimAgePC_Age_Acc_Tissue_cor$r[4,5] <- 0
GrimAgePC_Age_Acc_Tissue_cor$r[5,4] <- 0
GrimAgePC_Age_Acc_Tissue_cor$p[4,5] <- 1
GrimAgePC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(GrimAgePC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(GrimAgePC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(GrimAgePC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(GrimAgePC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(GrimAgePC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PC Skin and Blood ##
# Skin and Blood age acceleration #
SkinBloodPC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCHorvath2Resid.Buccal","PCHorvath2Resid.Saliva","PCHorvath2Resid.DBS",
                                                       "PCHorvath2Resid.Buffy Coat","PCHorvath2Resid.PBMC")])
SkinBloodPC_Age_Acc_Tissue_cor$r[4,5] <- 0
SkinBloodPC_Age_Acc_Tissue_cor$r[5,4] <- 0
SkinBloodPC_Age_Acc_Tissue_cor$p[4,5] <- 1
SkinBloodPC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(SkinBloodPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(SkinBloodPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(SkinBloodPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(SkinBloodPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(SkinBloodPC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## Stratified by Cohort ##

## Horvath (Adults Only) ##
# Horvath age acceleration #
HorvathTRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_Horvath1_AccMinus.Buccal","mcc_Horvath1_AccMinus.Saliva","mcc_Horvath1_AccMinus.DBS",
                                                                  "mcc_Horvath1_AccMinus.PBMC")])
rownames(HorvathTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HorvathTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(HorvathTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HorvathTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(HorvathTRN_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## Horvath (Children Only) ##
# Horvath age acceleration #
HorvathCHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_Horvath1_AccMinus.Buccal","mcc_Horvath1_AccMinus.Saliva","mcc_Horvath1_AccMinus.DBS",
                                                                  "mcc_Horvath1_AccMinus.Buffy Coat")])
rownames(HorvathCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HorvathCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(HorvathCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HorvathCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(HorvathCHS_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## Hannum (Adults Only) ##
# Hannum age acceleration #
HannumTRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_Hannum_AccMinus.Buccal","mcc_Hannum_AccMinus.Saliva","mcc_Hannum_AccMinus.DBS",
                                                                 "mcc_Hannum_AccMinus.PBMC")])
rownames(HannumTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HannumTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(HannumTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HannumTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(HannumTRN_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## Hannum (Children Only) ##
# Hannum age acceleration #
HannumCHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_Hannum_AccMinus.Buccal","mcc_Hannum_AccMinus.Saliva","mcc_Hannum_AccMinus.DBS",
                                                                 "mcc_Hannum_AccMinus.Buffy Coat")])
rownames(HannumCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HannumCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(HannumCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HannumCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(HannumCHS_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PhenoAge (Adults Only) ##
# PhenoAge acceleration #
PhenoAgeTRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_PhenoAge_AccMinus.Buccal","mcc_PhenoAge_AccMinus.Saliva","mcc_PhenoAge_AccMinus.DBS",
                                                           "mcc_PhenoAge_AccMinus.PBMC")])
rownames(PhenoAgeTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PhenoAgeTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(PhenoAgeTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PhenoAgeTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(PhenoAgeTRN_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PhenoAge (Children Only) ##
# PhenoAge age acceleration #
PhenoAgeCHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_PhenoAge_AccMinus.Buccal","mcc_PhenoAge_AccMinus.Saliva","mcc_PhenoAge_AccMinus.DBS",
                                                           "mcc_PhenoAge_AccMinus.Buffy Coat")])
rownames(PhenoAgeCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PhenoAgeCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(PhenoAgeCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PhenoAgeCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(PhenoAgeCHS_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## GrimAge2 (Adults Only) ##
# GrimAge2 acceleration #
GrimAge2TRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_GrimAge2_AccMinus.Buccal","mcc_GrimAge2_AccMinus.Saliva","mcc_GrimAge2_AccMinus.DBS",
                                                             "mcc_GrimAge2_AccMinus.PBMC")])
rownames(GrimAge2TRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(GrimAge2TRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(GrimAge2TRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(GrimAge2TRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(GrimAge2TRN_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge2_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## GrimAge2 (Children Only) ##
# GrimAge2 age acceleration #
GrimAge2CHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_GrimAge2_AccMinus.Buccal","mcc_GrimAge2_AccMinus.Saliva","mcc_GrimAge2_AccMinus.DBS",
                                                             "mcc_GrimAge2_AccMinus.Buffy Coat")])
rownames(GrimAge2CHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(GrimAge2CHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(GrimAge2CHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(GrimAge2CHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(GrimAge2CHS_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge2_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PACE (Adults Only) ##
# PACE acceleration #
PACETRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_DunedinPACE.Buccal","mcc_DunedinPACE.Saliva","mcc_DunedinPACE.DBS",
                                                             "mcc_DunedinPACE.PBMC")])
rownames(PACETRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PACETRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(PACETRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PACETRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(PACETRN_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PACE_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PACE (Children Only) ##
# PACE age acceleration #
PACECHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_DunedinPACE.Buccal","mcc_DunedinPACE.Saliva","mcc_DunedinPACE.DBS",
                                                             "mcc_DunedinPACE.Buffy Coat")])
rownames(PACECHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PACECHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(PACECHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PACECHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(PACECHS_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PACE_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## SkinBlood (Adults Only) ##
# SkinBlood acceleration #
SkinBloodTRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_Horvath2_AccMinus.Buccal","mcc_Horvath2_AccMinus.Saliva","mcc_Horvath2_AccMinus.DBS",
                                                         "mcc_Horvath2_AccMinus.PBMC")])
rownames(SkinBloodTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(SkinBloodTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(SkinBloodTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(SkinBloodTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(SkinBloodTRN_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = SkinBlood_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## SkinBlood (Children Only) ##
# SkinBlood age acceleration #
SkinBloodCHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_Horvath2_AccMinus.Buccal","mcc_Horvath2_AccMinus.Saliva","mcc_Horvath2_AccMinus.DBS",
                                                         "mcc_Horvath2_AccMinus.Buffy Coat")])
rownames(SkinBloodCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(SkinBloodCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(SkinBloodCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(SkinBloodCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(SkinBloodCHS_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = SkinBlood_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Horvath (Adults Only) ##
# PC Horvath age acceleration #
HorvathTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCHorvath1Resid.Buccal","PCHorvath1Resid.Saliva","PCHorvath1Resid.DBS",
                                                       "PCHorvath1Resid.PBMC")])
rownames(HorvathTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HorvathTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(HorvathTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HorvathTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(HorvathTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Horvath (Children Only) ##
# PC Horvath age acceleration #
HorvathCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCHorvath1Resid.Buccal","PCHorvath1Resid.Saliva","PCHorvath1Resid.DBS",
                                                            "PCHorvath1Resid.Buffy Coat")])
rownames(HorvathCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HorvathCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(HorvathCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HorvathCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(HorvathCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Hannum (Adults Only) ##
# PC Hannum age acceleration #
HannumTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCHannumResid.Buccal","PCHannumResid.Saliva","PCHannumResid.DBS",
                                                            "PCHannumResid.PBMC")])
rownames(HannumTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HannumTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(HannumTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HannumTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(HannumTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Hannum (Children Only) ##
# PC Hannum age acceleration #
HannumCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCHannumResid.Buccal","PCHannumResid.Saliva","PCHannumResid.DBS",
                                                            "PCHannumResid.Buffy Coat")])
rownames(HannumCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HannumCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(HannumCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HannumCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(HannumCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC PhenoAge (Adults Only) ##
# PC PhenoAge age acceleration #
PhenoAgeTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCPhenoAgeResid.Buccal","PCPhenoAgeResid.Saliva","PCPhenoAgeResid.DBS",
                                                               "PCPhenoAgeResid.PBMC")])
rownames(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC PhenoAge (Children Only) ##
# PC PhenoAge age acceleration #
PhenoAgeCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCPhenoAgeResid.Buccal","PCPhenoAgeResid.Saliva","PCPhenoAgeResid.DBS",
                                                               "PCPhenoAgeResid.Buffy Coat")])
rownames(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC GrimAge (Adults Only) ##
# PC GrimAge age acceleration #
GrimAgeTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCGrimAgeResid.Buccal","PCGrimAgeResid.Saliva","PCGrimAgeResid.DBS",
                                                                "PCGrimAgeResid.PBMC")])
rownames(GrimAgeTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(GrimAgeTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(GrimAgeTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(GrimAgeTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(GrimAgeTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC GrimAge (Children Only) ##
# PC GrimAge age acceleration #
GrimAgeCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCGrimAgeResid.Buccal","PCGrimAgeResid.Saliva","PCGrimAgeResid.DBS",
                                                                "PCGrimAgeResid.Buffy Coat")])
rownames(GrimAgeCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(GrimAgeCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(GrimAgeCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(GrimAgeCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(GrimAgeCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Skin and Blood (Adults Only) ##
# PC Skin and Blood age acceleration #
SkinBloodTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCHorvath2Resid.Buccal","PCHorvath2Resid.Saliva","PCHorvath2Resid.DBS",
                                                               "PCHorvath2Resid.PBMC")])
rownames(SkinBloodTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(SkinBloodTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(SkinBloodTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(SkinBloodTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(SkinBloodTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Skin and Blood (Children Only) ##
# PC Skin and Blood age acceleration #
SkinBloodCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCHorvath2Resid.Buccal","PCHorvath2Resid.Saliva","PCHorvath2Resid.DBS",
                                                               "PCHorvath2Resid.Buffy Coat")])
rownames(SkinBloodCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(SkinBloodCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(SkinBloodCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(SkinBloodCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(SkinBloodCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse",
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)


### Statistical Analyses - Paired T-Testing Across Tissues ### ----
## Horvath ##
# Buccal vs Saliva #
t.test(data_long$mcc_Horvath1.Buccal, data_long$mcc_Horvath1.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$mcc_Horvath1.Buccal, data_long$mcc_Horvath1.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$mcc_Horvath1.Buccal, data_long$`mcc_Horvath1.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$mcc_Horvath1.Buccal, data_long$mcc_Horvath1.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$mcc_Horvath1.Saliva, data_long$mcc_Horvath1.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$mcc_Horvath1.Saliva, data_long$`mcc_Horvath1.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$mcc_Horvath1.Saliva, data_long$mcc_Horvath1.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$mcc_Horvath1.DBS, data_long$`mcc_Horvath1.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$mcc_Horvath1.DBS, data_long$mcc_Horvath1.PBMC, paired = TRUE)

## Hannum ##
# Buccal vs Saliva #
t.test(data_long$mcc_Hannum.Buccal, data_long$mcc_Hannum.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$mcc_Hannum.Buccal, data_long$mcc_Hannum.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$mcc_Hannum.Buccal, data_long$`mcc_Hannum.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$mcc_Hannum.Buccal, data_long$mcc_Hannum.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$mcc_Hannum.Saliva, data_long$mcc_Hannum.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$mcc_Hannum.Saliva, data_long$`mcc_Hannum.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$mcc_Hannum.Saliva, data_long$mcc_Hannum.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$mcc_Hannum.DBS, data_long$`mcc_Hannum.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$mcc_Hannum.DBS, data_long$mcc_Hannum.PBMC, paired = TRUE)

## PhenoAge ##
# Buccal vs Saliva #
t.test(data_long$mcc_PhenoAge.Buccal, data_long$mcc_PhenoAge.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$mcc_PhenoAge.Buccal, data_long$mcc_PhenoAge.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$mcc_PhenoAge.Buccal, data_long$`mcc_PhenoAge.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$mcc_PhenoAge.Buccal, data_long$mcc_PhenoAge.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$mcc_PhenoAge.Saliva, data_long$mcc_PhenoAge.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$mcc_PhenoAge.Saliva, data_long$`mcc_PhenoAge.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$mcc_PhenoAge.Saliva, data_long$mcc_PhenoAge.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$mcc_PhenoAge.DBS, data_long$`mcc_PhenoAge.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$mcc_PhenoAge.DBS, data_long$mcc_PhenoAge.PBMC, paired = TRUE)

## GrimAge2 ##
# Buccal vs Saliva #
t.test(data_long$mcc_GrimAge2.Buccal, data_long$mcc_GrimAge2.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$mcc_GrimAge2.Buccal, data_long$mcc_GrimAge2.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$mcc_GrimAge2.Buccal, data_long$`mcc_GrimAge2.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$mcc_GrimAge2.Buccal, data_long$mcc_GrimAge2.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$mcc_GrimAge2.Saliva, data_long$mcc_GrimAge2.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$mcc_GrimAge2.Saliva, data_long$`mcc_GrimAge2.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$mcc_GrimAge2.Saliva, data_long$mcc_GrimAge2.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$mcc_GrimAge2.DBS, data_long$`mcc_GrimAge2.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$mcc_GrimAge2.DBS, data_long$mcc_GrimAge2.PBMC, paired = TRUE)

## PACE ##
# Buccal vs Saliva #
t.test(data_long$mcc_DunedinPACE.Buccal, data_long$mcc_DunedinPACE.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$mcc_DunedinPACE.Buccal, data_long$mcc_DunedinPACE.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$mcc_DunedinPACE.Buccal, data_long$`mcc_DunedinPACE.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$mcc_DunedinPACE.Buccal, data_long$mcc_DunedinPACE.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$mcc_DunedinPACE.Saliva, data_long$mcc_DunedinPACE.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$mcc_DunedinPACE.Saliva, data_long$`mcc_DunedinPACE.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$mcc_DunedinPACE.Saliva, data_long$mcc_DunedinPACE.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$mcc_DunedinPACE.DBS, data_long$`mcc_DunedinPACE.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$mcc_DunedinPACE.DBS, data_long$mcc_DunedinPACE.PBMC, paired = TRUE)

## Skin and Blood ##
# Buccal vs Saliva #
t.test(data_long$mcc_Horvath2.Buccal, data_long$mcc_Horvath2.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$mcc_Horvath2.Buccal, data_long$mcc_Horvath2.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$mcc_Horvath2.Buccal, data_long$`mcc_Horvath2.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$mcc_Horvath2.Buccal, data_long$mcc_Horvath2.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$mcc_Horvath2.Saliva, data_long$mcc_Horvath2.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$mcc_Horvath2.Saliva, data_long$`mcc_Horvath2.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$mcc_Horvath2.Saliva, data_long$mcc_Horvath2.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$mcc_Horvath2.DBS, data_long$`mcc_Horvath2.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$mcc_Horvath2.DBS, data_long$mcc_Horvath2.PBMC, paired = TRUE)

## PedBE ##
# Buccal vs Saliva #
t.test(data_long_CHS$mcc_PedBE.Buccal, data_long_CHS$mcc_PedBE.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long_CHS$mcc_PedBE.Buccal, data_long_CHS$mcc_PedBE.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long_CHS$mcc_PedBE.Buccal, data_long_CHS$`mcc_PedBE.Buffy Coat`, paired = TRUE)
# Saliva vs DBS #
t.test(data_long_CHS$mcc_PedBE.Saliva, data_long_CHS$mcc_PedBE.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long_CHS$mcc_PedBE.Saliva, data_long_CHS$`mcc_PedBE.Buffy Coat`, paired = TRUE)
# DBS vs BC #
t.test(data_long_CHS$mcc_PedBE.DBS, data_long_CHS$`mcc_PedBE.Buffy Coat`, paired = TRUE)

### Statistical Analyses - Paired T-Testing Across Standard and PC Clocks Within-Tissue ### ----
## Horvath ##
# Buccal  #
t.test(data_long$mcc_Horvath1.Buccal, data_long$PCHorvath1.Buccal, paired = TRUE)
# Saliva #
t.test(data_long$mcc_Horvath1.Saliva, data_long$PCHorvath1.Saliva, paired = TRUE)
# DBS #
t.test(data_long$mcc_Horvath1.DBS, data_long$PCHorvath1.DBS, paired = TRUE)
# Buffy Coat #
t.test(data_long$`mcc_Horvath1.Buffy Coat`, data_long$`PCHorvath1.Buffy Coat`, paired = TRUE)
# PBMC #
t.test(data_long$mcc_Horvath1.PBMC, data_long$PCHorvath1.PBMC, paired = TRUE)

## Hannum ##
# Buccal  #
t.test(data_long$mcc_Hannum.Buccal, data_long$PCHannum.Buccal, paired = TRUE)
# Saliva #
t.test(data_long$mcc_Hannum.Saliva, data_long$PCHannum.Saliva, paired = TRUE)
# DBS #
t.test(data_long$mcc_Hannum.DBS, data_long$PCHannum.DBS, paired = TRUE)
# Buffy Coat #
t.test(data_long$`mcc_Hannum.Buffy Coat`, data_long$`PCHannum.Buffy Coat`, paired = TRUE)
# PBMC #
t.test(data_long$mcc_Hannum.PBMC, data_long$PCHannum.PBMC, paired = TRUE)

## PhenoAge ##
# Buccal  #
t.test(data_long$mcc_PhenoAge.Buccal, data_long$PCPhenoAge.Buccal, paired = TRUE)
# Saliva #
t.test(data_long$mcc_PhenoAge.Saliva, data_long$PCPhenoAge.Saliva, paired = TRUE)
# DBS #
t.test(data_long$mcc_PhenoAge.DBS, data_long$PCPhenoAge.DBS, paired = TRUE)
# Buffy Coat #
t.test(data_long$`mcc_PhenoAge.Buffy Coat`, data_long$`PCPhenoAge.Buffy Coat`, paired = TRUE)
# PBMC #
t.test(data_long$mcc_PhenoAge.PBMC, data_long$PCPhenoAge.PBMC, paired = TRUE)

## GrimAge ##
# Buccal  #
t.test(data_long$mcc_GrimAge2.Buccal, data_long$PCGrimAge.Buccal, paired = TRUE)
# Saliva #
t.test(data_long$mcc_GrimAge2.Saliva, data_long$PCGrimAge.Saliva, paired = TRUE)
# DBS #
t.test(data_long$mcc_GrimAge2.DBS, data_long$PCGrimAge.DBS, paired = TRUE)
# Buffy Coat #
t.test(data_long$`mcc_GrimAge2.Buffy Coat`, data_long$`PCGrimAge.Buffy Coat`, paired = TRUE)
# PBMC #
t.test(data_long$mcc_GrimAge2.PBMC, data_long$PCGrimAge.PBMC, paired = TRUE)

## Skin and Blood ##
# Buccal  #
t.test(data_long$mcc_Horvath2.Buccal, data_long$PCHorvath2.Buccal, paired = TRUE)
# Saliva #
t.test(data_long$mcc_Horvath2.Saliva, data_long$PCHorvath2.Saliva, paired = TRUE)
# DBS #
t.test(data_long$mcc_Horvath2.DBS, data_long$PCHorvath2.DBS, paired = TRUE)
# Buffy Coat #
t.test(data_long$`mcc_Horvath2.Buffy Coat`, data_long$`PCHorvath2.Buffy Coat`, paired = TRUE)
# PBMC #
t.test(data_long$mcc_Horvath2.PBMC, data_long$PCHorvath2.PBMC, paired = TRUE)

### Statistical Analyses - Paired T-Testing Across Tissues (PC Clocks) ### ----
## Horvath ##
# Buccal vs Saliva #
t.test(data_long$PCHorvath1.Buccal, data_long$PCHorvath1.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$PCHorvath1.Buccal, data_long$PCHorvath1.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$PCHorvath1.Buccal, data_long$`PCHorvath1.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$PCHorvath1.Buccal, data_long$PCHorvath1.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$PCHorvath1.Saliva, data_long$PCHorvath1.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$PCHorvath1.Saliva, data_long$`PCHorvath1.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$PCHorvath1.Saliva, data_long$PCHorvath1.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$PCHorvath1.DBS, data_long$`PCHorvath1.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$PCHorvath1.DBS, data_long$PCHorvath1.PBMC, paired = TRUE)

## Hannum ##
# Buccal vs Saliva #
t.test(data_long$PCHannum.Buccal, data_long$PCHannum.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$PCHannum.Buccal, data_long$PCHannum.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$PCHannum.Buccal, data_long$`PCHannum.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$PCHannum.Buccal, data_long$PCHannum.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$PCHannum.Saliva, data_long$PCHannum.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$PCHannum.Saliva, data_long$`PCHannum.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$PCHannum.Saliva, data_long$PCHannum.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$PCHannum.DBS, data_long$`PCHannum.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$PCHannum.DBS, data_long$PCHannum.PBMC, paired = TRUE)

## PhenoAge ##
# Buccal vs Saliva #
t.test(data_long$PCPhenoAge.Buccal, data_long$PCPhenoAge.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$PCPhenoAge.Buccal, data_long$PCPhenoAge.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$PCPhenoAge.Buccal, data_long$`PCPhenoAge.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$PCPhenoAge.Buccal, data_long$PCPhenoAge.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$PCPhenoAge.Saliva, data_long$PCPhenoAge.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$PCPhenoAge.Saliva, data_long$`PCPhenoAge.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$PCPhenoAge.Saliva, data_long$PCPhenoAge.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$PCPhenoAge.DBS, data_long$`PCPhenoAge.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$PCPhenoAge.DBS, data_long$PCPhenoAge.PBMC, paired = TRUE)

## GrimAge2 ##
# Buccal vs Saliva #
t.test(data_long$PCGrimAge.Buccal, data_long$PCGrimAge.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$PCGrimAge.Buccal, data_long$PCGrimAge.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$PCGrimAge.Buccal, data_long$`PCGrimAge.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$PCGrimAge.Buccal, data_long$PCGrimAge.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$PCGrimAge.Saliva, data_long$PCGrimAge.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$PCGrimAge.Saliva, data_long$`PCGrimAge.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$PCGrimAge.Saliva, data_long$PCGrimAge.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$PCGrimAge.DBS, data_long$`PCGrimAge.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$PCGrimAge.DBS, data_long$PCGrimAge.PBMC, paired = TRUE)

## Skin and Blood ##
# Buccal vs Saliva #
t.test(data_long$PCHorvath2.Buccal, data_long$PCHorvath2.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long$PCHorvath2.Buccal, data_long$PCHorvath2.DBS, paired = TRUE)
# Buccal vs BC #
t.test(data_long$PCHorvath2.Buccal, data_long$`PCHorvath2.Buffy Coat`, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long$PCHorvath2.Buccal, data_long$PCHorvath2.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long$PCHorvath2.Saliva, data_long$PCHorvath2.DBS, paired = TRUE)
# Saliva vs BC #
t.test(data_long$PCHorvath2.Saliva, data_long$`PCHorvath2.Buffy Coat`, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long$PCHorvath2.Saliva, data_long$PCHorvath2.PBMC, paired = TRUE)
# DBS vs BC #
t.test(data_long$PCHorvath2.DBS, data_long$`PCHorvath2.Buffy Coat`, paired = TRUE)
# DBS vs PBMC #
t.test(data_long$PCHorvath2.DBS, data_long$PCHorvath2.PBMC, paired = TRUE)

### Statistical Analyses - Paired T-Testing Across Tissues (Adults) ### ----
## Horvath ##
# Buccal vs Saliva #
t.test(data_long_TRN$mcc_Horvath1.Buccal, data_long_TRN$mcc_Horvath1.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long_TRN$mcc_Horvath1.Buccal, data_long_TRN$mcc_Horvath1.DBS, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long_TRN$mcc_Horvath1.Buccal, data_long_TRN$mcc_Horvath1.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long_TRN$mcc_Horvath1.Saliva, data_long_TRN$mcc_Horvath1.DBS, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long_TRN$mcc_Horvath1.Saliva, data_long_TRN$mcc_Horvath1.PBMC, paired = TRUE)
# DBS vs PBMC #
t.test(data_long_TRN$mcc_Horvath1.DBS, data_long_TRN$mcc_Horvath1.PBMC, paired = TRUE)

## Hannum ##
# Buccal vs Saliva #
t.test(data_long_TRN$mcc_Hannum.Buccal, data_long_TRN$mcc_Hannum.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long_TRN$mcc_Hannum.Buccal, data_long_TRN$mcc_Hannum.DBS, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long_TRN$mcc_Hannum.Buccal, data_long_TRN$mcc_Hannum.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long_TRN$mcc_Hannum.Saliva, data_long_TRN$mcc_Hannum.DBS, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long_TRN$mcc_Hannum.Saliva, data_long_TRN$mcc_Hannum.PBMC, paired = TRUE)
# DBS vs PBMC #
t.test(data_long_TRN$mcc_Hannum.DBS, data_long_TRN$mcc_Hannum.PBMC, paired = TRUE)

## PhenoAge ##
# Buccal vs Saliva #
t.test(data_long_TRN$mcc_PhenoAge.Buccal, data_long_TRN$mcc_PhenoAge.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long_TRN$mcc_PhenoAge.Buccal, data_long_TRN$mcc_PhenoAge.DBS, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long_TRN$mcc_PhenoAge.Buccal, data_long_TRN$mcc_PhenoAge.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long_TRN$mcc_PhenoAge.Saliva, data_long_TRN$mcc_PhenoAge.DBS, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long_TRN$mcc_PhenoAge.Saliva, data_long_TRN$mcc_PhenoAge.PBMC, paired = TRUE)
# DBS vs PBMC #
t.test(data_long_TRN$mcc_PhenoAge.DBS, data_long_TRN$mcc_PhenoAge.PBMC, paired = TRUE)

## GrimAge2 ##
# Buccal vs Saliva #
t.test(data_long_TRN$mcc_GrimAge2.Buccal, data_long_TRN$mcc_GrimAge2.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long_TRN$mcc_GrimAge2.Buccal, data_long_TRN$mcc_GrimAge2.DBS, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long_TRN$mcc_GrimAge2.Buccal, data_long_TRN$mcc_GrimAge2.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long_TRN$mcc_GrimAge2.Saliva, data_long_TRN$mcc_GrimAge2.DBS, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long_TRN$mcc_GrimAge2.Saliva, data_long_TRN$mcc_GrimAge2.PBMC, paired = TRUE)
# DBS vs PBMC #
t.test(data_long_TRN$mcc_GrimAge2.DBS, data_long_TRN$mcc_GrimAge2.PBMC, paired = TRUE)

## PACE ##
# Buccal vs Saliva #
t.test(data_long_TRN$mcc_DunedinPACE.Buccal, data_long_TRN$mcc_DunedinPACE.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long_TRN$mcc_DunedinPACE.Buccal, data_long_TRN$mcc_DunedinPACE.DBS, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long_TRN$mcc_DunedinPACE.Buccal, data_long_TRN$mcc_DunedinPACE.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long_TRN$mcc_DunedinPACE.Saliva, data_long_TRN$mcc_DunedinPACE.DBS, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long_TRN$mcc_DunedinPACE.Saliva, data_long_TRN$mcc_DunedinPACE.PBMC, paired = TRUE)
# DBS vs PBMC #
t.test(data_long_TRN$mcc_DunedinPACE.DBS, data_long_TRN$mcc_DunedinPACE.PBMC, paired = TRUE)

## Skin and Blood ##
# Buccal vs Saliva #
t.test(data_long_TRN$mcc_Horvath2.Buccal, data_long_TRN$mcc_Horvath2.Saliva, paired = TRUE)
# Buccal vs DBS #
t.test(data_long_TRN$mcc_Horvath2.Buccal, data_long_TRN$mcc_Horvath2.DBS, paired = TRUE)
# Buccal vs PBMC #
t.test(data_long_TRN$mcc_Horvath2.Buccal, data_long_TRN$mcc_Horvath2.PBMC, paired = TRUE)
# Saliva vs DBS #
t.test(data_long_TRN$mcc_Horvath2.Saliva, data_long_TRN$mcc_Horvath2.DBS, paired = TRUE)
# Saliva vs PBMC #
t.test(data_long_TRN$mcc_Horvath2.Saliva, data_long_TRN$mcc_Horvath2.PBMC, paired = TRUE)
# DBS vs PBMC #
t.test(data_long_TRN$mcc_Horvath2.DBS, data_long_TRN$mcc_Horvath2.PBMC, paired = TRUE)


### Statistical Analyses - Paired T-Testing Across Tissues (Children) ### ----
## Horvath ##
# Buccal vs Saliva #
t.test(data_long_CHS$mcc_Horvath1.Buccal, data_long_CHS$mcc_Horvath1.Saliva, paired = TRUE)$p.value
# Buccal vs DBS #
t.test(data_long_CHS$mcc_Horvath1.Buccal, data_long_CHS$mcc_Horvath1.DBS, paired = TRUE)$p.value
# Buccal vs BC #
t.test(data_long_CHS$mcc_Horvath1.Buccal, data_long_CHS$`mcc_Horvath1.Buffy Coat`, paired = TRUE)$p.value
# Saliva vs DBS #
t.test(data_long_CHS$mcc_Horvath1.Saliva, data_long_CHS$mcc_Horvath1.DBS, paired = TRUE)$p.value
# Saliva vs BC #
t.test(data_long_CHS$mcc_Horvath1.Saliva, data_long_CHS$`mcc_Horvath1.Buffy Coat`, paired = TRUE)$p.value
# DBS vs BC #
t.test(data_long_CHS$mcc_Horvath1.DBS, data_long_CHS$`mcc_Horvath1.Buffy Coat`, paired = TRUE)$p.value

## Hannum ##
# Buccal vs Saliva #
t.test(data_long_CHS$mcc_Hannum.Buccal, data_long_CHS$mcc_Hannum.Saliva, paired = TRUE)$p.value
# Buccal vs DBS #
t.test(data_long_CHS$mcc_Hannum.Buccal, data_long_CHS$mcc_Hannum.DBS, paired = TRUE)$p.value
# Buccal vs BC #
t.test(data_long_CHS$mcc_Hannum.Buccal, data_long_CHS$`mcc_Hannum.Buffy Coat`, paired = TRUE)$p.value
# Saliva vs DBS #
t.test(data_long_CHS$mcc_Hannum.Saliva, data_long_CHS$mcc_Hannum.DBS, paired = TRUE)$p.value
# Saliva vs BC #
t.test(data_long_CHS$mcc_Hannum.Saliva, data_long_CHS$`mcc_Hannum.Buffy Coat`, paired = TRUE)$p.value
# DBS vs BC #
t.test(data_long_CHS$mcc_Hannum.DBS, data_long_CHS$`mcc_Hannum.Buffy Coat`, paired = TRUE)$p.value

## PhenoAge ##
# Buccal vs Saliva #
t.test(data_long_CHS$mcc_PhenoAge.Buccal, data_long_CHS$mcc_PhenoAge.Saliva, paired = TRUE)$p.value
# Buccal vs DBS #
t.test(data_long_CHS$mcc_PhenoAge.Buccal, data_long_CHS$mcc_PhenoAge.DBS, paired = TRUE)$p.value
# Buccal vs BC #
t.test(data_long_CHS$mcc_PhenoAge.Buccal, data_long_CHS$`mcc_PhenoAge.Buffy Coat`, paired = TRUE)$p.value
# Saliva vs DBS #
t.test(data_long_CHS$mcc_PhenoAge.Saliva, data_long_CHS$mcc_PhenoAge.DBS, paired = TRUE)$p.value
# Saliva vs BC #
t.test(data_long_CHS$mcc_PhenoAge.Saliva, data_long_CHS$`mcc_PhenoAge.Buffy Coat`, paired = TRUE)$p.value
# DBS vs BC #
t.test(data_long_CHS$mcc_PhenoAge.DBS, data_long_CHS$`mcc_PhenoAge.Buffy Coat`, paired = TRUE)$p.value

## GrimAge2 ##
# Buccal vs Saliva #
t.test(data_long_CHS$mcc_GrimAge2.Buccal, data_long_CHS$mcc_GrimAge2.Saliva, paired = TRUE)$p.value
# Buccal vs DBS #
t.test(data_long_CHS$mcc_GrimAge2.Buccal, data_long_CHS$mcc_GrimAge2.DBS, paired = TRUE)$p.value
# Buccal vs BC #
t.test(data_long_CHS$mcc_GrimAge2.Buccal, data_long_CHS$`mcc_GrimAge2.Buffy Coat`, paired = TRUE)$p.value
# Saliva vs DBS #
t.test(data_long_CHS$mcc_GrimAge2.Saliva, data_long_CHS$mcc_GrimAge2.DBS, paired = TRUE)$p.value
# Saliva vs BC #
t.test(data_long_CHS$mcc_GrimAge2.Saliva, data_long_CHS$`mcc_GrimAge2.Buffy Coat`, paired = TRUE)$p.value
# DBS vs BC #
t.test(data_long_CHS$mcc_GrimAge2.DBS, data_long_CHS$`mcc_GrimAge2.Buffy Coat`, paired = TRUE)$p.value

## PACE ##
# Buccal vs Saliva #
t.test(data_long_CHS$mcc_DunedinPACE.Buccal, data_long_CHS$mcc_DunedinPACE.Saliva, paired = TRUE)$p.value
# Buccal vs DBS #
t.test(data_long_CHS$mcc_DunedinPACE.Buccal, data_long_CHS$mcc_DunedinPACE.DBS, paired = TRUE)$p.value
# Buccal vs BC #
t.test(data_long_CHS$mcc_DunedinPACE.Buccal, data_long_CHS$`mcc_DunedinPACE.Buffy Coat`, paired = TRUE)$p.value
# Saliva vs DBS #
t.test(data_long_CHS$mcc_DunedinPACE.Saliva, data_long_CHS$mcc_DunedinPACE.DBS, paired = TRUE)$p.value
# Saliva vs BC #
t.test(data_long_CHS$mcc_DunedinPACE.Saliva, data_long_CHS$`mcc_DunedinPACE.Buffy Coat`, paired = TRUE)$p.value
# DBS vs BC #
t.test(data_long_CHS$mcc_DunedinPACE.DBS, data_long_CHS$`mcc_DunedinPACE.Buffy Coat`, paired = TRUE)$p.value

## Skin and Blood ##
# Buccal vs Saliva #
t.test(data_long_CHS$mcc_Horvath2.Buccal, data_long_CHS$mcc_Horvath2.Saliva, paired = TRUE)$p.value
# Buccal vs DBS #
t.test(data_long_CHS$mcc_Horvath2.Buccal, data_long_CHS$mcc_Horvath2.DBS, paired = TRUE)$p.value
# Buccal vs BC #
t.test(data_long_CHS$mcc_Horvath2.Buccal, data_long_CHS$`mcc_Horvath2.Buffy Coat`, paired = TRUE)$p.value
# Saliva vs DBS #
t.test(data_long_CHS$mcc_Horvath2.Saliva, data_long_CHS$mcc_Horvath2.DBS, paired = TRUE)$p.value
# Saliva vs BC #
t.test(data_long_CHS$mcc_Horvath2.Saliva, data_long_CHS$`mcc_Horvath2.Buffy Coat`, paired = TRUE)$p.value
# DBS vs BC #
t.test(data_long_CHS$mcc_Horvath2.DBS, data_long_CHS$`mcc_Horvath2.Buffy Coat`, paired = TRUE)$p.value

## PedBE ##
# Buccal vs Saliva #
t.test(data_long_CHS$mcc_PedBE.Buccal, data_long_CHS_CHS$mcc_PedBE.Saliva, paired = TRUE)$p.value
# Buccal vs DBS #
t.test(data_long_CHS$mcc_PedBE.Buccal, data_long_CHS_CHS$mcc_PedBE.DBS, paired = TRUE)$p.value
# Buccal vs BC #
t.test(data_long_CHS$mcc_PedBE.Buccal, data_long_CHS_CHS$`mcc_PedBE.Buffy Coat`, paired = TRUE)$p.value
# Saliva vs DBS #
t.test(data_long_CHS$mcc_PedBE.Saliva, data_long_CHS_CHS$mcc_PedBE.DBS, paired = TRUE)$p.value
# Saliva vs BC #
t.test(data_long_CHS$mcc_PedBE.Saliva, data_long_CHS_CHS$`mcc_PedBE.Buffy Coat`, paired = TRUE)$p.value
# DBS vs BC #
t.test(data_long_CHS$mcc_PedBE.DBS, data_long_CHS_CHS$`mcc_PedBE.Buffy Coat`, paired = TRUE)$p.value

###############################################################
### Statistical Analyses - Tissue EpiClock Correlations Controlling for Batch and Cell Compositions ### ----
## Residualize Values by Batch and Cell Composition ##

## Define color palette ##
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
## Horvath ##
# Horvath age acceleration #
Horvath_Age_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_Horvath1_AccMinus_BatchCell_Res.Buccal","mcc_Horvath1_AccMinus_BatchCell_Res.Saliva","mcc_Horvath1_AccMinus_BatchCell_Res.DBS",
                                                     "mcc_Horvath1_AccMinus_BatchCell_Res.Buffy Coat","mcc_Horvath1_AccMinus_BatchCell_Res.PBMC")])
Horvath_Age_Acc_Tissue_cor$r[4,5] <- 0
Horvath_Age_Acc_Tissue_cor$r[5,4] <- 0
Horvath_Age_Acc_Tissue_cor$p[4,5] <- 1
Horvath_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(Horvath_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(Horvath_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(Horvath_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(Horvath_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(Horvath_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)
         #title = "Horvath Pan-Tissue Acceleration (Batch and Cell Composition Res.)") 

## Hannum ##
# Hannum age acceleration #
Hannum_Age_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_Hannum_AccMinus_BatchCell_Res.Buccal","mcc_Hannum_AccMinus_BatchCell_Res.Saliva","mcc_Hannum_AccMinus_BatchCell_Res.DBS",
                                                    "mcc_Hannum_AccMinus_BatchCell_Res.Buffy Coat","mcc_Hannum_AccMinus_BatchCell_Res.PBMC")])
Hannum_Age_Acc_Tissue_cor$r[4,5] <- 0
Hannum_Age_Acc_Tissue_cor$r[5,4] <- 0
Hannum_Age_Acc_Tissue_cor$p[4,5] <- 1
Hannum_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(Hannum_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(Hannum_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(Hannum_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(Hannum_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(Hannum_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)
         #title = "Hannum Acceleration (Batch and Cell Composition Res.)")

## Skin and Blood ##
# SkinBlood age acceleration #
SkinBlood_Age_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_Horvath2_AccMinus_BatchCell_Res.Buccal","mcc_Horvath2_AccMinus_BatchCell_Res.Saliva","mcc_Horvath2_AccMinus_BatchCell_Res.DBS",
                                                       "mcc_Horvath2_AccMinus_BatchCell_Res.Buffy Coat","mcc_Horvath2_AccMinus_BatchCell_Res.PBMC")])
SkinBlood_Age_Acc_Tissue_cor$r[4,5] <- 0
SkinBlood_Age_Acc_Tissue_cor$r[5,4] <- 0
SkinBlood_Age_Acc_Tissue_cor$p[4,5] <- 1
SkinBlood_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(SkinBlood_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(SkinBlood_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(SkinBlood_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(SkinBlood_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(SkinBlood_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = SkinBlood_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)
         #title = "Skin and Blood Acceleration (Batch and Cell Composition Res.)")  

## PedBE ##
# PedBE age acceleration #
PedBE_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_PedBE_AccMinus_BatchCell_Res.Buccal","mcc_PedBE_AccMinus_BatchCell_Res.Saliva","mcc_PedBE_AccMinus_BatchCell_Res.DBS",
                                                       "mcc_PedBE_AccMinus_BatchCell_Res.Buffy Coat")])
rownames(PedBE_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PedBE_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(PedBE_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PedBE_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(PedBE_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)
         #title = "PedBE Acceleration (Batch and Cell Composition Res.)")

## PhenoAge ##
# PhenoAge acceleration #
PhenoAge_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_PhenoAge_AccMinus_BatchCell_Res.Buccal","mcc_PhenoAge_AccMinus_BatchCell_Res.Saliva","mcc_PhenoAge_AccMinus_BatchCell_Res.DBS",
                                                  "mcc_PhenoAge_AccMinus_BatchCell_Res.Buffy Coat","mcc_PhenoAge_AccMinus_BatchCell_Res.PBMC")])
PhenoAge_Acc_Tissue_cor$r[4,5] <- 0
PhenoAge_Acc_Tissue_cor$r[5,4] <- 0
PhenoAge_Acc_Tissue_cor$p[4,5] <- 1
PhenoAge_Acc_Tissue_cor$p[5,4] <- 1
rownames(PhenoAge_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PhenoAge_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(PhenoAge_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PhenoAge_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(PhenoAge_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)
         #title = "PhenoAge Acceleration (Batch and Cell Composition Res.)") 

## GrimAge2 ##
# GrimAge2 acceleration #
GrimAge2_Acc_Tissue_cor <- corr.test(data_long[,c("mcc_GrimAge2_AccMinus_BatchCell_Res.Buccal","mcc_GrimAge2_AccMinus_BatchCell_Res.Saliva","mcc_GrimAge2_AccMinus_BatchCell_Res.DBS",
                                                  "mcc_GrimAge2_AccMinus_BatchCell_Res.Buffy Coat","mcc_GrimAge2_AccMinus_BatchCell_Res.PBMC")])
GrimAge2_Acc_Tissue_cor$r[4,5] <- 0
GrimAge2_Acc_Tissue_cor$r[5,4] <- 0
GrimAge2_Acc_Tissue_cor$p[4,5] <- 1
GrimAge2_Acc_Tissue_cor$p[5,4] <- 1
rownames(GrimAge2_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(GrimAge2_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(GrimAge2_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(GrimAge2_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(GrimAge2_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)
         #title = "GrimAge2 Acceleration (Batch and Cell Composition Res.)") 

## PACE ##
# PACE #
PACE_Tissue_cor <- corr.test(data_long[,c("mcc_DunedinPACE_BatchCell_Res.Buccal","mcc_DunedinPACE_BatchCell_Res.Saliva","mcc_DunedinPACE_BatchCell_Res.DBS",
                                          "mcc_DunedinPACE_BatchCell_Res.Buffy Coat","mcc_DunedinPACE_BatchCell_Res.PBMC")])
PACE_Tissue_cor$r[4,5] <- 0
PACE_Tissue_cor$r[5,4] <- 0
PACE_Tissue_cor$p[4,5] <- 1
PACE_Tissue_cor$p[5,4] <- 1
rownames(PACE_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PACE_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(PACE_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PACE_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(PACE_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)
         #title = "DunedinPACE (Batch and Cell Composition Res.)") 

## DNAmTL ##
# DNAmTL #
DNAmTL_Tissue_cor <- corr.test(data_long[,c("DNAmTL_BatchCell_Res.Buccal","DNAmTL_BatchCell_Res.Saliva","DNAmTL_BatchCell_Res.DBS",
                                            "DNAmTL_BatchCell_Res.Buffy Coat","DNAmTL_BatchCell_Res.PBMC")])
DNAmTL_Tissue_cor$r[4,5] <- 0
DNAmTL_Tissue_cor$r[5,4] <- 0
DNAmTL_Tissue_cor$p[4,5] <- 1
DNAmTL_Tissue_cor$p[5,4] <- 1
rownames(DNAmTL_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(DNAmTL_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(DNAmTL_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(DNAmTL_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(DNAmTL_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)

# DNAmTL age adjusted #
DNAmTL_Adj_Tissue_cor <- corr.test(data_long[,c("DNAmTLAdjAge_BatchCell_Res.Buccal","DNAmTLAdjAge_BatchCell_Res.Saliva","DNAmTLAdjAge_BatchCell_Res.DBS",
                                                "DNAmTLAdjAge_BatchCell_Res.Buffy Coat","DNAmTLAdjAge_BatchCell_Res.PBMC")])
DNAmTL_Adj_Tissue_cor$r[4,5] <- 0
DNAmTL_Adj_Tissue_cor$r[5,4] <- 0
DNAmTL_Adj_Tissue_cor$p[4,5] <- 1
DNAmTL_Adj_Tissue_cor$p[5,4] <- 1
rownames(DNAmTL_Adj_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(DNAmTL_Adj_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(DNAmTL_Adj_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(DNAmTL_Adj_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(DNAmTL_Adj_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)

## PC Horvath ##
# Horvath age acceleration #
HorvathPC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCHorvath1Resid.Buccal","PCHorvath1Resid.Saliva","PCHorvath1Resid.DBS",
                                                       "PCHorvath1Resid.Buffy Coat","PCHorvath1Resid.PBMC")])
HorvathPC_Age_Acc_Tissue_cor$r[4,5] <- 0
HorvathPC_Age_Acc_Tissue_cor$r[5,4] <- 0
HorvathPC_Age_Acc_Tissue_cor$p[4,5] <- 1
HorvathPC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(HorvathPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(HorvathPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(HorvathPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(HorvathPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(HorvathPC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PC Hannum ##
# Hannum age acceleration #
HannumPC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCHannumResid.Buccal","PCHannumResid.Saliva","PCHannumResid.DBS",
                                                      "PCHannumResid.Buffy Coat","PCHannumResid.PBMC")])
HannumPC_Age_Acc_Tissue_cor$r[4,5] <- 0
HannumPC_Age_Acc_Tissue_cor$r[5,4] <- 0
HannumPC_Age_Acc_Tissue_cor$p[4,5] <- 1
HannumPC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(HannumPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(HannumPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(HannumPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(HannumPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(HannumPC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PC PhenoAge ##
# PhenoAge age acceleration #
PhenoAgePC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCPhenoAgeResid.Buccal","PCPhenoAgeResid.Saliva","PCPhenoAgeResid.DBS",
                                                        "PCPhenoAgeResid.Buffy Coat","PCPhenoAgeResid.PBMC")])
PhenoAgePC_Age_Acc_Tissue_cor$r[4,5] <- 0
PhenoAgePC_Age_Acc_Tissue_cor$r[5,4] <- 0
PhenoAgePC_Age_Acc_Tissue_cor$p[4,5] <- 1
PhenoAgePC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(PhenoAgePC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PhenoAgePC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(PhenoAgePC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(PhenoAgePC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(PhenoAgePC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PC GrimAge ##
# GrimAge age acceleration #
GrimAgePC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCGrimAgeResid.Buccal","PCGrimAgeResid.Saliva","PCGrimAgeResid.DBS",
                                                       "PCGrimAgeResid.Buffy Coat","PCGrimAgeResid.PBMC")])
GrimAgePC_Age_Acc_Tissue_cor$r[4,5] <- 0
GrimAgePC_Age_Acc_Tissue_cor$r[5,4] <- 0
GrimAgePC_Age_Acc_Tissue_cor$p[4,5] <- 1
GrimAgePC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(GrimAgePC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(GrimAgePC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(GrimAgePC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(GrimAgePC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(GrimAgePC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## PC Skin and Blood ##
# Skin and Blood age acceleration #
SkinBloodPC_Age_Acc_Tissue_cor <- corr.test(data_long[,c("PCHorvath2Resid.Buccal","PCHorvath2Resid.Saliva","PCHorvath2Resid.DBS",
                                                         "PCHorvath2Resid.Buffy Coat","PCHorvath2Resid.PBMC")])
SkinBloodPC_Age_Acc_Tissue_cor$r[4,5] <- 0
SkinBloodPC_Age_Acc_Tissue_cor$r[5,4] <- 0
SkinBloodPC_Age_Acc_Tissue_cor$p[4,5] <- 1
SkinBloodPC_Age_Acc_Tissue_cor$p[5,4] <- 1
rownames(SkinBloodPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(SkinBloodPC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
rownames(SkinBloodPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
colnames(SkinBloodPC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat","PBMC")
# Correlation plot #
corrplot(SkinBloodPC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1) 

## Stratified by Cohort ##

## Horvath (Adults Only) ##
# Horvath age acceleration #
HorvathTRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_Horvath1_AccMinus.Buccal","mcc_Horvath1_AccMinus.Saliva","mcc_Horvath1_AccMinus.DBS",
                                                            "mcc_Horvath1_AccMinus.PBMC")])
rownames(HorvathTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HorvathTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(HorvathTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HorvathTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(HorvathTRN_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## Horvath (Children Only) ##
# Horvath age acceleration #
HorvathCHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_Horvath1_AccMinus.Buccal","mcc_Horvath1_AccMinus.Saliva","mcc_Horvath1_AccMinus.DBS",
                                                            "mcc_Horvath1_AccMinus.Buffy Coat")])
rownames(HorvathCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HorvathCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(HorvathCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HorvathCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(HorvathCHS_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## Hannum (Adults Only) ##
# Hannum age acceleration #
HannumTRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_Hannum_AccMinus.Buccal","mcc_Hannum_AccMinus.Saliva","mcc_Hannum_AccMinus.DBS",
                                                           "mcc_Hannum_AccMinus.PBMC")])
rownames(HannumTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HannumTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(HannumTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HannumTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(HannumTRN_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## Hannum (Children Only) ##
# Hannum age acceleration #
HannumCHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_Hannum_AccMinus.Buccal","mcc_Hannum_AccMinus.Saliva","mcc_Hannum_AccMinus.DBS",
                                                           "mcc_Hannum_AccMinus.Buffy Coat")])
rownames(HannumCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HannumCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(HannumCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HannumCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(HannumCHS_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PhenoAge (Adults Only) ##
# PhenoAge acceleration #
PhenoAgeTRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_PhenoAge_AccMinus.Buccal","mcc_PhenoAge_AccMinus.Saliva","mcc_PhenoAge_AccMinus.DBS",
                                                             "mcc_PhenoAge_AccMinus.PBMC")])
rownames(PhenoAgeTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PhenoAgeTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(PhenoAgeTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PhenoAgeTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(PhenoAgeTRN_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PhenoAge (Children Only) ##
# PhenoAge age acceleration #
PhenoAgeCHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_PhenoAge_AccMinus.Buccal","mcc_PhenoAge_AccMinus.Saliva","mcc_PhenoAge_AccMinus.DBS",
                                                             "mcc_PhenoAge_AccMinus.Buffy Coat")])
rownames(PhenoAgeCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PhenoAgeCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(PhenoAgeCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PhenoAgeCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(PhenoAgeCHS_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## GrimAge2 (Adults Only) ##
# GrimAge2 acceleration #
GrimAge2TRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_GrimAge2_AccMinus.Buccal","mcc_GrimAge2_AccMinus.Saliva","mcc_GrimAge2_AccMinus.DBS",
                                                             "mcc_GrimAge2_AccMinus.PBMC")])
rownames(GrimAge2TRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(GrimAge2TRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(GrimAge2TRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(GrimAge2TRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(GrimAge2TRN_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge2_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## GrimAge2 (Children Only) ##
# GrimAge2 age acceleration #
GrimAge2CHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_GrimAge2_AccMinus.Buccal","mcc_GrimAge2_AccMinus.Saliva","mcc_GrimAge2_AccMinus.DBS",
                                                             "mcc_GrimAge2_AccMinus.Buffy Coat")])
rownames(GrimAge2CHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(GrimAge2CHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(GrimAge2CHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(GrimAge2CHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(GrimAge2CHS_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge2_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PACE (Adults Only) ##
# PACE acceleration #
PACETRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_DunedinPACE.Buccal","mcc_DunedinPACE.Saliva","mcc_DunedinPACE.DBS",
                                                         "mcc_DunedinPACE.PBMC")])
rownames(PACETRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PACETRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(PACETRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PACETRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(PACETRN_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PACE_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PACE (Children Only) ##
# PACE age acceleration #
PACECHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_DunedinPACE.Buccal","mcc_DunedinPACE.Saliva","mcc_DunedinPACE.DBS",
                                                         "mcc_DunedinPACE.Buffy Coat")])
rownames(PACECHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PACECHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(PACECHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PACECHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(PACECHS_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PACE_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## SkinBlood (Adults Only) ##
# SkinBlood acceleration #
SkinBloodTRN_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("mcc_Horvath2_AccMinus.Buccal","mcc_Horvath2_AccMinus.Saliva","mcc_Horvath2_AccMinus.DBS",
                                                              "mcc_Horvath2_AccMinus.PBMC")])
rownames(SkinBloodTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(SkinBloodTRN_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(SkinBloodTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(SkinBloodTRN_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(SkinBloodTRN_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = SkinBlood_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## SkinBlood (Children Only) ##
# SkinBlood age acceleration #
SkinBloodCHS_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("mcc_Horvath2_AccMinus.Buccal","mcc_Horvath2_AccMinus.Saliva","mcc_Horvath2_AccMinus.DBS",
                                                              "mcc_Horvath2_AccMinus.Buffy Coat")])
rownames(SkinBloodCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(SkinBloodCHS_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(SkinBloodCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(SkinBloodCHS_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(SkinBloodCHS_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = SkinBlood_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Horvath (Adults Only) ##
# PC Horvath age acceleration #
HorvathTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCHorvath1Resid.Buccal","PCHorvath1Resid.Saliva","PCHorvath1Resid.DBS",
                                                               "PCHorvath1Resid.PBMC")])
rownames(HorvathTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HorvathTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(HorvathTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HorvathTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(HorvathTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Horvath (Children Only) ##
# PC Horvath age acceleration #
HorvathCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCHorvath1Resid.Buccal","PCHorvath1Resid.Saliva","PCHorvath1Resid.DBS",
                                                               "PCHorvath1Resid.Buffy Coat")])
rownames(HorvathCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HorvathCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(HorvathCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HorvathCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(HorvathCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Hannum (Adults Only) ##
# PC Hannum age acceleration #
HannumTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCHannumResid.Buccal","PCHannumResid.Saliva","PCHannumResid.DBS",
                                                              "PCHannumResid.PBMC")])
rownames(HannumTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HannumTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(HannumTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(HannumTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(HannumTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Hannum (Children Only) ##
# PC Hannum age acceleration #
HannumCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCHannumResid.Buccal","PCHannumResid.Saliva","PCHannumResid.DBS",
                                                              "PCHannumResid.Buffy Coat")])
rownames(HannumCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HannumCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(HannumCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(HannumCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(HannumCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Hannum_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC PhenoAge (Adults Only) ##
# PC PhenoAge age acceleration #
PhenoAgeTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCPhenoAgeResid.Buccal","PCPhenoAgeResid.Saliva","PCPhenoAgeResid.DBS",
                                                                "PCPhenoAgeResid.PBMC")])
rownames(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(PhenoAgeTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC PhenoAge (Children Only) ##
# PC PhenoAge age acceleration #
PhenoAgeCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCPhenoAgeResid.Buccal","PCPhenoAgeResid.Saliva","PCPhenoAgeResid.DBS",
                                                                "PCPhenoAgeResid.Buffy Coat")])
rownames(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(PhenoAgeCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = PhenoAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC GrimAge (Adults Only) ##
# PC GrimAge age acceleration #
GrimAgeTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCGrimAgeResid.Buccal","PCGrimAgeResid.Saliva","PCGrimAgeResid.DBS",
                                                               "PCGrimAgeResid.PBMC")])
rownames(GrimAgeTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(GrimAgeTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(GrimAgeTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(GrimAgeTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(GrimAgeTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC GrimAge (Children Only) ##
# PC GrimAge age acceleration #
GrimAgeCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCGrimAgeResid.Buccal","PCGrimAgeResid.Saliva","PCGrimAgeResid.DBS",
                                                               "PCGrimAgeResid.Buffy Coat")])
rownames(GrimAgeCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(GrimAgeCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(GrimAgeCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(GrimAgeCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(GrimAgeCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = GrimAge_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Skin and Blood (Adults Only) ##
# PC Skin and Blood age acceleration #
SkinBloodTRN_PC_Age_Acc_Tissue_cor <- corr.test(data_long_TRN[,c("PCHorvath2Resid.Buccal","PCHorvath2Resid.Saliva","PCHorvath2Resid.DBS",
                                                                 "PCHorvath2Resid.PBMC")])
rownames(SkinBloodTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
colnames(SkinBloodTRN_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","PBMC")
rownames(SkinBloodTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
colnames(SkinBloodTRN_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","PBMC")
# Correlation plot #
corrplot(SkinBloodTRN_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)

## PC Skin and Blood (Children Only) ##
# PC Skin and Blood age acceleration #
SkinBloodCHS_PC_Age_Acc_Tissue_cor <- corr.test(data_long_CHS[,c("PCHorvath2Resid.Buccal","PCHorvath2Resid.Saliva","PCHorvath2Resid.DBS",
                                                                 "PCHorvath2Resid.Buffy Coat")])
rownames(SkinBloodCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(SkinBloodCHS_PC_Age_Acc_Tissue_cor$r) <- c("Buccal","Saliva","DBS","Buffy Coat")
rownames(SkinBloodCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
colnames(SkinBloodCHS_PC_Age_Acc_Tissue_cor$p) <- c("Buccal","Saliva","DBS","Buffy Coat")
# Correlation plot #
corrplot(SkinBloodCHS_PC_Age_Acc_Tissue_cor$r, method="ellipse", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=0, #Text label color and rotation
         cl.pos = "n",
         # Combine with significance
         #p.mat = Horvath_Age_Tissue_cor$p, insig = "pch", sig.level = 0.01, pch = "*", pch.cex = 0.5,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         # Change font size of text labels
         tl.cex = 1.1)


###############################################################
