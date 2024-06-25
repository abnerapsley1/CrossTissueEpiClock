### Manual Computation of Epigenetic Clock Estimates ### ----

### Load Packages ### ----
library(minfi)
library(wateRmelon)
library(tibble)
library(dplyr)
library(methylCIPHER)

##### COMPUTE AGES #######
### Loading in Data ### ----
CpGs_batch1 <- readRDS("bVals_SampFilter_ProbeFilt_BMIQNorm_Batch1.RDS")
CpGs_batch2 <- readRDS("bVals_SampFilter_ProbeFilt_BMIQNorm_Batch2.RDS")
CpGs_batch3 <- readRDS("bVals_SampFilter_ProbeFilt_BMIQNorm_Batch3.RDS")
CpGs_batch4 <- readRDS("bVals_SampFilter_ProbeFilt_BMIQNorm_Batch4.RDS")
# Remove problem probes across all batches #
Probe_Remove_Batch1 <- read.csv("Batch1_pvalFilt_Probes.csv")[,-1]
Probe_Remove_Batch1 <- substr(Probe_Remove_Batch1,1,nchar(Probe_Remove_Batch1)-5)
Probe_Remove_Batch2 <- read.csv("Batch2_pvalFilt_Probes.csv")[,-1]
Probe_Remove_Batch2 <- substr(Probe_Remove_Batch2,1,nchar(Probe_Remove_Batch2)-5)
Probe_Remove_Batch3 <- read.csv("Batch3_pvalFilt_Probes.csv")[,-1]
Probe_Remove_Batch3 <- substr(Probe_Remove_Batch3,1,nchar(Probe_Remove_Batch3)-5)
Probe_Remove_Batch4 <- read.csv("Batch4_pvalFilt_Probes.csv")[,-1]
Probe_Remove_Batch4 <- substr(Probe_Remove_Batch4,1,nchar(Probe_Remove_Batch4)-5)
# Keep unique probes #
ProbeRemove <- unique(c(Probe_Remove_Batch1, Probe_Remove_Batch2, Probe_Remove_Batch3, Probe_Remove_Batch4))
# Remove probes from all batches #
CpGs_batch1 <- CpGs_batch1[!(rownames(CpGs_batch1) %in% ProbeRemove),]
CpGs_batch2 <- CpGs_batch2[!(rownames(CpGs_batch2) %in% ProbeRemove),]
CpGs_batch3 <- CpGs_batch3[!(rownames(CpGs_batch3) %in% ProbeRemove),]
CpGs_batch4 <- CpGs_batch4[!(rownames(CpGs_batch4) %in% ProbeRemove),]
# Combine methylation data #
data <- as.data.frame(cbind(CpGs_batch1, CpGs_batch2, CpGs_batch3, CpGs_batch4))
CpGs_batch1 <- NULL
CpGs_batch2 <- NULL
CpGs_batch3 <- NULL
CpGs_batch4 <- NULL
# Remove problem samples #
keep <- !(colnames(data) %in% c("2023314", "2025214", "504", "54"))
data <- data[,keep]
colnames(data)[278] <- "2029215"
## Golden Standard Datasets ##
GoldenStandard_BuccalAdults_Vector <- readRDS("GoldenStandard_Buccal_Adults_Mean_v3.Rds")
GoldenStandard_BuccalChildren_Vector <- readRDS("GoldenStandard_Buccal_Children_Mean_v3.Rds")
GoldenStandard_SalivaAdults_Vector <- readRDS("GoldenStandard_Saliva_Adults_Mean_v3.Rds")
GoldenStandard_SalivaChildren_Vector <- readRDS("GoldenStandard_Saliva_Children_Mean_v3.Rds")
GoldenStandard_WholeBloodAdults_Vector <- readRDS("GoldenStandard_WholeBlood_Adults_Mean_v3.Rds")
GoldenStandard_WholeBloodChildren_Vector <- readRDS("GoldenStandard_WholeBlood_Children_Mean_v3.Rds")
GoldenStandard_PBMCAdults_Vector <- readRDS("GoldenStandard_PBMC_Adults_Mean_v3.Rds")
GoldenStandard_PBMCChildren_Vector <- readRDS("GoldenStandard_PBMC_Children_Mean_v3.Rds")


### Split Data by Tissue and Age ###
## Buccal - Adults ##
# Keep only buccal tissue #
meth_buccal_adults <- data %>% select(ends_with("3"))
# Keep only adults #
meth_buccal_adults <- meth_buccal_adults[,as.numeric(colnames(meth_buccal_adults)) < 10000]
## Buccal - Children ##
# Keep only buccal tissue #
meth_buccal_children <- data %>% select(ends_with("3"))
# Keep only children #
meth_buccal_children <- meth_buccal_children[,as.numeric(colnames(meth_buccal_children)) > 10000]
## Saliva - Adults ##
# Keep only saliva tissue #
meth_saliva_adults <- data %>% select(ends_with("4"))
# Keep only adults #
meth_saliva_adults <- meth_saliva_adults[,as.numeric(colnames(meth_saliva_adults)) < 10000]
## Saliva - Children ##
# Keep only saliva tissue #
meth_saliva_children <- data %>% select(ends_with("4"))
# Keep only children #
meth_saliva_children <- meth_saliva_children[,as.numeric(colnames(meth_saliva_children)) > 10000]
## DBS - Adults ##
# Keep only DBS tissue #
meth_dbs_adults <- data %>% select(ends_with("2"))
# Keep only adults #
meth_dbs_adults <- meth_dbs_adults[,as.numeric(colnames(meth_dbs_adults)) < 10000]
## DBS - Children ##
# Keep only DBS tissue #
meth_dbs_children <- data %>% select(ends_with("2"))
# Keep only children #
meth_dbs_children <- meth_dbs_children[,as.numeric(colnames(meth_dbs_children)) > 10000]
## Buffy Coat - Adults ##
# Keep only Buffy Coat tissue #
meth_buffycoat_adults <- data %>% select(ends_with("5"))
# Keep only adults #
meth_buffycoat_adults <- meth_buffycoat_adults[,as.numeric(colnames(meth_buffycoat_adults)) < 10000]
## Buffy Coat - Children ##
# Keep only Buffy Coat tissue #
meth_buffycoat_children <- data %>% select(ends_with("5"))
# Keep only children #
meth_buffycoat_children <- meth_buffycoat_children[,as.numeric(colnames(meth_buffycoat_children)) > 10000]
## PBMC - Adults ##
# Keep only PBMC tissue #
meth_pbmc_adults <- data %>% select(ends_with("1"))
# Keep only adults #
meth_pbmc_adults <- meth_pbmc_adults[,as.numeric(colnames(meth_pbmc_adults)) < 10000]
## PBMC - Children ##
# Keep only PBMC tissue #
meth_pbmc_children <- data %>% select(ends_with("1"))
# Keep only children #
meth_pbmc_children <- meth_pbmc_children[,as.numeric(colnames(meth_pbmc_children)) > 10000]

### Load in Manual Clock Computation Functions ### ----
source("RewritingManualClockComputations.R")

### Define Path to PC Clocks Directory ### ----
### This needs to be downloaded from the PC Clocks GitHub (https://github.com/MorganLevineLab/PC-Clocks) ###
clocksDir <- "~/PC-Clocks-main/"

### Load in Phenotype Information ### ----
PhenotypeData <- read.csv("PhenotypeInformationReduced.csv")
# Edit one SampleID #
PhenotypeData$SampleID[116] <- 2029215
## Extract tissue-specific ages and sex ##
# Buccal - Adult #
Ages_Buccal_Adult <- PhenotypeData$Age[PhenotypeData$Tissue == "Buccal" & PhenotypeData$Age >= 18]
names(Ages_Buccal_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Buccal" & PhenotypeData$Age >= 18]
Ages_Buccal_Adult <- Ages_Buccal_Adult[order(match(names(Ages_Buccal_Adult), colnames(meth_buccal_adults)))]
Sex_Buccal_Adult <- PhenotypeData$Female[PhenotypeData$Tissue == "Buccal" & PhenotypeData$Age >= 18]
names(Sex_Buccal_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Buccal" & PhenotypeData$Age >= 18]
Sex_Buccal_Adult <- Sex_Buccal_Adult[order(match(names(Sex_Buccal_Adult), colnames(meth_buccal_adults)))]
Phenotypes_Buccal_Adult <- data.frame(SampleID = names(Sex_Buccal_Adult),
                                      Age = Ages_Buccal_Adult,
                                      Female = Sex_Buccal_Adult)
# Buccal - Children #
Ages_Buccal_Children <- PhenotypeData$Age[PhenotypeData$Tissue == "Buccal" & PhenotypeData$Age < 18]
names(Ages_Buccal_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Buccal" & PhenotypeData$Age < 18]
Ages_Buccal_Children <- Ages_Buccal_Children[order(match(names(Ages_Buccal_Children), colnames(meth_buccal_children)))]
Sex_Buccal_Children <- PhenotypeData$Female[PhenotypeData$Tissue == "Buccal" & PhenotypeData$Age < 18]
names(Sex_Buccal_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Buccal" & PhenotypeData$Age < 18]
Sex_Buccal_Children <- Sex_Buccal_Children[order(match(names(Sex_Buccal_Children), colnames(meth_buccal_children)))]
Phenotypes_Buccal_Children <- data.frame(SampleID = names(Sex_Buccal_Children),
                                      Age = Ages_Buccal_Children,
                                      Female = Sex_Buccal_Children)
# Saliva - Adult #
Ages_Saliva_Adult <- PhenotypeData$Age[PhenotypeData$Tissue == "Saliva" & PhenotypeData$Age >= 18]
names(Ages_Saliva_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Saliva" & PhenotypeData$Age >= 18]
Ages_Saliva_Adult <- Ages_Saliva_Adult[order(match(names(Ages_Saliva_Adult), colnames(meth_saliva_adults)))]
Sex_Saliva_Adult <- PhenotypeData$Female[PhenotypeData$Tissue == "Saliva" & PhenotypeData$Age >= 18]
names(Sex_Saliva_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Saliva" & PhenotypeData$Age >= 18]
Sex_Saliva_Adult <- Sex_Saliva_Adult[order(match(names(Sex_Saliva_Adult), colnames(meth_saliva_adults)))]
Phenotypes_Saliva_Adult <- data.frame(SampleID = names(Sex_Saliva_Adult),
                                      Age = Ages_Saliva_Adult,
                                      Female = Sex_Saliva_Adult)
# Saliva - Children #
Ages_Saliva_Children <- PhenotypeData$Age[PhenotypeData$Tissue == "Saliva" & PhenotypeData$Age < 18]
names(Ages_Saliva_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Saliva" & PhenotypeData$Age < 18]
Ages_Saliva_Children <- Ages_Saliva_Children[order(match(names(Ages_Saliva_Children), colnames(meth_saliva_children)))]
Sex_Saliva_Children <- PhenotypeData$Female[PhenotypeData$Tissue == "Saliva" & PhenotypeData$Age < 18]
names(Sex_Saliva_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Saliva" & PhenotypeData$Age < 18]
Sex_Saliva_Children <- Sex_Saliva_Children[order(match(names(Sex_Saliva_Children), colnames(meth_saliva_children)))]
Phenotypes_Saliva_Children <- data.frame(SampleID = names(Sex_Saliva_Children),
                                      Age = Ages_Saliva_Children,
                                      Female = Sex_Saliva_Children)
# DBS - Adult #
Ages_DBS_Adult <- PhenotypeData$Age[PhenotypeData$Tissue == "DBS" & PhenotypeData$Age >= 18]
names(Ages_DBS_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "DBS" & PhenotypeData$Age >= 18]
Ages_DBS_Adult <- Ages_DBS_Adult[order(match(names(Ages_DBS_Adult), colnames(meth_dbs_adults)))]
Sex_DBS_Adult <- PhenotypeData$Female[PhenotypeData$Tissue == "DBS" & PhenotypeData$Age >= 18]
names(Sex_DBS_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "DBS" & PhenotypeData$Age >= 18]
Sex_DBS_Adult <- Sex_DBS_Adult[order(match(names(Sex_DBS_Adult), colnames(meth_dbs_adults)))]
Phenotypes_DBS_Adult <- data.frame(SampleID = names(Sex_DBS_Adult),
                                      Age = Ages_DBS_Adult,
                                      Female = Sex_DBS_Adult)
# DBS - Children #
Ages_DBS_Children <- PhenotypeData$Age[PhenotypeData$Tissue == "DBS" & PhenotypeData$Age < 18]
names(Ages_DBS_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "DBS" & PhenotypeData$Age < 18]
Ages_DBS_Children <- Ages_DBS_Children[order(match(names(Ages_DBS_Children), colnames(meth_dbs_children)))]
Sex_DBS_Children <- PhenotypeData$Female[PhenotypeData$Tissue == "DBS" & PhenotypeData$Age < 18]
names(Sex_DBS_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "DBS" & PhenotypeData$Age < 18]
Sex_DBS_Children <- Sex_DBS_Children[order(match(names(Sex_DBS_Children), colnames(meth_dbs_children)))]
Phenotypes_DBS_Children <- data.frame(SampleID = names(Sex_DBS_Children),
                                      Age = Ages_DBS_Children,
                                      Female = Sex_DBS_Children)
# BuffyCoat - Adult #
Ages_BuffyCoat_Adult <- PhenotypeData$Age[PhenotypeData$Tissue == "Buffy Coat" & PhenotypeData$Age >= 18]
names(Ages_BuffyCoat_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Buffy Coat" & PhenotypeData$Age >= 18]
Ages_BuffyCoat_Adult <- Ages_BuffyCoat_Adult[order(match(names(Ages_BuffyCoat_Adult), colnames(meth_buffycoat_adults)))]
Sex_BuffyCoat_Adult <- PhenotypeData$Female[PhenotypeData$Tissue == "Buffy Coat" & PhenotypeData$Age >= 18]
names(Sex_BuffyCoat_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Buffy Coat" & PhenotypeData$Age >= 18]
Sex_BuffyCoat_Adult <- Sex_BuffyCoat_Adult[order(match(names(Sex_BuffyCoat_Adult), colnames(meth_buffycoat_adults)))]
Phenotypes_BuffyCoat_Adult <- data.frame(SampleID = names(Sex_BuffyCoat_Adult),
                                      Age = Ages_BuffyCoat_Adult,
                                      Female = Sex_BuffyCoat_Adult)
# BuffyCoat - Children #
Ages_BuffyCoat_Children <- PhenotypeData$Age[PhenotypeData$Tissue == "Buffy Coat" & PhenotypeData$Age < 18]
names(Ages_BuffyCoat_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Buffy Coat" & PhenotypeData$Age < 18]
Ages_BuffyCoat_Children <- Ages_BuffyCoat_Children[order(match(names(Ages_BuffyCoat_Children), colnames(meth_buffycoat_children)))]
Sex_BuffyCoat_Children <- PhenotypeData$Female[PhenotypeData$Tissue == "Buffy Coat" & PhenotypeData$Age < 18]
names(Sex_BuffyCoat_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "Buffy Coat" & PhenotypeData$Age < 18]
Sex_BuffyCoat_Children <- Sex_BuffyCoat_Children[order(match(names(Sex_BuffyCoat_Children), colnames(meth_buffycoat_children)))]
Phenotypes_BuffyCoat_Children <- data.frame(SampleID = names(Sex_BuffyCoat_Children),
                                         Age = Ages_BuffyCoat_Children,
                                         Female = Sex_BuffyCoat_Children)
# PBMC - Adult #
Ages_PBMC_Adult <- PhenotypeData$Age[PhenotypeData$Tissue == "PBMC" & PhenotypeData$Age >= 18]
names(Ages_PBMC_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "PBMC" & PhenotypeData$Age >= 18]
Ages_PBMC_Adult <- Ages_PBMC_Adult[order(match(names(Ages_PBMC_Adult), colnames(meth_pbmc_adults)))]
Sex_PBMC_Adult <- PhenotypeData$Female[PhenotypeData$Tissue == "PBMC" & PhenotypeData$Age >= 18]
names(Sex_PBMC_Adult) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "PBMC" & PhenotypeData$Age >= 18]
Sex_PBMC_Adult <- Sex_PBMC_Adult[order(match(names(Sex_PBMC_Adult), colnames(meth_pbmc_adults)))]
Phenotypes_PBMC_Adult <- data.frame(SampleID = names(Sex_PBMC_Adult),
                                         Age = Ages_PBMC_Adult,
                                         Female = Sex_PBMC_Adult)
# PBMC - Children #
Ages_PBMC_Children <- PhenotypeData$Age[PhenotypeData$Tissue == "PBMC" & PhenotypeData$Age < 18]
names(Ages_PBMC_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "PBMC" & PhenotypeData$Age < 18]
Ages_PBMC_Children <- Ages_PBMC_Children[order(match(names(Ages_PBMC_Children), colnames(meth_pbmc_children)))]
Sex_PBMC_Children <- PhenotypeData$Female[PhenotypeData$Tissue == "PBMC" & PhenotypeData$Age < 18]
names(Sex_PBMC_Children) <- PhenotypeData$SampleID[PhenotypeData$Tissue == "PBMC" & PhenotypeData$Age < 18]
Sex_PBMC_Children <- Sex_PBMC_Children[order(match(names(Sex_PBMC_Children), colnames(meth_pbmc_children)))]
Phenotypes_PBMC_Children <- data.frame(SampleID = names(Sex_PBMC_Children),
                                    Age = Ages_PBMC_Children,
                                    Female = Sex_PBMC_Children)

### Compute Clocks for Buccal Tissue - Adults ### ----
## Horvath Pan-Tissue Estimations ##
Adult_Buccal_Horvath1 <- calcHorvath1Edit(DNAm = as.data.frame(t(meth_buccal_adults)),
                                          CpGImputation = GoldenStandard_BuccalAdults_Vector,
                                          imputation = TRUE)
## Hannum Estimations ##
Adult_Buccal_Hannum <- calcHannumEdit(DNAm = as.data.frame(t(meth_buccal_adults)),
                                          CpGImputation = GoldenStandard_BuccalAdults_Vector,
                                          imputation = TRUE)
## PhenoAge Estimations ##
Adult_Buccal_PhenoAge <- calcPhenoAgeEdit(DNAm = as.data.frame(t(meth_buccal_adults)),
                                      CpGImputation = GoldenStandard_BuccalAdults_Vector,
                                      imputation = TRUE)
## GrimAge2 Estimations ##
Adult_Buccal_GrimAge2 <- calcGrimAge2Edit(DNAm = as.data.frame(t(meth_buccal_adults)),
                                          CpGImputation = GoldenStandard_BuccalAdults_Vector,
                                          Ages = Ages_Buccal_Adult,
                                          Sex = Sex_Buccal_Adult,
                                          imputation = TRUE)
## DunedinPACE Estimations ##
Adult_Buccal_PACE <- calcPACEEdit(betas = meth_buccal_adults,
                                  proportionOfProbesRequired=0.5,
                                  GoldenStandard = GoldenStandard_BuccalAdults_Vector)
## Horvath2 Estimations ##
Adult_Buccal_Horvath2 <- calcHorvath2Edit(DNAm = as.data.frame(t(meth_buccal_adults)),
                                    CpGImputation = GoldenStandard_BuccalAdults_Vector,
                                    imputation = TRUE)
## PedBE Estimations ##
Adult_Buccal_PedBE <- calcPEDBEEdit(DNAm = as.data.frame(t(meth_buccal_adults)),
                                          CpGImputation = GoldenStandard_BuccalAdults_Vector,
                                          imputation = TRUE)
## PC Clock Estimations ##
Adult_Buccal_PCClock <- calcPCClocksEdit(path_to_PCClocks_directory = clocksDir,
                                         datMeth = as.data.frame(t(meth_buccal_adults)),
                                         datPheno = Phenotypes_Buccal_Adult,
                                         GoldenStandard = GoldenStandard_BuccalAdults_Vector)
## DNAmTL Estimations ##
Adult_Buccal_DNAmTL <- calcDNAmTLEdit(DNAm = as.data.frame(t(meth_buccal_adults)),
                                      CpGImputation = GoldenStandard_BuccalAdults_Vector,
                                      imputation = TRUE)

## Combine Clock Estimates ##
ClockEstimates_Adult_Buccal <- data.frame(SampleID = names(meth_buccal_adults),
                                          mcc_Horvath1 = Adult_Buccal_Horvath1,
                                          mcc_Hannum = Adult_Buccal_Hannum,
                                          mcc_PhenoAge = Adult_Buccal_PhenoAge,
                                          mcc_GrimAge2 = Adult_Buccal_GrimAge2,
                                          mcc_DunedinPACE = Adult_Buccal_PACE$DunedinPACE,
                                          mcc_Horvath2 = Adult_Buccal_Horvath2,
                                          mcc_PedBE = Adult_Buccal_PedBE,
                                          mcc_PCHorvath1 = Adult_Buccal_PCClock$PCHorvath1,
                                          mcc_PCHorvath2 = Adult_Buccal_PCClock$PCHorvath2,
                                          mcc_PCHannum = Adult_Buccal_PCClock$PCHannum,
                                          mcc_PCPhenoAge = Adult_Buccal_PCClock$PCPhenoAge,
                                          mcc_PCGrimAge = Adult_Buccal_PCClock$PCGrimAge,
                                          mcc_PCDNAmTL = Adult_Buccal_PCClock$PCDNAmTL,
                                          mcc_DNAmTL = Adult_Buccal_DNAmTL)

### Compute Clocks for Buccal Tissue - Children ### ----
## Horvath Pan-Tissue Estimations ##
Children_Buccal_Horvath1 <- calcHorvath1Edit(DNAm = as.data.frame(t(meth_buccal_children)),
                                          CpGImputation = GoldenStandard_BuccalChildren_Vector,
                                          imputation = TRUE)
## Hannum Estimations ##
Children_Buccal_Hannum <- calcHannumEdit(DNAm = as.data.frame(t(meth_buccal_children)),
                                      CpGImputation = GoldenStandard_BuccalChildren_Vector,
                                      imputation = TRUE)
## PhenoAge Estimations ##
Children_Buccal_PhenoAge <- calcPhenoAgeEdit(DNAm = as.data.frame(t(meth_buccal_children)),
                                          CpGImputation = GoldenStandard_BuccalChildren_Vector,
                                          imputation = TRUE)
## GrimAge2 Estimations ##
Children_Buccal_GrimAge2 <- calcGrimAge2Edit(DNAm = as.data.frame(t(meth_buccal_children)),
                                          CpGImputation = GoldenStandard_BuccalChildren_Vector,
                                          Ages = Ages_Buccal_Children,
                                          Sex = Sex_Buccal_Children,
                                          imputation = TRUE)
## DunedinPACE Estimations ##
Children_Buccal_PACE <- calcPACEEdit(betas = meth_buccal_children,
                                  proportionOfProbesRequired=0.5,
                                  GoldenStandard = GoldenStandard_BuccalChildren_Vector)
## Horvath2 Estimations ##
Children_Buccal_Horvath2 <- calcHorvath2Edit(DNAm = as.data.frame(t(meth_buccal_children)),
                                          CpGImputation = GoldenStandard_BuccalChildren_Vector,
                                          imputation = TRUE)
## PedBE Estimations ##
Children_Buccal_PedBE <- calcPEDBEEdit(DNAm = as.data.frame(t(meth_buccal_children)),
                                    CpGImputation = GoldenStandard_BuccalChildren_Vector,
                                    imputation = TRUE)
## PC Clock Estimations ##
Children_Buccal_PCClock <- calcPCClocksEdit(path_to_PCClocks_directory = clocksDir,
                                         datMeth = as.data.frame(t(meth_buccal_children)),
                                         datPheno = Phenotypes_Buccal_Children,
                                         GoldenStandard = GoldenStandard_BuccalChildren_Vector)
## DNAmTL Estimations ##
Children_Buccal_DNAmTL <- calcDNAmTLEdit(DNAm = as.data.frame(t(meth_buccal_children)),
                                      CpGImputation = GoldenStandard_BuccalChildren_Vector,
                                      imputation = TRUE)

## Combine Clock Estimates ##
ClockEstimates_Children_Buccal <- data.frame(SampleID = names(meth_buccal_children),
                                          mcc_Horvath1 = Children_Buccal_Horvath1,
                                          mcc_Hannum = Children_Buccal_Hannum,
                                          mcc_PhenoAge = Children_Buccal_PhenoAge,
                                          mcc_GrimAge2 = Children_Buccal_GrimAge2,
                                          mcc_DunedinPACE = Children_Buccal_PACE$DunedinPACE,
                                          mcc_Horvath2 = Children_Buccal_Horvath2,
                                          mcc_PedBE = Children_Buccal_PedBE,
                                          mcc_PCHorvath1 = Children_Buccal_PCClock$PCHorvath1,
                                          mcc_PCHorvath2 = Children_Buccal_PCClock$PCHorvath2,
                                          mcc_PCHannum = Children_Buccal_PCClock$PCHannum,
                                          mcc_PCPhenoAge = Children_Buccal_PCClock$PCPhenoAge,
                                          mcc_PCGrimAge = Children_Buccal_PCClock$PCGrimAge,
                                          mcc_PCDNAmTL = Children_Buccal_PCClock$PCDNAmTL,
                                          mcc_DNAmTL = Children_Buccal_DNAmTL)

### Compute Clocks for Saliva Tissue - Adults ### ----
## Horvath Pan-Tissue Estimations ##
Adult_Saliva_Horvath1 <- calcHorvath1Edit(DNAm = as.data.frame(t(meth_saliva_adults)),
                                          CpGImputation = GoldenStandard_SalivaAdults_Vector,
                                          imputation = TRUE)
## Hannum Estimations ##
Adult_Saliva_Hannum <- calcHannumEdit(DNAm = as.data.frame(t(meth_saliva_adults)),
                                      CpGImputation = GoldenStandard_SalivaAdults_Vector,
                                      imputation = TRUE)
## PhenoAge Estimations ##
Adult_Saliva_PhenoAge <- calcPhenoAgeEdit(DNAm = as.data.frame(t(meth_saliva_adults)),
                                          CpGImputation = GoldenStandard_SalivaAdults_Vector,
                                          imputation = TRUE)
## GrimAge2 Estimations ##
Adult_Saliva_GrimAge2 <- calcGrimAge2Edit(DNAm = as.data.frame(t(meth_saliva_adults)),
                                          CpGImputation = GoldenStandard_SalivaAdults_Vector,
                                          Ages = Ages_Saliva_Adult,
                                          Sex = Sex_Saliva_Adult,
                                          imputation = TRUE)
## DunedinPACE Estimations ##
Adult_Saliva_PACE <- calcPACEEdit(betas = meth_saliva_adults,
                                  proportionOfProbesRequired=0.5,
                                  GoldenStandard = GoldenStandard_SalivaAdults_Vector)
## Horvath2 Estimations ##
Adult_Saliva_Horvath2 <- calcHorvath2Edit(DNAm = as.data.frame(t(meth_saliva_adults)),
                                          CpGImputation = GoldenStandard_SalivaAdults_Vector,
                                          imputation = TRUE)
## PedBE Estimations ##
Adult_Saliva_PedBE <- calcPEDBEEdit(DNAm = as.data.frame(t(meth_saliva_adults)),
                                    CpGImputation = GoldenStandard_SalivaAdults_Vector,
                                    imputation = TRUE)
## PC Clock Estimations ##
Adult_Saliva_PCClock <- calcPCClocksEdit(path_to_PCClocks_directory = clocksDir,
                                         datMeth = as.data.frame(t(meth_saliva_adults)),
                                         datPheno = Phenotypes_Saliva_Adult,
                                         GoldenStandard = GoldenStandard_SalivaAdults_Vector)
## DNAmTL Estimations ##
Adult_Saliva_DNAmTL <- calcDNAmTLEdit(DNAm = as.data.frame(t(meth_saliva_adults)),
                                      CpGImputation = GoldenStandard_SalivaAdults_Vector,
                                      imputation = TRUE)

## Combine Clock Estimates ##
ClockEstimates_Adult_Saliva <- data.frame(SampleID = names(meth_saliva_adults),
                                          mcc_Horvath1 = Adult_Saliva_Horvath1,
                                          mcc_Hannum = Adult_Saliva_Hannum,
                                          mcc_PhenoAge = Adult_Saliva_PhenoAge,
                                          mcc_GrimAge2 = Adult_Saliva_GrimAge2,
                                          mcc_DunedinPACE = Adult_Saliva_PACE$DunedinPACE,
                                          mcc_Horvath2 = Adult_Saliva_Horvath2,
                                          mcc_PedBE = Adult_Saliva_PedBE,
                                          mcc_PCHorvath1 = Adult_Saliva_PCClock$PCHorvath1,
                                          mcc_PCHorvath2 = Adult_Saliva_PCClock$PCHorvath2,
                                          mcc_PCHannum = Adult_Saliva_PCClock$PCHannum,
                                          mcc_PCPhenoAge = Adult_Saliva_PCClock$PCPhenoAge,
                                          mcc_PCGrimAge = Adult_Saliva_PCClock$PCGrimAge,
                                          mcc_PCDNAmTL = Adult_Saliva_PCClock$PCDNAmTL,
                                          mcc_DNAmTL = Adult_Saliva_DNAmTL)

### Compute Clocks for Saliva Tissue - Children ### ----
## Horvath Pan-Tissue Estimations ##
Children_Saliva_Horvath1 <- calcHorvath1Edit(DNAm = as.data.frame(t(meth_saliva_children)),
                                             CpGImputation = GoldenStandard_SalivaChildren_Vector,
                                             imputation = TRUE)
## Hannum Estimations ##
Children_Saliva_Hannum <- calcHannumEdit(DNAm = as.data.frame(t(meth_saliva_children)),
                                         CpGImputation = GoldenStandard_SalivaChildren_Vector,
                                         imputation = TRUE)
## PhenoAge Estimations ##
Children_Saliva_PhenoAge <- calcPhenoAgeEdit(DNAm = as.data.frame(t(meth_saliva_children)),
                                             CpGImputation = GoldenStandard_SalivaChildren_Vector,
                                             imputation = TRUE)
## GrimAge2 Estimations ##
Children_Saliva_GrimAge2 <- calcGrimAge2Edit(DNAm = as.data.frame(t(meth_saliva_children)),
                                             CpGImputation = GoldenStandard_SalivaChildren_Vector,
                                             Ages = Ages_Saliva_Children,
                                             Sex = Sex_Saliva_Children,
                                             imputation = TRUE)
## DunedinPACE Estimations ##
Children_Saliva_PACE <- calcPACEEdit(betas = meth_saliva_children,
                                     proportionOfProbesRequired=0.5,
                                     GoldenStandard = GoldenStandard_SalivaChildren_Vector)
## Horvath2 Estimations ##
Children_Saliva_Horvath2 <- calcHorvath2Edit(DNAm = as.data.frame(t(meth_saliva_children)),
                                             CpGImputation = GoldenStandard_SalivaChildren_Vector,
                                             imputation = TRUE)
## PedBE Estimations ##
Children_Saliva_PedBE <- calcPEDBEEdit(DNAm = as.data.frame(t(meth_saliva_children)),
                                       CpGImputation = GoldenStandard_SalivaChildren_Vector,
                                       imputation = TRUE)
## PC Clock Estimations ##
Children_Saliva_PCClock <- calcPCClocksEdit(path_to_PCClocks_directory = clocksDir,
                                            datMeth = as.data.frame(t(meth_saliva_children)),
                                            datPheno = Phenotypes_Saliva_Children,
                                            GoldenStandard = GoldenStandard_SalivaChildren_Vector)
## DNAmTL Estimations ##
Children_Saliva_DNAmTL <- calcDNAmTLEdit(DNAm = as.data.frame(t(meth_saliva_children)),
                                         CpGImputation = GoldenStandard_SalivaChildren_Vector,
                                         imputation = TRUE)

## Combine Clock Estimates ##
ClockEstimates_Children_Saliva <- data.frame(SampleID = names(meth_saliva_children),
                                             mcc_Horvath1 = Children_Saliva_Horvath1,
                                             mcc_Hannum = Children_Saliva_Hannum,
                                             mcc_PhenoAge = Children_Saliva_PhenoAge,
                                             mcc_GrimAge2 = Children_Saliva_GrimAge2,
                                             mcc_DunedinPACE = Children_Saliva_PACE$DunedinPACE,
                                             mcc_Horvath2 = Children_Saliva_Horvath2,
                                             mcc_PedBE = Children_Saliva_PedBE,
                                             mcc_PCHorvath1 = Children_Saliva_PCClock$PCHorvath1,
                                             mcc_PCHorvath2 = Children_Saliva_PCClock$PCHorvath2,
                                             mcc_PCHannum = Children_Saliva_PCClock$PCHannum,
                                             mcc_PCPhenoAge = Children_Saliva_PCClock$PCPhenoAge,
                                             mcc_PCGrimAge = Children_Saliva_PCClock$PCGrimAge,
                                             mcc_PCDNAmTL = Children_Saliva_PCClock$PCDNAmTL,
                                             mcc_DNAmTL = Children_Saliva_DNAmTL)

### Compute Clocks for DBS Tissue - Adults ### ----
## Horvath Pan-Tissue Estimations ##
Adult_DBS_Horvath1 <- calcHorvath1Edit(DNAm = as.data.frame(t(meth_dbs_adults)),
                                          CpGImputation = GoldenStandard_WholeBloodAdults_Vector,
                                          imputation = TRUE)
## Hannum Estimations ##
Adult_DBS_Hannum <- calcHannumEdit(DNAm = as.data.frame(t(meth_dbs_adults)),
                                      CpGImputation = GoldenStandard_WholeBloodAdults_Vector,
                                      imputation = TRUE)
## PhenoAge Estimations ##
Adult_DBS_PhenoAge <- calcPhenoAgeEdit(DNAm = as.data.frame(t(meth_dbs_adults)),
                                          CpGImputation = GoldenStandard_WholeBloodAdults_Vector,
                                          imputation = TRUE)
## GrimAge2 Estimations ##
Adult_DBS_GrimAge2 <- calcGrimAge2Edit(DNAm = as.data.frame(t(meth_dbs_adults)),
                                          CpGImputation = GoldenStandard_WholeBloodAdults_Vector,
                                          Ages = Ages_DBS_Adult,
                                          Sex = Sex_DBS_Adult,
                                          imputation = TRUE)
## DunedinPACE Estimations ##
Adult_DBS_PACE <- calcPACEEdit(betas = meth_dbs_adults,
                                  proportionOfProbesRequired=0.5,
                                  GoldenStandard = GoldenStandard_WholeBloodAdults_Vector)
## Horvath2 Estimations ##
Adult_DBS_Horvath2 <- calcHorvath2Edit(DNAm = as.data.frame(t(meth_dbs_adults)),
                                          CpGImputation = GoldenStandard_WholeBloodAdults_Vector,
                                          imputation = TRUE)
## PedBE Estimations ##
Adult_DBS_PedBE <- calcPEDBEEdit(DNAm = as.data.frame(t(meth_dbs_adults)),
                                    CpGImputation = GoldenStandard_WholeBloodAdults_Vector,
                                    imputation = TRUE)
## PC Clock Estimations ##
Adult_DBS_PCClock <- calcPCClocksEdit(path_to_PCClocks_directory = clocksDir,
                                         datMeth = as.data.frame(t(meth_dbs_adults)),
                                         datPheno = Phenotypes_DBS_Adult,
                                         GoldenStandard = GoldenStandard_WholeBloodAdults_Vector)
## DNAmTL Estimations ##
Adult_DBS_DNAmTL <- calcDNAmTLEdit(DNAm = as.data.frame(t(meth_dbs_adults)),
                                      CpGImputation = GoldenStandard_WholeBloodAdults_Vector,
                                      imputation = TRUE)

## Combine Clock Estimates ##
ClockEstimates_Adult_DBS <- data.frame(SampleID = names(meth_dbs_adults),
                                          mcc_Horvath1 = Adult_DBS_Horvath1,
                                          mcc_Hannum = Adult_DBS_Hannum,
                                          mcc_PhenoAge = Adult_DBS_PhenoAge,
                                          mcc_GrimAge2 = Adult_DBS_GrimAge2,
                                          mcc_DunedinPACE = Adult_DBS_PACE$DunedinPACE,
                                          mcc_Horvath2 = Adult_DBS_Horvath2,
                                          mcc_PedBE = Adult_DBS_PedBE,
                                          mcc_PCHorvath1 = Adult_DBS_PCClock$PCHorvath1,
                                          mcc_PCHorvath2 = Adult_DBS_PCClock$PCHorvath2,
                                          mcc_PCHannum = Adult_DBS_PCClock$PCHannum,
                                          mcc_PCPhenoAge = Adult_DBS_PCClock$PCPhenoAge,
                                          mcc_PCGrimAge = Adult_DBS_PCClock$PCGrimAge,
                                          mcc_PCDNAmTL = Adult_DBS_PCClock$PCDNAmTL,
                                          mcc_DNAmTL = Adult_DBS_DNAmTL)

### Compute Clocks for DBS Tissue - Children ### ----
## Horvath Pan-Tissue Estimations ##
Children_DBS_Horvath1 <- calcHorvath1Edit(DNAm = as.data.frame(t(meth_dbs_children)),
                                             CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                             imputation = TRUE)
## Hannum Estimations ##
Children_DBS_Hannum <- calcHannumEdit(DNAm = as.data.frame(t(meth_dbs_children)),
                                         CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                         imputation = TRUE)
## PhenoAge Estimations ##
Children_DBS_PhenoAge <- calcPhenoAgeEdit(DNAm = as.data.frame(t(meth_dbs_children)),
                                             CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                             imputation = TRUE)
## GrimAge2 Estimations ##
Children_DBS_GrimAge2 <- calcGrimAge2Edit(DNAm = as.data.frame(t(meth_dbs_children)),
                                             CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                             Ages = Ages_DBS_Children,
                                             Sex = Sex_DBS_Children,
                                             imputation = TRUE)
## DunedinPACE Estimations ##
Children_DBS_PACE <- calcPACEEdit(betas = meth_dbs_children,
                                     proportionOfProbesRequired=0.5,
                                     GoldenStandard = GoldenStandard_WholeBloodChildren_Vector)
## Horvath2 Estimations ##
Children_DBS_Horvath2 <- calcHorvath2Edit(DNAm = as.data.frame(t(meth_dbs_children)),
                                             CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                             imputation = TRUE)
## PedBE Estimations ##
Children_DBS_PedBE <- calcPEDBEEdit(DNAm = as.data.frame(t(meth_dbs_children)),
                                       CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                       imputation = TRUE)
## PC Clock Estimations ##
Children_DBS_PCClock <- calcPCClocksEdit(path_to_PCClocks_directory = clocksDir,
                                            datMeth = as.data.frame(t(meth_dbs_children)),
                                            datPheno = Phenotypes_DBS_Children,
                                            GoldenStandard = GoldenStandard_WholeBloodChildren_Vector)
## DNAmTL Estimations ##
Children_DBS_DNAmTL <- calcDNAmTLEdit(DNAm = as.data.frame(t(meth_dbs_children)),
                                         CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                         imputation = TRUE)

## Combine Clock Estimates ##
ClockEstimates_Children_DBS <- data.frame(SampleID = names(meth_dbs_children),
                                             mcc_Horvath1 = Children_DBS_Horvath1,
                                             mcc_Hannum = Children_DBS_Hannum,
                                             mcc_PhenoAge = Children_DBS_PhenoAge,
                                             mcc_GrimAge2 = Children_DBS_GrimAge2,
                                             mcc_DunedinPACE = Children_DBS_PACE$DunedinPACE,
                                             mcc_Horvath2 = Children_DBS_Horvath2,
                                             mcc_PedBE = Children_DBS_PedBE,
                                             mcc_PCHorvath1 = Children_DBS_PCClock$PCHorvath1,
                                             mcc_PCHorvath2 = Children_DBS_PCClock$PCHorvath2,
                                             mcc_PCHannum = Children_DBS_PCClock$PCHannum,
                                             mcc_PCPhenoAge = Children_DBS_PCClock$PCPhenoAge,
                                             mcc_PCGrimAge = Children_DBS_PCClock$PCGrimAge,
                                             mcc_PCDNAmTL = Children_DBS_PCClock$PCDNAmTL,
                                             mcc_DNAmTL = Children_DBS_DNAmTL)

### Compute Clocks for BuffyCoat Tissue - Children ### ----
## Horvath Pan-Tissue Estimations ##
Children_BuffyCoat_Horvath1 <- calcHorvath1Edit(DNAm = as.data.frame(t(meth_buffycoat_children)),
                                          CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                          imputation = TRUE)
## Hannum Estimations ##
Children_BuffyCoat_Hannum <- calcHannumEdit(DNAm = as.data.frame(t(meth_buffycoat_children)),
                                      CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                      imputation = TRUE)
## PhenoAge Estimations ##
Children_BuffyCoat_PhenoAge <- calcPhenoAgeEdit(DNAm = as.data.frame(t(meth_buffycoat_children)),
                                          CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                          imputation = TRUE)
## GrimAge2 Estimations ##
Children_BuffyCoat_GrimAge2 <- calcGrimAge2Edit(DNAm = as.data.frame(t(meth_buffycoat_children)),
                                          CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                          Ages = Ages_BuffyCoat_Children,
                                          Sex = Sex_BuffyCoat_Children,
                                          imputation = TRUE)
## DunedinPACE Estimations ##
Children_BuffyCoat_PACE <- calcPACEEdit(betas = meth_buffycoat_children,
                                  proportionOfProbesRequired=0.5,
                                  GoldenStandard = GoldenStandard_WholeBloodChildren_Vector)
## Horvath2 Estimations ##
Children_BuffyCoat_Horvath2 <- calcHorvath2Edit(DNAm = as.data.frame(t(meth_buffycoat_children)),
                                          CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                          imputation = TRUE)
## PedBE Estimations ##
Children_BuffyCoat_PedBE <- calcPEDBEEdit(DNAm = as.data.frame(t(meth_buffycoat_children)),
                                    CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                    imputation = TRUE)
## PC Clock Estimations ##
Children_BuffyCoat_PCClock <- calcPCClocksEdit(path_to_PCClocks_directory = clocksDir,
                                         datMeth = as.data.frame(t(meth_buffycoat_children)),
                                         datPheno = Phenotypes_BuffyCoat_Children,
                                         GoldenStandard = GoldenStandard_WholeBloodChildren_Vector)
## DNAmTL Estimations ##
Children_BuffyCoat_DNAmTL <- calcDNAmTLEdit(DNAm = as.data.frame(t(meth_buffycoat_children)),
                                      CpGImputation = GoldenStandard_WholeBloodChildren_Vector,
                                      imputation = TRUE)

## Combine Clock Estimates ##
ClockEstimates_Children_BuffyCoat <- data.frame(SampleID = names(meth_buffycoat_children),
                                          mcc_Horvath1 = Children_BuffyCoat_Horvath1,
                                          mcc_Hannum = Children_BuffyCoat_Hannum,
                                          mcc_PhenoAge = Children_BuffyCoat_PhenoAge,
                                          mcc_GrimAge2 = Children_BuffyCoat_GrimAge2,
                                          mcc_DunedinPACE = Children_BuffyCoat_PACE$DunedinPACE,
                                          mcc_Horvath2 = Children_BuffyCoat_Horvath2,
                                          mcc_PedBE = Children_BuffyCoat_PedBE,
                                          mcc_PCHorvath1 = Children_BuffyCoat_PCClock$PCHorvath1,
                                          mcc_PCHorvath2 = Children_BuffyCoat_PCClock$PCHorvath2,
                                          mcc_PCHannum = Children_BuffyCoat_PCClock$PCHannum,
                                          mcc_PCPhenoAge = Children_BuffyCoat_PCClock$PCPhenoAge,
                                          mcc_PCGrimAge = Children_BuffyCoat_PCClock$PCGrimAge,
                                          mcc_PCDNAmTL = Children_BuffyCoat_PCClock$PCDNAmTL,
                                          mcc_DNAmTL = Children_BuffyCoat_DNAmTL)

### Compute Clocks for PBMC Tissue - Adults ### ----
## Horvath Pan-Tissue Estimations ##
Adult_PBMC_Horvath1 <- calcHorvath1Edit(DNAm = as.data.frame(t(meth_pbmc_adults)),
                                       CpGImputation = GoldenStandard_PBMCAdults_Vector,
                                       imputation = TRUE)
## Hannum Estimations ##
Adult_PBMC_Hannum <- calcHannumEdit(DNAm = as.data.frame(t(meth_pbmc_adults)),
                                   CpGImputation = GoldenStandard_PBMCAdults_Vector,
                                   imputation = TRUE)
## PhenoAge Estimations ##
Adult_PBMC_PhenoAge <- calcPhenoAgeEdit(DNAm = as.data.frame(t(meth_pbmc_adults)),
                                       CpGImputation = GoldenStandard_PBMCAdults_Vector,
                                       imputation = TRUE)
## GrimAge2 Estimations ##
Adult_PBMC_GrimAge2 <- calcGrimAge2Edit(DNAm = as.data.frame(t(meth_pbmc_adults)),
                                       CpGImputation = GoldenStandard_PBMCAdults_Vector,
                                       Ages = Ages_PBMC_Adult,
                                       Sex = Sex_PBMC_Adult,
                                       imputation = TRUE)
## DunedinPACE Estimations ##
Adult_PBMC_PACE <- calcPACEEdit(betas = meth_pbmc_adults,
                               proportionOfProbesRequired=0.5,
                               GoldenStandard = GoldenStandard_PBMCAdults_Vector)
## Horvath2 Estimations ##
Adult_PBMC_Horvath2 <- calcHorvath2Edit(DNAm = as.data.frame(t(meth_pbmc_adults)),
                                       CpGImputation = GoldenStandard_PBMCAdults_Vector,
                                       imputation = TRUE)
## PedBE Estimations ##
Adult_PBMC_PedBE <- calcPEDBEEdit(DNAm = as.data.frame(t(meth_pbmc_adults)),
                                 CpGImputation = GoldenStandard_PBMCAdults_Vector,
                                 imputation = TRUE)
## PC Clock Estimations ##
Adult_PBMC_PCClock <- calcPCClocksEdit(path_to_PCClocks_directory = clocksDir,
                                      datMeth = as.data.frame(t(meth_pbmc_adults)),
                                      datPheno = Phenotypes_PBMC_Adult,
                                      GoldenStandard = GoldenStandard_PBMCAdults_Vector)
## DNAmTL Estimations ##
Adult_PBMC_DNAmTL <- calcDNAmTLEdit(DNAm = as.data.frame(t(meth_pbmc_adults)),
                                   CpGImputation = GoldenStandard_PBMCAdults_Vector,
                                   imputation = TRUE)

## Combine Clock Estimates ##
ClockEstimates_Adult_PBMC <- data.frame(SampleID = names(meth_pbmc_adults),
                                       mcc_Horvath1 = Adult_PBMC_Horvath1,
                                       mcc_Hannum = Adult_PBMC_Hannum,
                                       mcc_PhenoAge = Adult_PBMC_PhenoAge,
                                       mcc_GrimAge2 = Adult_PBMC_GrimAge2,
                                       mcc_DunedinPACE = Adult_PBMC_PACE$DunedinPACE,
                                       mcc_Horvath2 = Adult_PBMC_Horvath2,
                                       mcc_PedBE = Adult_PBMC_PedBE,
                                       mcc_PCHorvath1 = Adult_PBMC_PCClock$PCHorvath1,
                                       mcc_PCHorvath2 = Adult_PBMC_PCClock$PCHorvath2,
                                       mcc_PCHannum = Adult_PBMC_PCClock$PCHannum,
                                       mcc_PCPhenoAge = Adult_PBMC_PCClock$PCPhenoAge,
                                       mcc_PCGrimAge = Adult_PBMC_PCClock$PCGrimAge,
                                       mcc_PCDNAmTL = Adult_PBMC_PCClock$PCDNAmTL,
                                       mcc_DNAmTL = Adult_PBMC_DNAmTL)

### Combine All Data ### ----
# Combine all clock data #
ClockEstimates_All <- rbind(ClockEstimates_Adult_Buccal,
                            ClockEstimates_Children_Buccal,
                            ClockEstimates_Adult_Saliva,
                            ClockEstimates_Children_Saliva,
                            ClockEstimates_Adult_DBS,
                            ClockEstimates_Children_DBS,
                            ClockEstimates_Adult_PBMC,
                            ClockEstimates_Children_BuffyCoat)
# Merge clock data with phenotypes #
data_all <- merge(PhenotypeData, ClockEstimates_All, by = "SampleID")
data_all <- data_all[,-c(2)]

### Calculate Age Acceleration Values (Subtraction Method) ### ----
data_all$mcc_Horvath1_AccMinus <- data_all$mcc_Horvath1 - data_all$Age
data_all$mcc_Hannum_AccMinus <- data_all$mcc_Hannum - data_all$Age
data_all$mcc_PhenoAge_AccMinus <- data_all$mcc_PhenoAge - data_all$Age
data_all$mcc_GrimAge2_AccMinus <- data_all$mcc_GrimAge2 - data_all$Age
data_all$mcc_Horvath2_AccMinus <- data_all$mcc_Horvath2 - data_all$Age
data_all$mcc_PedBE_AccMinus <- data_all$mcc_PedBE - data_all$Age
data_all$mcc_PCHorvath1_AccMinus <- data_all$mcc_PCHorvath1 - data_all$Age
data_all$mcc_PCHorvath2_AccMinus <- data_all$mcc_PCHorvath2 - data_all$Age
data_all$mcc_PCHannum_AccMinus <- data_all$mcc_PCHannum - data_all$Age
data_all$mcc_PCPhenoAge_AccMinus <- data_all$mcc_PCPhenoAge - data_all$Age
data_all$mcc_PCGrimAge_AccMinus <- data_all$mcc_PCGrimAge - data_all$Age

### Write as File ### ----
write.csv(data_all, "ManualClockComputationData_v3.csv")





