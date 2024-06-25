### EPIC Array Preprocessing ### ----
# This code is for converting all idat files from the EPIC arrays into Beta-value matrices #
# The Beta-value matrices are then used for various downstream analyses such as predicting mDNA ages and CBC immune cell compositions #

### Load Packages ### ----
library(knitr)
library(limma)
library(minfi)
# Change for type of array being used #
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Change for type of array being used #
library(IlluminaHumanMethylationEPICv2manifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(ENmix)
library(FlowSorted.Blood.EPIC)
library(FlowSorted.BloodExtended.EPIC)
library(sqldf)
library(tidyverse)
library(DT)
library(plotly)
library(gt)
library(BiocManager)
library(wateRmelon)
library(RPMM)
require(minfiData)
library(tibble)
library(DunedinPoAm38)

### Load in Data ### ----
# Create annotation matrix for EPIC probe sites #
# If you are using a different array (i.e., 450k), this should be changed #
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Edit probe names #
#annEPIC@rownames = substr(annEPIC@rownames,1,nchar(annEPIC@rownames)-5)
# Set directory for data #
# This is the path from you working directory to the directory containing all iDAT and phenotype information #
dataDirectory <- "./RawData_Batch1"
# When reading in sample sheets below, ensure that the SentrixID is not in scientific notation #
# Once the csv file is ready, open it in text-edit and add a line break at the end of the document, otherwise this will not run correctly #
# Read in iDAT file names and phenotypes from EPIC array 1 #
# Pattern should be the name of the sample sheet #
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet_Batch1.csv")

### Create an RGChannel Set for the EPIC Array ### ----
# Define the RGChannel Set #
rgSet <- read.metharray.exp(targets=targets, force=TRUE, extended = TRUE)
rgSet@annotation=c(array='IlluminaHumanMethylationEPICv2', annotation='20a1.hg38')
# Save the rgSet #
#saveRDS(rgSet, "rgSet_raw.Rds")
# Read the rgSet #
#rgSet <- readRDS("rgSet_raw.Rds")
# Calculate detection p-values for each probe and each sample #
detP <- detectionP(rgSet)
# Rename Samples in all rgSets #
# Ensure sample names are in a column named Sample_Name #
sampleNames(rgSet) <- targets$Sample_Name


### Quality Control - Pt. 1 ### ----
## Remove samples with an average detection p-value above 0.05 ##
# Create a keep vector that meets the exclusion criteria #
keep <- colMeans(detP) < 0.05
# Filter RGChannel Set for bad samples #
rgSet <- rgSet[,keep]
# Remove poor quality samples from targets/phenotype data #
targets <- targets[keep,]
# Remove poor quality samples from detection p-value data #
detP <- detP[,keep]
## Check Bisulfite Conversion Rates ##
# Compute bisulfite conversion rates #
rgSet_bsconv <- bscon(rgSet)
# Create keep mask #
keep <- rgSet_bsconv > 0.80
# Filter out samples not passing bisulfite conversion threshold (80%) #
rgSet <- rgSet[,keep]
# Remove poor quality samples from targets/phenotype data #
targets <- targets[keep,]
# Remove poor quality samples from detection p-value data #
detP <- detP[,keep]

### Extract out Beta Values of Control Probes ### ----
# Get manifest for array #
manifest <- getManifest(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
# Define control probe addresses #
ControlProbeAddress <- manifest@data$TypeControl$Address
# Define keep vector for control probes #
keep_control <- rgSet@NAMES %in% ControlProbeAddress
# Define methylated signal matrix for control probes #
ControlMeth <- rgSet@assays@data@listData$Red[keep_control,]
# Define unmethylated signal matrix for control probes #
ControlUnmeth <- rgSet@assays@data@listData$Green[keep_control,]
# Create empty beta matrix #
ControlBeta <- matrix(nrow = nrow(ControlMeth), ncol = ncol(ControlMeth))
# Add sample names to matrix #
colnames(ControlBeta) <- targets$Sample_Name
# Add control probe addresses to matrix #
rownames(ControlBeta) <- rownames(ControlMeth)
# Compute beta values from meth and unmeth signals #
for (i in 1:ncol(ControlMeth)){
  for (j in 1:nrow(ControlMeth)){
    ControlBeta[j,i] <- max(ControlMeth[j,i],0) / (max(ControlMeth[j,i],0) + max(ControlUnmeth[j,i],0) + 100)
  }
}
# Save beta value matrix for further computations #
saveRDS(ControlBeta, "Batch1_ControlBeta.RDS")
# Delete unnecessary files #
manifest <- NULL
ControlProbeAddress <- NULL
keep_control <- NULL
ControlMeth <- NULL
ControlUnmeth <- NULL
ControlBeta <- NULL

### Quality Control - Pt. 2 ### ----
## Check Bead Counts ##
# Create matrix of bead counts #
rgSet_bc <- beadcount(rgSet)
# Create mask to exclude probes with bead counts less than 3 in greater than 5% of samples #
keep <- rowSums(rgSet_bc < 3, na.rm = TRUE) <= ceiling(0.05 * ncol(rgSet_bc))
table(keep)

### Normalization ### ----
# Normalize the data #
mSetSq <- preprocessRaw(rgSet)
mSetSqBMIQ <- BMIQ(mSetSq)

### P-Value Probe Filtering ### ----
## Filter for Bead Count ##
# Remove poor quality probes from detection p-value data #
detP <- detP[keep,]
# Remove poor quality probes from beta values #
mSetSqBMIQ <- mSetSqBMIQ[keep,]
## Filter for probe p-values ##
# Create a keep vector that meets the exclusion criteria #
keep <- rowMeans(detP) < 0.05
write.csv(rownames(detP[!keep,]), "Batch1_pvalFilt_Probes.csv")
# Remove poor quality probes from detection p-value data #
detP <- detP[keep,]
# Remove poor quality probes from beta values #
mSetSqBMIQ <- mSetSqBMIQ[keep,]

### Remove Probe Suffixes and Combine Duplicate Probes ### ----
mSetSqBMIQ <- rm.cgsuffix(mSetSqBMIQ)

### Generate Beta or M Value Matrices for Further Analyses ### ----
# Write Beta-values to RDS #
saveRDS(mSetSqBMIQ, "bVals_SampFilter_ProbeFilt_BMIQNorm_Batch1.RDS")
