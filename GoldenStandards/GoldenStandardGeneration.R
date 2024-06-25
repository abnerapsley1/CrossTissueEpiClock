### Publicly Available Data - Golden Standard Dataset Generation ### ----

### Load Packages ### ----
library(minfi)
library(wateRmelon)
library(tibble)
library(dplyr)
library(methylCIPHER)

### Load in Adult Buccal Data ### ----
data_GSE50586 <- read.csv("./PubliclyAvailableData/Buccal_Adults/GSE50586/GSE50586_ClockProbes_v3.csv")[,-1]
data_GSE166844 <- read.csv("./PubliclyAvailableData/Buccal_Adults/GSE166844/GSE166844_ClockProbes_v3.csv")[,-1]

### Load in Children Buccal Data ### ----
data_GSE50759 <- read.csv("./PubliclyAvailableData/Buccal_Children/GSE50759/GSE50759_ClockProbes_v3.csv")[,-1]
data_GSE147058 <- read.csv("./PubliclyAvailableData/Buccal_Children/GSE147058/GSE147058_ClockProbes_v3.csv")[,-1]

### Load in Adult Saliva Data ### ----
data_GSE130153 <- read.csv("./PubliclyAvailableData/Saliva_Adults/GSE130153/GSE130153_ClockProbes_v3.csv")[,-1]
data_GSE111165 <- read.csv("./PubliclyAvailableData/Saliva_Adults/GSE111165/GSE111165_ClockProbes_v3.csv")[,-1]

### Load in Children Saliva Data ### ----
data_GSE147318 <- read.csv("./PubliclyAvailableData/Saliva_Children/GSE147318/GSE147318_ClockProbes_v3.csv")[,-1]
data_GSE112314 <- read.csv("./PubliclyAvailableData/Saliva_Children/GSE112314/GSE112314_ClockProbes_v3.csv")[,-1]

### Load in Adult Whole Blood Data ### ----
data_GSE201287 <- read.csv("./PubliclyAvailableData/WholeBlood_Adults/GSE201287/GSE201287_ClockProbes_v3.csv")[,-1]
data_GSE218186 <- read.csv("./PubliclyAvailableData/WholeBlood_Adults/GSE218186/GSE218186_ClockProbes_v3.csv")[,-1]

### Load in Children Whole Blood Data ### ----
data_GSE174555 <- read.csv("./PubliclyAvailableData/WholeBlood_Children/GSE174555/GSE174555_ClockProbes_v3.csv")[,-1]
data_GSE221864 <- read.csv("./PubliclyAvailableData/WholeBlood_Children/GSE221864/GSE221864_ClockProbes_v3.csv")[,-1]

### Load in Adult PBMC Data ### ----
data_GSE117929 <- read.csv("./PubliclyAvailableData/PBMC_Adults/GSE117929/GSE117929_ClockProbes_v3.csv")[,-1]
data_GSE245924 <- read.csv("./PubliclyAvailableData/PBMC_Adults/GSE245924/GSE245924_ClockProbes_v3.csv")[,-1]

### Load in Children PBMC Data ### ----
data_GSE40576 <- read.csv("./PubliclyAvailableData/PBMC_Children/GSE40576/GSE40576_ClockProbes_v3.csv")[,-1]
data_GSE132181 <- read.csv("./PubliclyAvailableData/PBMC_Children/GSE132181/GSE132181_ClockProbes_v3.csv")[,-1]

### Combine Data ### ----
data_buccal_adults <- cbind(data_GSE50586, data_GSE166844[,-1])
data_buccal_children <- cbind(data_GSE50759, data_GSE147058[,-1])
data_saliva_adults <- cbind(data_GSE130153, data_GSE111165[,-1])
data_saliva_children <- cbind(data_GSE147318, data_GSE112314[,-1])
data_wholeblood_adults <- cbind(data_GSE201287, data_GSE218186[,-1])
data_wholeblood_children <- cbind(data_GSE174555, data_GSE221864[,-1])
data_pbmc_adults <- cbind(data_GSE117929, data_GSE245924[,-1])
data_pbmc_children <- cbind(data_GSE40576, data_GSE132181[,-1])

### Calculate Golden Standard Means - Buccal Adults ### ----
GoldenStandard_BuccalAdults <- data.frame(ProbeID = data_buccal_adults$ProbeID,
                                          MeanBeta = rowMeans(data_buccal_adults[,-1]))
# Create named vector #
GoldenStandard_BuccalAdults_Vector <- as.vector(GoldenStandard_BuccalAdults$MeanBeta)
names(GoldenStandard_BuccalAdults_Vector) <- GoldenStandard_BuccalAdults$ProbeID

### Calculate Golden Standard Means - Buccal Children ### ----
GoldenStandard_BuccalChildren <- data.frame(ProbeID = data_buccal_children$ProbeID,
                                          MeanBeta = rowMeans(data_buccal_children[,-1]))
# Create named vector #
GoldenStandard_BuccalChildren_Vector <- as.vector(GoldenStandard_BuccalChildren$MeanBeta)
names(GoldenStandard_BuccalChildren_Vector) <- GoldenStandard_BuccalChildren$ProbeID

### Calculate Golden Standard Means - Saliva Adults ### ----
GoldenStandard_SalivaAdults <- data.frame(ProbeID = data_saliva_adults$ProbeID,
                                          MeanBeta = rowMeans(data_saliva_adults[,-1]))
# Create named vector #
GoldenStandard_SalivaAdults_Vector <- as.vector(GoldenStandard_SalivaAdults$MeanBeta)
names(GoldenStandard_SalivaAdults_Vector) <- GoldenStandard_SalivaAdults$ProbeID

### Calculate Golden Standard Means - Saliva Children ### ----
GoldenStandard_SalivaChildren <- data.frame(ProbeID = data_saliva_children$ProbeID,
                                            MeanBeta = rowMeans(data_saliva_children[,-1]))
# Create named vector #
GoldenStandard_SalivaChildren_Vector <- as.vector(GoldenStandard_SalivaChildren$MeanBeta)
names(GoldenStandard_SalivaChildren_Vector) <- GoldenStandard_SalivaChildren$ProbeID

### Calculate Golden Standard Means - Whole Blood Adults ### ----
GoldenStandard_WholeBloodAdults <- data.frame(ProbeID = data_wholeblood_adults$ProbeID,
                                          MeanBeta = rowMeans(data_wholeblood_adults[,-1]))
# Create named vector #
GoldenStandard_WholeBloodAdults_Vector <- as.vector(GoldenStandard_WholeBloodAdults$MeanBeta)
names(GoldenStandard_WholeBloodAdults_Vector) <- GoldenStandard_WholeBloodAdults$ProbeID

### Calculate Golden Standard Means - WholeBlood Children ### ----
GoldenStandard_WholeBloodChildren <- data.frame(ProbeID = data_wholeblood_children$ProbeID,
                                            MeanBeta = rowMeans(data_wholeblood_children[,-1]))
# Create named vector #
GoldenStandard_WholeBloodChildren_Vector <- as.vector(GoldenStandard_WholeBloodChildren$MeanBeta)
names(GoldenStandard_WholeBloodChildren_Vector) <- GoldenStandard_WholeBloodChildren$ProbeID

### Calculate Golden Standard Means - PBMC Adults ### ----
GoldenStandard_PBMCAdults <- data.frame(ProbeID = data_pbmc_adults$ProbeID,
                                          MeanBeta = rowMeans(data_pbmc_adults[,-1]))
# Create named vector #
GoldenStandard_PBMCAdults_Vector <- as.vector(GoldenStandard_PBMCAdults$MeanBeta)
names(GoldenStandard_PBMCAdults_Vector) <- GoldenStandard_PBMCAdults$ProbeID

### Calculate Golden Standard Means - PBMC Children ### ----
GoldenStandard_PBMCChildren <- data.frame(ProbeID = data_pbmc_children$ProbeID,
                                            MeanBeta = rowMeans(data_pbmc_children[,-1]))
# Create named vector #
GoldenStandard_PBMCChildren_Vector <- as.vector(GoldenStandard_PBMCChildren$MeanBeta)
names(GoldenStandard_PBMCChildren_Vector) <- GoldenStandard_PBMCChildren$ProbeID

### Save Golden Standard Means as RDS Files ### ----
saveRDS(GoldenStandard_BuccalAdults_Vector, "GoldenStandard_Buccal_Adults_Mean_v3.Rds")
saveRDS(GoldenStandard_BuccalChildren_Vector, "GoldenStandard_Buccal_Children_Mean_v3.Rds")
saveRDS(GoldenStandard_SalivaAdults_Vector, "GoldenStandard_Saliva_Adults_Mean_v3.Rds")
saveRDS(GoldenStandard_SalivaChildren_Vector, "GoldenStandard_Saliva_Children_Mean_v3.Rds")
saveRDS(GoldenStandard_WholeBloodAdults_Vector, "GoldenStandard_WholeBlood_Adults_Mean_v3.Rds")
saveRDS(GoldenStandard_WholeBloodChildren_Vector, "GoldenStandard_WholeBlood_Children_Mean_v3.Rds")
saveRDS(GoldenStandard_PBMCAdults_Vector, "GoldenStandard_PBMC_Adults_Mean_v3.Rds")
saveRDS(GoldenStandard_PBMCChildren_Vector, "GoldenStandard_PBMC_Children_Mean_v3.Rds")

### Delete Unnecessary Environmental Variables ### ----
data_buccal_adults <- NULL
data_buccal_children <- NULL
data_saliva_adults <- NULL
data_saliva_children <- NULL
data_wholeblood_adults <- NULL
data_wholeblood_children <- NULL
data_pbmc_adults <- NULL
data_pbmc_children <- NULL
data_GSE111165 <- NULL
data_GSE112314 <- NULL
data_GSE117929 <- NULL
data_GSE130153 <- NULL
data_GSE132181 <- NULL
data_GSE147058 <- NULL
data_GSE147318 <- NULL
data_GSE166844 <- NULL
data_GSE174555 <- NULL
data_GSE201287 <- NULL
data_GSE218186 <- NULL
data_GSE221864 <- NULL
data_GSE245924 <- NULL
data_GSE40576 <- NULL
data_GSE50586 <- NULL
data_GSE50759 <- NULL
GoldenStandard_BuccalAdults <- NULL
GoldenStandard_BuccalChildren <- NULL
GoldenStandard_PBMCAdults <- NULL
GoldenStandard_PBMCChildren <- NULL
GoldenStandard_SalivaAdults <- NULL
GoldenStandard_SalivaChildren <- NULL
GoldenStandard_WholeBloodAdults <- NULL
GoldenStandard_WholeBloodChildren <- NULL

