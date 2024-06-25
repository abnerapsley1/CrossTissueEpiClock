### Bootstrapping Reference-Free Cell Deconvolution ### ----

### RefFreeEWAS files were downloaded from the depreciated R package RefFreeEWAS_2.2 (https://cran.r-project.org/src/contrib/Archive/RefFreeEWAS/) ###
### Loading in RefFreeEWAS Data ### ----
load("./RefFreeEWAS/data/HNSCC.RData")
load("./RefFreeEWAS/data/RefFreeEWAS.RData")
library(quadprog)

### Reading in RefFreeEWAS Functions ### ----
source("./RefFreeEWAS/R/PairsBootRefFreeEwasModel.R")
source("./RefFreeEWAS/R/RefFreeCellMix.R")
source("./RefFreeEWAS/R/RefFreeEWAS.R")

### Reading in Data ### ----
Bval_Saliva <- read.csv("Bval_Saliva_TopVar.csv")[,-1]
Bval_Saliva <- Bval_Saliva[complete.cases(Bval_Saliva),]
Bval_Saliva <- Bval_Saliva[,2:(length(Bval_Saliva)-1)]

# Get PCs (without standardization)
svSel = svd(Bval_Saliva)

# Initial deconvolution #
cellmixArray  <- RefFreeCellMixArrayWithCustomStart(Bval_Saliva, 
                                                    mu.start = svSel$u,  # Initial methylome matrix 
                                                    Klist=1:15,          # List of K values to try (# constituent cell types)
                                                    iters=25             # Number of iterations per value of K
)

########################
# Bootstrap
# Do the bootstrap for selecting the K parameter (# assumed cell types)
cellmixArrayBoot <- RefFreeCellMixArrayDevianceBoots(
  cellmixArray,            # Array object
  Y=Bval_Saliva,                  # Data on which array was based
  R=1000,                   # 100 bootstraps for example (should really use 500 or so)
  bootstrapIterations=10)  # Num iterations per bootstrap (can be smaller than above)

########################
# Selecting K
# Show winsorized mean deviance per K
wnsrMeanDev <-apply(cellmixArrayBoot[-1,], 2, mean, trim=0.25)
# Write results #
write.csv(wnsrMeanDev, "Saliva_Winsorized_Deviance.csv")

# Write as R data #
saveRDS(cellmixArrayBoot, "Saliva.data.RDS")


