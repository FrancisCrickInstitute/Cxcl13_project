library(tidyverse)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(knitr)

## Set documentation parameters ##

VersionPdfExt <- paste0(".V", gsub("-", "", Sys.Date()), ".pdf")

if (dir.exists("/Volumes/babs/working/boeings/")){
  hpc.mount <- "/Volumes/babs/working/boeings/"
} else if (dir.exists("Y:/working/boeings/")){
  hpc.mount <- "Y:/working/boeings/"
} else if (dir.exists("/camp/stp/babs/working/boeings/")){
  hpc.mount <- "/camp/stp/babs/working/boeings/"
} else {
  hpc.mount <- ""
}

## Load custom packages specific for this analysis ##
source("assets/scTools.r")
source("assets/SBwebtools.pckg.r")

FN <- paste0(hpc.mount, "Projects/reference_data/documentation/BC.parameters.txt")
dbTable <- read.delim(
  FN, 
  sep = "\t",
  stringsAsFactors = F
)

db.pwd <- as.vector(dbTable[1,1])

## Set filename for temp pdf files ##

if (length(.libPaths()) > 2){
  .libPaths(.libPaths()[2:3])
}


###############################################################################
## Load Obio object from step A                                              ##
ObioFN <- paste0(
  "../", 
  list.files("..")[grep(".bioLOGIC.Robj", list.files(".."))]
)

load(ObioFN)


###############################################################################
## Filter cells 

filterGenes <- c(
  "Lyve1", 
  "Hba-a1", 
  "Hba-a2", 
  "Krt18", 
  "Trac", 
  "Cd3d", 
  "Cldn5", 
  "Ly6c1", 
  "Egfl7", 
  "Ptprc", 
  "S100b", 
  "Cd79a", 
  "Cd79b"
)


library(Seurat)

FN <- "/camp/stp/babs/working/boeings/Projects/tolarp/ana.martinez.riano/423B_scRNAseq_CITE_Cxcl13_SC20209/workdir/"

load("/camp/stp/babs/working/boeings/Projects/tolarp/ana.martinez.riano/423B_scRNAseq_CITE_Cxcl13_SC20209/workdir/SC20209B.Seurat.Robj")

## find filterGenes in this experiment
allGenes <- rownames(x = OsC@assays[["RNA"]])


filterGenes <- filterGenes[filterGenes %in% allGenes]
notPresent <- filterGenes[!(filterGenes %in% allGenes)]

dfSel <- OsC[["RNA"]]@counts
dfSel <- data.frame(dfSel[filterGenes, ])

h <- colSums(dfSel)
h1 <- h[h == 0]
h2 <- h[h > 0]

inCells <- names(h1)
outCells <- names(h2)

dfIn <- data.frame(cellID = inCells, meta_GeneFilter = "Retained")
dfOut <- data.frame(cellID = outCells, meta_GeneFilter = "Removed")
dfF <- rbind(dfIn, dfOut)
row.names(dfF) <- dfF$cellID
dfF$cellID <- NULL



OsC <- addDf2seuratMetaData(
  obj = OsC,
  dfAdd = dfF
)

save(OsC,
     file = paste0(
       Obio@parameterList$localWorkDir,
       Obio@parameterList$project_id,
       ".Seurat.Robj"
     )
)


## Subset data and save ##
OsC_filt <- subset(x = OsC, subset = meta_GeneFilter == "Retained")


# "CXCL13ICPE24h" "Biotin24h"     "Biotin7d"
OsC_CXCL13ICPE24h <- subset(x = OsC_filt, subset = sampleID == "CXCL13ICPE24h")
dfCXCL13ICPE24h <- OsC_CXCL13ICPE24h[["RNA"]]@counts
FN1 <- paste0(
  Obio@parameterList$localWorkDir,
  Obio@parameterList$project_id,
  ".dfCXCL13ICPE24h.matrix.meta_GeneFilter.retained",
  ".txt"
)

write.table(
  dfCXCL13ICPE24h,
  FN1,
  sep = "\t"
)

OsC_Biotin24h <- subset(x = OsC_filt, subset = sampleID == "Biotin24h")
dfBiotin24h <-  OsC_Biotin24h[["RNA"]]@counts
FN2 <- paste0(
  Obio@parameterList$localWorkDir,
  Obio@parameterList$project_id,
  ".dfBiotin24h.matrix.meta_GeneFilter.retained",
  ".txt"
)

write.table(
  dfBiotin24h,
  FN2,
  sep = "\t"
)

OsC_Biotin7d <- subset(x = OsC_filt, subset = sampleID == "Biotin7d")
dfBiotin7d <-  OsC_Biotin7d[["RNA"]]@counts
FN3 <- paste0(
  Obio@parameterList$localWorkDir,
  Obio@parameterList$project_id,
  ".dfBiotin7d.matrix.meta_GeneFilter.retained",
  ".txt"
)

write.table(
  dfBiotin7d,
  FN3,
  sep = "\t"
)

dfMfilt <-  OsC_filt[["RNA"]]@counts

save(OsC_filt,
     file = paste0(
       Obio@parameterList$localWorkDir,
       Obio@parameterList$project_id,
       ".matrix.meta_GeneFilter.retained",
       ".txt"
     )
)

write.table(
    dfMfilt,
    file,
    sep = "\t"
)

# /camp/stp/babs/working/boeings/Projects/tolarp/ana.martinez.riano/423B_scRNAseq_CITE_Cxcl13_SC20209/workdir/SC20209B.matrix.meta_GeneFilter.retained.txt
/camp/stp/babs/working/boeings/Projects/tolarp/ana.martinez.riano/423B_scRNAseq_CITE_Cxcl13_SC20209/workdir/SC20209B.dfCXCL13ICPE24h.matrix.meta_GeneFilter.retained.txt
/camp/stp/babs/working/boeings/Projects/tolarp/ana.martinez.riano/423B_scRNAseq_CITE_Cxcl13_SC20209/workdir/SC20209B.dfBiotin24h.matrix.meta_GeneFilter.retained.txt
/camp/stp/babs/working/boeings/Projects/tolarp/ana.martinez.riano/423B_scRNAseq_CITE_Cxcl13_SC20209/workdir/SC20209B.dfBiotin7d.matrix.meta_GeneFilter.retained.txt

## Done