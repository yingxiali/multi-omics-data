# OpenML dataset ids to querry
#install.packages("OpenML")
#install.packages("mlr")
#install.packages("ParamHelpers")
#install.packages("farff")
library(OpenML)
library(mlr)
library(ParamHelpers)
library(farff)



# --> This file can be downloaded from the Github page in
# the subfolder "data":
load("E:/18_omics_datas/1_download/datset_ids.RData")
setwd("E:/18_omics_datas/1_download")
getwd()

nams <- c("LAML", "BLCA", "LGG",  "COAD", "ESCA", 
          "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", 
           "OV", "PAAD", "SARC", "SKCM", "STAD", "UCEC",
          "PCPG", "PRAD", "TGCT",  "THCA", "THYM")

#nam <- "UCEC"
for(nam in nams){
  # download dataset
  dat_part1 <- getOMLDataSet(datset_ids[[nam]][[1]])
  dat_part2 <- getOMLDataSet(datset_ids[[nam]][[2]])
  
  dat <- cbind.data.frame(dat_part1, dat_part2)
  
  
  blocknames <- c("clinical", "cnv", "mirna", "mutation", "rna")
  
  blockinds <- lapply(paste0("_", blocknames), function(x) grep(x, names(dat)))
  # --> blockinds is list of length 5, where the first list element contains the indices
  # of the clinical data, the second list element that of the cnv data
  # and so on (see "blocknames" above).
  
  
  clindata <- dat[,blockinds[[1]]]
  cnvdata <- dat[,blockinds[[2]]]
  mirnadata <- dat[,blockinds[[3]]]
  mutationdata <- dat[,blockinds[[4]]]
  rnadata <- dat[,blockinds[[5]]]
  
  #head(names(dat))
  # --> "bcr_patient_barcode", "time", und "status" do not belong to the covariates
  # and have to be removed.
  
  save(clindata, cnvdata, mirnadata,mutationdata, rnadata,
       file = paste(nam,".RData", sep = ""))
}

#load("E:/18_omics_datas/UCEC.RData")

#########download "BRCA", due to "BRCA" has 3 parts

nam <- "BRCA"


dat_part1 <- getOMLDataSet(datset_ids[[nam]][[1]])
dat_part2 <- getOMLDataSet(datset_ids[[nam]][[2]])

dat <- cbind.data.frame(dat_part1, dat_part2)

dat_part3 <- getOMLDataSet(datset_ids[[nam]][[3]])
dat <- cbind.data.frame(dat, dat_part3)


blocknames <- c("clinical", "cnv", "mirna", "mutation", "rna")

blockinds <- lapply(paste0("_", blocknames), function(x) grep(x, names(dat)))
# --> blockinds is list of length 5, where the first list element contains the indices
# of the clinical data, the second list element that of the cnv data
# and so on (see "blocknames" above).


clindata <- dat[,blockinds[[1]]]
cnvdata <- dat[,blockinds[[2]]]
mirnadata <- dat[,blockinds[[3]]]
mutationdata <- dat[,blockinds[[4]]]
rnadata <- dat[,blockinds[[5]]]

#head(names(dat))
# --> "bcr_patient_barcode", "time", und "status" do not belong to the covariates
# and have to be removed.

save(clindata, cnvdata, mirnadata,mutationdata, rnadata,
     file = paste(nam,".RData", sep = ""))
