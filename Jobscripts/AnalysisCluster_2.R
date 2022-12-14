########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data")

########################################################

#### 1
# Set working directory:

setwd("~/")

#### 2
# Make table of settings:

dat <- c("BLCA.RData","BRCA.RData", "COAD.RData", "ESCA.RData", "HNSC.RData", 
         "LGG.RData", "LIHC.RData", "LUAD.RData", "LUSC.RData", "PAAD.RData", 
          "SARC.RData", "SKCM.RData", "STAD.RData", "UCEC.RData")

comb <- c( "rna_mirna","rna_cnv","rna_methy","rna_mutation","mirna_methy", 
           "mirna_mutation", "miran_cnv", "methy_mutation", "methy_cnv", "mutation_cnv")

cvind <- 1:5

cvfoldind <- 1:5

scenariogrid <- expand.grid( cvind=cvind,  cvfoldind=cvfoldind,
                             dat=dat, comb=comb, stringsAsFactors = FALSE)
scenariogrid <- scenariogrid[,ncol(scenariogrid):1]


set.seed(1234)
seeds <- sample(1000:10000000, size=length(dat)*length(cvind)) #repeat 3 times

scenariogrid$seed <- rep(seeds, times=length(comb)*length(cvfoldind))


#### 3
# Randomly permute rows of the table containing the settings.
# This is performed to ensure a comparable computational burden for
# the jobs to be performed in parallel:

set.seed(1234)
reorderind <- sample(1:nrow(scenariogrid))
scenariogrid <- scenariogrid[reorderind,]

#### 4
# Save scenariogrid, needed in evaluation of the results:

save(scenariogrid, file="./Results/rda_files/scenariogrid2.Rda") 


#### 5
# Source the functions that are used in performing the calculations 
# on the cluster:

source("./Functions/Functions_AnalysisCluster_2.R")

#### 6
# Start the cluster:

# NOTE: This syntax requires the use of the RMPISNOW script, see the README file
# contained in the root folder "Additional_file_2_HornungWright".

library(snow)

cl <- makeCluster()


#### 7
# Export the objects in the workspace to the
# parallel jobs:

clusterExport(cl, list=ls())

#### 8
# Perform the calculations:

Results <- parLapply(cl, 1:nrow(scenariogrid), function(z)
  try({evaluatesetting(z)}))


#save(Results, file="./VarselCompStudy/Additional_file/Results/Results2/Results_AnalysisCluster.Rda")

#### 9
# Stop the cluster:

stopCluster(cl)
