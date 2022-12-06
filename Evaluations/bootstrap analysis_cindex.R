########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data/Data")

########################################################

##### load data ####
rm(list = ls())
library(bootstrap)
load("./resultsumsum.RData")
library(dplyr)
library(Rmisc)
library(xlsx)
manes <- c("rna","mirna","methy","mutation","cnv",
           
           "rna_mirna","rna_cnv","rna_methy","rna_mutation","mirna_methy", 
           "mirna_mutation", "miran_cnv", "methy_mutation", "methy_cnv", "mutation_cnv",
           
           "miran_methy_cnv", "miran_mutation_cnv", "rna_mirna_methy","rna_mirna_mutation","rna_mirna_cnv",
           "rna_methy_mutation","rna_methy_cnv", "rna_mutation_cnv",
           "mirna_methy_mutation", "methy_mutation_cnv",
           
           "rna_mirna_methy_mutation","rna_mirna_methy_cnv", "rna_mirna_mutation_cnv", 
           "rna_methy_mutation_cnv", "mirna_methy_mutation_cnv",
           
           "rna_mirna_methy_mutation_cnv")

##### block forest ####
resultscindex_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_bf)
colnames(resultscindex_bf) <- c("comb", "dat","cindex_bf")
resultswide_bf <- reshape(resultscindex_bf , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_bf) <- gsub("cindex_bf.", "", colnames(resultswide_bf))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_bf[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))
  resultrankstemp <- data.frame(resultranks)
  resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                          v.names="rank", 
                          timevar="combin", times=colnames(resultranks),
                          direction="long")
  means <- aggregate(rank ~  combin, resultranks2, mean)
  means <- means[order(means$rank),]
  means$sort <- seq(1,length(means$rank),1)
  means <- means[order(means$combin),]
  Ranks[[i]] <- means[,3]
  
}

timeend <-Sys.time()
runningtime <-timeend-timestart
allranks = data.frame(do.call(rbind, Ranks))

colnames(allranks) <- means[,1]
temp <- apply(allranks, 2, function(x) quantile(x, c(0.025, 0.975))) 
temp <- rbind(apply(allranks, 2, mean), temp)
CI_bf <- temp[,match(manes,colnames(temp))]


##### random forest ####
resultscindex_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_rf)
colnames(resultscindex_rf) <- c("comb", "dat","cindex_rf")
resultswide_rf <- reshape(resultscindex_rf , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_rf) <- gsub("cindex_rf.", "", colnames(resultswide_rf))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_rf[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))
  resultrankstemp <- data.frame(resultranks)
  resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                          v.names="rank", 
                          timevar="combin", times=colnames(resultranks),
                          direction="long")
  means <- aggregate(rank ~  combin, resultranks2, mean)
  means <- means[order(means$rank),]
  means$sort <- seq(1,length(means$rank),1)
  means <- means[order(means$combin),]
  Ranks[[i]] <- means[,3]
  
}

timeend <-Sys.time()
runningtime <-timeend-timestart
allranks = data.frame(do.call(rbind, Ranks))

colnames(allranks) <- means[,1]
temp <- apply(allranks, 2, function(x) quantile(x, c(0.025, 0.975))) 
temp <- rbind(apply(allranks, 2, mean), temp)
CI_rf <- temp[,match(manes,colnames(temp))]


##### lasso ####
resultscindex_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_lasso)
colnames(resultscindex_lasso) <- c("comb", "dat","cindex_lasso")
resultswide_lasso <- reshape(resultscindex_lasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_lasso) <- gsub("cindex_lasso.", "", colnames(resultswide_lasso))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_lasso[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))
  resultrankstemp <- data.frame(resultranks)
  resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                          v.names="rank", 
                          timevar="combin", times=colnames(resultranks),
                          direction="long")
  means <- aggregate(rank ~  combin, resultranks2, mean)
  means <- means[order(means$rank),]
  means$sort <- seq(1,length(means$rank),1)
  means <- means[order(means$combin),]
  Ranks[[i]] <- means[,3]
  
}

timeend <-Sys.time()
runningtime <-timeend-timestart
allranks = data.frame(do.call(rbind, Ranks))
colnames(allranks) <- means[,1]
temp <- apply(allranks, 2, function(x) quantile(x, c(0.025, 0.975))) 
temp <- rbind(apply(allranks, 2, mean), temp)
CI_lasso <- temp[,match(manes,colnames(temp))]


##### ipflasso ####

resultscindex_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_ipflasso)
colnames(resultscindex_ipflasso) <- c("comb", "dat","cindex_ipflasso")
resultswide_ipflasso <- reshape(resultscindex_ipflasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_ipflasso) <- gsub("cindex_ipflasso.", "", colnames(resultswide_ipflasso))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_ipflasso[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))
  resultrankstemp <- data.frame(resultranks)
  resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                          v.names="rank", 
                          timevar="combin", times=colnames(resultranks),
                          direction="long")
  means <- aggregate(rank ~  combin, resultranks2, mean)
  means <- means[order(means$rank),]
  means$sort <- seq(1,length(means$rank),1)
  means <- means[order(means$combin),]
  Ranks[[i]] <- means[,3]
  
}

timeend <-Sys.time()
runningtime <-timeend-timestart
allranks = data.frame(do.call(rbind, Ranks))
colnames(allranks) <- means[,1]
temp <- apply(allranks, 2, function(x) quantile(x, c(0.025, 0.975))) 
temp <- rbind(apply(allranks, 2, mean), temp)
CI_ipflasso <- temp[,match(manes,colnames(temp))]

##### prioritylasso ####
resultscindex_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_prioritylasso)
colnames(resultscindex_prioritylasso) <- c("comb", "dat","cindex_prioritylasso")
resultswide_prioritylasso <- reshape(resultscindex_prioritylasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_prioritylasso) <- gsub("cindex_prioritylasso.", "", colnames(resultswide_prioritylasso))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_prioritylasso[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))
  resultrankstemp <- data.frame(resultranks)
  resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                          v.names="rank", 
                          timevar="combin", times=colnames(resultranks),
                          direction="long")
  means <- aggregate(rank ~  combin, resultranks2, mean)
  means <- means[order(means$rank),]
  means$sort <- seq(1,length(means$rank),1)
  means <- means[order(means$combin),]
  Ranks[[i]] <- means[,3]
  
}

timeend <-Sys.time()
runningtime <-timeend-timestart
allranks = data.frame(do.call(rbind, Ranks))

colnames(allranks) <- means[,1]
temp <- apply(allranks, 2, function(x) quantile(x, c(0.025, 0.975))) 
temp <- rbind(apply(allranks, 2, mean), temp)
CI_priority <- temp[,match(manes,colnames(temp))]


##combination
CI_cindex <- rbind(CI_bf, CI_rf, CI_lasso, CI_ipflasso, CI_priority)
rownames(CI_cindex) <- c("bf_mean","bf_lower", "bf_upper",
                         "rf_mean","rf_lower", "rf_upper",
                         "lasso_mean","lasso_lower", "lasso_upper",
                         "ipflasso_mean","ipflasso_lower", "ipflasso_upper",
                         "prioritylasso_mean","prioritylasso_lower", "prioritylasso_upper")

save(CI_cindex, file = "./CI_cindex.RData")
