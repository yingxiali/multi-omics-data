########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data/Results/rda_files")

########################################################

##### load data ####
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
resultsibrier_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_bf)
colnames(resultsibrier_bf) <- c("comb", "dat","ibrier_bf")
resultswide_bf <- reshape(resultsibrier_bf , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_bf) <- gsub("ibrier_bf.", "", colnames(resultswide_bf))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_bf[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))
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
#save(CI_bf, file = "C:/Users/yingxiali/Desktop/paper3/3_rcode/CI_ibrier_bf.RData")
#write.xlsx2(CI_bf, file = "C:/Users/yingxiali/Desktop/paper3/3_rcode/CI_ibrier_bf.xlsx",
#col.names = TRUE, row.names = TRUE, append = FALSE)



##### random forest ####
resultsibrier_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_rf)
colnames(resultsibrier_rf) <- c("comb", "dat","ibrier_rf")
resultswide_rf <- reshape(resultsibrier_rf , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_rf) <- gsub("ibrier_rf.", "", colnames(resultswide_rf))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_rf[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))
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
#save(CI_rf, file = "C:/Users/yingxiali/Desktop/paper3/3_rcode/CI_ibrier_rf.RData")
#write.xlsx2(CI_rf, file = "C:/Users/yingxiali/Desktop/paper3/3_rcode/CI_ibrier_rf.xlsx",
#col.names = TRUE, row.names = TRUE, append = FALSE)

##### lasso ####
resultsibrier_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_lasso)
colnames(resultsibrier_lasso) <- c("comb", "dat","ibrier_lasso")
resultswide_lasso <- reshape(resultsibrier_lasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_lasso) <- gsub("ibrier_lasso.", "", colnames(resultswide_lasso))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_lasso[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))
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
#save(CI_lasso, file = "C:/Users/yingxiali/Desktop/paper3/3_rcode/CI_ibrier_lasso.RData")
#write.xlsx2(CI_lasso, file = "C:/Users/yingxiali/Desktop/paper3/3_rcode/CI_ibrier_lasso.xlsx",
#col.names = TRUE, row.names = TRUE, append = FALSE)

##### ipflasso ####

resultsibrier_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_ipflasso)
colnames(resultsibrier_ipflasso) <- c("comb", "dat","ibrier_ipflasso")
resultswide_ipflasso <- reshape(resultsibrier_ipflasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_ipflasso) <- gsub("ibrier_ipflasso.", "", colnames(resultswide_ipflasso))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_ipflasso[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))
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
#save(CI_ipflasso, file = "C:/Users/yingxiali/Desktop/paper3/3_rcode/CI_ibrier_ipflasso.RData")
#write.xlsx2(CI_ipflasso, file = "C:/Users/yingxiali/Desktop/paper3/3_rcode/CI_ibrier_ipflasso.xlsx",
#col.names = TRUE, row.names = TRUE, append = FALSE)

##### prioritylasso ####
resultsibrier_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_prioritylasso)
colnames(resultsibrier_prioritylasso) <- c("comb", "dat","ibrier_prioritylasso")
resultswide_prioritylasso <- reshape(resultsibrier_prioritylasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
colnames(resultswide_prioritylasso) <- gsub("ibrier_prioritylasso.", "", colnames(resultswide_prioritylasso))

set.seed(1234)
timestart <-Sys.time()
Ranks <- list()# the results will be stored
i=1
for (i in 1:5000) {
  #randomly select the indices
  n <- sample(1:14, size = 14, replace = TRUE)
  resultswide1 <- resultswide_prioritylasso[n,]
  # Calculate ranks of the methods:
  resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))
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
CI_ibrier <- rbind(CI_bf, CI_rf, CI_lasso, CI_ipflasso, CI_priority)
rownames(CI_ibrier) <- c("bf_mean","bf_lower", "bf_upper",
                         "rf_mean","rf_lower", "rf_upper",
                         "lasso_mean","lasso_lower", "lasso_upper",
                         "ipflasso_mean","ipflasso_lower", "ipflasso_upper",
                         "prioritylasso_mean","prioritylasso_lower", "prioritylasso_upper")


write.xlsx2(CI_ibrier, file = "./CI_ibrier.xlsx",
            col.names = TRUE, row.names = TRUE, append = FALSE)

