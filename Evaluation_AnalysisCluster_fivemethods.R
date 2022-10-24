rm(list = ls())

###### get the data of five method. ####

setwd("C:/Users/yingxiali/Desktop/paper3/LRZ_Jul_Results")

##### prepared the results for 1-4-5-omics data 
# part1. Read results 
load("./scenariogrid1.Rda")
scenariogrid1 <- scenariogrid

allinds <- 1:nrow(scenariogrid1) # a vector

allresults <- list.files("./Results1") # a vector
indres <- gsub("res", "", allresults) # delete "res"
indres1 <- sort(as.numeric(gsub(".Rda", "", indres))) # delete ".Rda"

(anymissing <- length(setdiff(allinds, indres1)) > 0) # setdiff()

scenariogrid1 <- scenariogrid1[indres1,] # finished iteration

# part2. got the missing iterations from the part 1
#note: we need to rerun the scenmiss on the LRZ

(scenmiss1 <- scenariogrid[setdiff(allinds, indres1),]) # get the remaining cases

# part3. read iteration results

Results <- list()

count <- 1

counttemp <- 1
for(i in 1:nrow(scenariogrid1)) {
  load(paste("./Results1/res", indres1[counttemp], ".Rda", sep=""))
  Results[[count]] <- res
  count <- count + 1
  counttemp <- counttemp + 1
}


ibrier_rf <- lapply(rapply(sapply(Results, function(x) x$ibrier_rf), enquote, how="unlist"), eval)
cindex_rf <- lapply(rapply(sapply(Results, function(x) x$cindex_rf), enquote, how="unlist"), eval)
ibrier_bf <- lapply(rapply(sapply(Results, function(x) x$ibrier_bf), enquote, how="unlist"), eval)
cindex_bf <- lapply(rapply(sapply(Results, function(x) x$cindex_bf), enquote, how="unlist"), eval)
ibrier_prioritylasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_prioritylasso), enquote, how="unlist"), eval)
cindex_prioritylasso <- lapply(rapply(sapply(Results, function(x) x$cindex_prioritylasso), enquote, how="unlist"), eval)
ibrier_ipflasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_ipflasso), enquote, how="unlist"), eval)
cindex_ipflasso <- lapply(rapply(sapply(Results, function(x) x$cindex_ipflasso), enquote, how="unlist"), eval)
ibrier_lasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_lasso), enquote, how="unlist"), eval)
cindex_lasso <- lapply(rapply(sapply(Results, function(x) x$cindex_lasso), enquote, how="unlist"), eval)



ibrier_bf <-unlist(ibrier_bf)
cindex_bf <-unlist(cindex_bf)
ibrier_rf <-unlist(ibrier_rf)
cindex_rf <-unlist(cindex_rf)
ibrier_prioritylasso <-unlist(ibrier_prioritylasso)
cindex_prioritylasso <-unlist(cindex_prioritylasso)
ibrier_ipflasso <-unlist(ibrier_ipflasso)
cindex_ipflasso <-unlist(cindex_ipflasso)
ibrier_lasso <-unlist(ibrier_lasso)
cindex_lasso <-unlist(cindex_lasso)


results <- cbind(scenariogrid1,ibrier_rf,cindex_rf,ibrier_bf,cindex_bf,
                 ibrier_lasso, cindex_lasso, ibrier_prioritylasso, cindex_prioritylasso,
                 ibrier_ipflasso, cindex_ipflasso)#,ibrier_prioritylasso, cindex_prioritylasso
results$seed <- NULL
table(is.na(results))
table(is.na(results$ibrier_prioritylasso))
table(is.na(results$cindex_prioritylasso))
table(is.na(results$ibrier_bf))
table(is.na(results$cindex_bf))
results[is.na(results$ibrier_prioritylasso),which(colnames(results) == 'ibrier_prioritylasso')] <- 0.25
results[is.na(results$cindex_prioritylasso),which(colnames(results) == 'cindex_prioritylasso')] <- 0.5
results[is.na(results$ibrier_ipflasso),which(colnames(results) == 'ibrier_ipflasso')] <- 0.25
results[is.na(results$cindex_ipflasso),which(colnames(results) == 'cindex_ipflasso')] <- 0.5
results[is.na(results$ibrier_lasso),which(colnames(results) == 'ibrier_lasso')] <- 0.25
results[is.na(results$cindex_lasso),which(colnames(results) == 'cindex_lasso')] <- 0.5
results[is.na(results$ibrier_bf),which(colnames(results) == 'ibrier_bf')] <- 0.25
results[is.na(results$cindex_bf),which(colnames(results) == 'cindex_bf')] <- 0.5
results[is.na(results$ibrier_rf),which(colnames(results) == 'ibrier_rf')] <- 0.25
results[is.na(results$cindex_rf),which(colnames(results) == 'cindex_rf')] <- 0.5
table(is.na(results))
#save(results, file="./results.Rda")


# Calculate mean values per data set:

library("plyr")
resultsum145 <- ddply(results, .variables=c("comb", "dat", "cvind"), .fun=summarise,
                   ibrier_bf=mean(ibrier_bf, na.rm=TRUE),cindex_bf=mean(cindex_bf, na.rm=TRUE),
                   ibrier_rf=mean(ibrier_rf, na.rm=TRUE),cindex_rf=mean(cindex_rf, na.rm=TRUE),
                   ibrier_lasso=mean(ibrier_lasso, na.rm=TRUE),cindex_lasso=mean(cindex_lasso, na.rm=TRUE),
                   ibrier_ipflasso=mean(ibrier_ipflasso, na.rm=TRUE),cindex_ipflasso=mean(cindex_ipflasso, na.rm=TRUE),
                   ibrier_prioritylasso=mean(ibrier_prioritylasso, na.rm=TRUE),cindex_prioritylasso=mean(cindex_prioritylasso, na.rm=TRUE))

resultsum145$comb <- factor(resultsum145$comb, levels=c("rna","mirna","methy","mutation","cnv",
                                                  "rna_mirna_methy_mutation","rna_mirna_methy_cnv", "rna_mirna_mutation_cnv", "rna_methy_mutation_cnv", "mirna_methy_mutation_cnv",
                                                  "rna_mirna_methy_mutation_cnv"))


# Calculate the mean values across data sets:

resultsumsum145 <- ddply(resultsum145, .variables=c("comb", "dat"), .fun=summarise, 
                      ibrier_bf=mean(ibrier_bf, na.rm=TRUE),cindex_bf=mean(cindex_bf, na.rm=TRUE),
                      ibrier_rf=mean(ibrier_rf, na.rm=TRUE),cindex_rf=mean(cindex_rf, na.rm=TRUE),
                      ibrier_lasso=mean(ibrier_lasso, na.rm=TRUE),cindex_lasso=mean(cindex_lasso, na.rm=TRUE),
                      ibrier_ipflasso=mean(ibrier_ipflasso, na.rm=TRUE),cindex_ipflasso=mean(cindex_ipflasso, na.rm=TRUE),
                      ibrier_prioritylasso=mean(ibrier_prioritylasso, na.rm=TRUE),cindex_prioritylasso=mean(cindex_prioritylasso, na.rm=TRUE))

resultswide145 <- reshape(resultsumsum145, idvar = c("dat"), timevar = "comb", direction = "wide")

##### prepared the results for 2-omics data 

# part1. Read results 
#rm(list = ls())
load("./scenariogrid2.Rda")
scenariogrid2 <- scenariogrid

allinds <- 1:nrow(scenariogrid2) # a vector

allresults <- list.files("./Results2") # a vector
indres <- gsub("res", "", allresults) # delete "res"
indres1 <- sort(as.numeric(gsub(".Rda", "", indres))) # delete ".Rda"

(anymissing <- length(setdiff(allinds, indres1)) > 0) # setdiff()

scenariogrid2 <- scenariogrid2[indres1,] # finished iteration

# part2. got the missing iterations from the part 1
#note: we need to rerun the scenmiss on the LRZ

(scenmiss2 <- scenariogrid[setdiff(allinds, indres1),]) # get the remaining cases
#save(scenmiss, file="./scenmiss3.Rda")

# part3. read iteration results

Results <- list()

count <- 1

counttemp <- 1
for(i in 1:nrow(scenariogrid2)) {
  load(paste("./Results2/res", indres1[counttemp], ".Rda", sep=""))
  Results[[count]] <- res
  count <- count + 1
  counttemp <- counttemp + 1
}

#scenariogrid$method[scenariogrid$method=="randomsurvivalforest"] <- "RSF"

ibrier_rf <- lapply(rapply(sapply(Results, function(x) x$ibrier_rf), enquote, how="unlist"), eval)
cindex_rf <- lapply(rapply(sapply(Results, function(x) x$cindex_rf), enquote, how="unlist"), eval)
ibrier_bf <- lapply(rapply(sapply(Results, function(x) x$ibrier_bf), enquote, how="unlist"), eval)
cindex_bf <- lapply(rapply(sapply(Results, function(x) x$cindex_bf), enquote, how="unlist"), eval)
ibrier_prioritylasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_prioritylasso), enquote, how="unlist"), eval)
cindex_prioritylasso <- lapply(rapply(sapply(Results, function(x) x$cindex_prioritylasso), enquote, how="unlist"), eval)
ibrier_ipflasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_ipflasso), enquote, how="unlist"), eval)
cindex_ipflasso <- lapply(rapply(sapply(Results, function(x) x$cindex_ipflasso), enquote, how="unlist"), eval)
ibrier_lasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_lasso), enquote, how="unlist"), eval)
cindex_lasso <- lapply(rapply(sapply(Results, function(x) x$cindex_lasso), enquote, how="unlist"), eval)



ibrier_bf <-unlist(ibrier_bf)
cindex_bf <-unlist(cindex_bf)
ibrier_rf <-unlist(ibrier_rf)
cindex_rf <-unlist(cindex_rf)
ibrier_prioritylasso <-unlist(ibrier_prioritylasso)
cindex_prioritylasso <-unlist(cindex_prioritylasso)
ibrier_ipflasso <-unlist(ibrier_ipflasso)
cindex_ipflasso <-unlist(cindex_ipflasso)
ibrier_lasso <-unlist(ibrier_lasso)
cindex_lasso <-unlist(cindex_lasso)


results <- cbind(scenariogrid2,ibrier_rf,cindex_rf,ibrier_bf,cindex_bf,
                 ibrier_lasso, cindex_lasso, ibrier_prioritylasso, cindex_prioritylasso,
                 ibrier_ipflasso, cindex_ipflasso)#,ibrier_prioritylasso, cindex_prioritylasso
results$seed <- NULL
table(is.na(results))
table(is.na(results$ibrier_prioritylasso))
table(is.na(results$cindex_prioritylasso))
table(is.na(results$ibrier_bf))
table(is.na(results$cindex_bf))
results[is.na(results$ibrier_prioritylasso),which(colnames(results) == 'ibrier_prioritylasso')] <- 0.25
results[is.na(results$cindex_prioritylasso),which(colnames(results) == 'cindex_prioritylasso')] <- 0.5
results[is.na(results$ibrier_ipflasso),which(colnames(results) == 'ibrier_ipflasso')] <- 0.25
results[is.na(results$cindex_ipflasso),which(colnames(results) == 'cindex_ipflasso')] <- 0.5
results[is.na(results$ibrier_lasso),which(colnames(results) == 'ibrier_lasso')] <- 0.25
results[is.na(results$cindex_lasso),which(colnames(results) == 'cindex_lasso')] <- 0.5
results[is.na(results$ibrier_bf),which(colnames(results) == 'ibrier_bf')] <- 0.25
results[is.na(results$cindex_bf),which(colnames(results) == 'cindex_bf')] <- 0.5
results[is.na(results$ibrier_rf),which(colnames(results) == 'ibrier_rf')] <- 0.25
results[is.na(results$cindex_rf),which(colnames(results) == 'cindex_rf')] <- 0.5
table(is.na(results))

#save(results, file="./results.Rda")


# Calculate mean values per data set:

#library("plyr")
resultsum2 <- ddply(results, .variables=c("comb", "dat", "cvind"), .fun=summarise,
                   ibrier_bf=mean(ibrier_bf, na.rm=TRUE),cindex_bf=mean(cindex_bf, na.rm=TRUE),
                   ibrier_rf=mean(ibrier_rf, na.rm=TRUE),cindex_rf=mean(cindex_rf, na.rm=TRUE),
                   ibrier_lasso=mean(ibrier_lasso, na.rm=TRUE),cindex_lasso=mean(cindex_lasso, na.rm=TRUE),
                   ibrier_ipflasso=mean(ibrier_ipflasso, na.rm=TRUE),cindex_ipflasso=mean(cindex_ipflasso, na.rm=TRUE),
                   ibrier_prioritylasso=mean(ibrier_prioritylasso, na.rm=TRUE),cindex_prioritylasso=mean(cindex_prioritylasso, na.rm=TRUE))

resultsum2$comb <- factor(resultsum2$comb, levels=c("rna_mirna","rna_cnv","rna_methy","rna_mutation","mirna_methy", 
                                                  "mirna_mutation", "miran_cnv", "methy_mutation", "methy_cnv", "mutation_cnv"))


# Calculate the mean values across data sets:

resultsumsum2 <- ddply(resultsum2, .variables=c("comb", "dat"), .fun=summarise, 
                      ibrier_bf=mean(ibrier_bf, na.rm=TRUE),cindex_bf=mean(cindex_bf, na.rm=TRUE),
                      ibrier_rf=mean(ibrier_rf, na.rm=TRUE),cindex_rf=mean(cindex_rf, na.rm=TRUE),
                      ibrier_lasso=mean(ibrier_lasso, na.rm=TRUE),cindex_lasso=mean(cindex_lasso, na.rm=TRUE),
                      ibrier_ipflasso=mean(ibrier_ipflasso, na.rm=TRUE),cindex_ipflasso=mean(cindex_ipflasso, na.rm=TRUE),
                      ibrier_prioritylasso=mean(ibrier_prioritylasso, na.rm=TRUE),cindex_prioritylasso=mean(cindex_prioritylasso, na.rm=TRUE))

resultswide2 <- reshape(resultsumsum2, idvar = c("dat"), timevar = "comb", direction = "wide")



##### prepared the results for 3-omics data 
# part1. Read results 
#rm(list = ls())
load("./scenariogrid3.Rda")
scenariogrid3 <- scenariogrid

allinds <- 1:nrow(scenariogrid3) # a vector

allresults <- list.files("./Results3") # a vector
indres <- gsub("res", "", allresults) # delete "res"
indres1 <- sort(as.numeric(gsub(".Rda", "", indres))) # delete ".Rda"

(anymissing <- length(setdiff(allinds, indres1)) > 0) # setdiff()

scenariogrid3 <- scenariogrid3[indres1,] # finished iteration

# part2. got the missing iterations from the part 1
#note: we need to rerun the scenmiss on the LRZ

(scenmiss3 <- scenariogrid[setdiff(allinds, indres1),]) # get the remaining cases
#save(scenmiss3, file="./scenmiss3.Rda")

# part3. read iteration results

Results <- list()

count <- 1

counttemp <- 1
for(i in 1:nrow(scenariogrid3)) {
  load(paste("./Results3/res", indres1[counttemp], ".Rda", sep=""))
  Results[[count]] <- res
  count <- count + 1
  counttemp <- counttemp + 1
}

ibrier_rf <- lapply(rapply(sapply(Results, function(x) x$ibrier_rf), enquote, how="unlist"), eval)
cindex_rf <- lapply(rapply(sapply(Results, function(x) x$cindex_rf), enquote, how="unlist"), eval)
ibrier_bf <- lapply(rapply(sapply(Results, function(x) x$ibrier_bf), enquote, how="unlist"), eval)
cindex_bf <- lapply(rapply(sapply(Results, function(x) x$cindex_bf), enquote, how="unlist"), eval)
ibrier_prioritylasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_prioritylasso), enquote, how="unlist"), eval)
cindex_prioritylasso <- lapply(rapply(sapply(Results, function(x) x$cindex_prioritylasso), enquote, how="unlist"), eval)
ibrier_ipflasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_ipflasso), enquote, how="unlist"), eval)
cindex_ipflasso <- lapply(rapply(sapply(Results, function(x) x$cindex_ipflasso), enquote, how="unlist"), eval)
ibrier_lasso <- lapply(rapply(sapply(Results, function(x) x$ibrier_lasso), enquote, how="unlist"), eval)
cindex_lasso <- lapply(rapply(sapply(Results, function(x) x$cindex_lasso), enquote, how="unlist"), eval)

ibrier_bf <-unlist(ibrier_bf)
cindex_bf <-unlist(cindex_bf)
ibrier_rf <-unlist(ibrier_rf)
cindex_rf <-unlist(cindex_rf)
ibrier_prioritylasso <-unlist(ibrier_prioritylasso)
cindex_prioritylasso <-unlist(cindex_prioritylasso)
ibrier_ipflasso <-unlist(ibrier_ipflasso)
cindex_ipflasso <-unlist(cindex_ipflasso)
ibrier_lasso <-unlist(ibrier_lasso)
cindex_lasso <-unlist(cindex_lasso)


results <- cbind(scenariogrid3,ibrier_rf,cindex_rf,ibrier_bf,cindex_bf,
                 ibrier_lasso, cindex_lasso, ibrier_prioritylasso, cindex_prioritylasso,
                 ibrier_ipflasso, cindex_ipflasso)#,ibrier_prioritylasso, cindex_prioritylasso
results$seed <- NULL
table(is.na(results))
table(is.na(results$ibrier_prioritylasso))
table(is.na(results$cindex_prioritylasso))
table(is.na(results$ibrier_bf))
table(is.na(results$cindex_bf))
results[is.na(results$ibrier_prioritylasso),which(colnames(results) == 'ibrier_prioritylasso')] <- 0.25
results[is.na(results$cindex_prioritylasso),which(colnames(results) == 'cindex_prioritylasso')] <- 0.5
results[is.na(results$ibrier_ipflasso),which(colnames(results) == 'ibrier_ipflasso')] <- 0.25
results[is.na(results$cindex_ipflasso),which(colnames(results) == 'cindex_ipflasso')] <- 0.5
results[is.na(results$ibrier_lasso),which(colnames(results) == 'ibrier_lasso')] <- 0.25
results[is.na(results$cindex_lasso),which(colnames(results) == 'cindex_lasso')] <- 0.5
results[is.na(results$ibrier_bf),which(colnames(results) == 'ibrier_bf')] <- 0.25
results[is.na(results$cindex_bf),which(colnames(results) == 'cindex_bf')] <- 0.5
results[is.na(results$ibrier_rf),which(colnames(results) == 'ibrier_rf')] <- 0.25
results[is.na(results$cindex_rf),which(colnames(results) == 'cindex_rf')] <- 0.5
table(is.na(results))



# Calculate mean values per data set:

#library("plyr")
resultsum3 <- ddply(results, .variables=c("comb", "dat", "cvind"), .fun=summarise,
                   ibrier_bf=mean(ibrier_bf, na.rm=TRUE),cindex_bf=mean(cindex_bf, na.rm=TRUE),
                   ibrier_rf=mean(ibrier_rf, na.rm=TRUE),cindex_rf=mean(cindex_rf, na.rm=TRUE),
                   ibrier_lasso=mean(ibrier_lasso, na.rm=TRUE),cindex_lasso=mean(cindex_lasso, na.rm=TRUE),
                   ibrier_ipflasso=mean(ibrier_ipflasso, na.rm=TRUE),cindex_ipflasso=mean(cindex_ipflasso, na.rm=TRUE),
                   ibrier_prioritylasso=mean(ibrier_prioritylasso, na.rm=TRUE),cindex_prioritylasso=mean(cindex_prioritylasso, na.rm=TRUE))

resultsum3$comb <- factor(resultsum3$comb, levels=c("miran_methy_cnv", "miran_mutation_cnv", "rna_mirna_methy","rna_mirna_mutation","rna_mirna_cnv",
                                                  "rna_methy_mutation","rna_methy_cnv", "rna_mutation_cnv",
                                                  "mirna_methy_mutation", "methy_mutation_cnv"))


# Calculate the mean values across data sets:

resultsumsum3 <- ddply(resultsum3, .variables=c("comb", "dat"), .fun=summarise, 
                      ibrier_bf=mean(ibrier_bf, na.rm=TRUE),cindex_bf=mean(cindex_bf, na.rm=TRUE),
                      ibrier_rf=mean(ibrier_rf, na.rm=TRUE),cindex_rf=mean(cindex_rf, na.rm=TRUE),
                      ibrier_lasso=mean(ibrier_lasso, na.rm=TRUE),cindex_lasso=mean(cindex_lasso, na.rm=TRUE),
                      ibrier_ipflasso=mean(ibrier_ipflasso, na.rm=TRUE),cindex_ipflasso=mean(cindex_ipflasso, na.rm=TRUE),
                      ibrier_prioritylasso=mean(ibrier_prioritylasso, na.rm=TRUE),cindex_prioritylasso=mean(cindex_prioritylasso, na.rm=TRUE))

resultswide3 <- reshape(resultsumsum3, idvar = c("dat"), timevar = "comb", direction = "wide")



resultsum <- rbind(resultsum145,resultsum2,resultsum3)
save(resultsum, file = "./resultsum.RData")

resultsumsum <- rbind(resultsumsum145,resultsumsum2,resultsumsum3)
save(resultsumsum, file = "./resultsumsum.RData")



