########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data/Results")

########################################################

##### load data ####
library("grid")
library("gridExtra")
library("ggplot2")
library("ggpubr")
library("patchwork")
library(stringr)
library(RColorBrewer)
load("./rda_files/resultsumsum.RData")

#### prepare heatmap 1 ####
mnames_bf <- c("bf.rna","bf.mirna","bf.methy","bf.mutation","bf.cnv",
            
            "bf.rna_mirna","bf.rna_cnv","bf.rna_methy","bf.rna_mutation","bf.mirna_methy", 
            "bf.mirna_mutation", "bf.miran_cnv", "bf.methy_mutation", "bf.methy_cnv", "bf.mutation_cnv",
            
            "bf.miran_methy_cnv", "bf.miran_mutation_cnv", "bf.rna_mirna_methy","bf.rna_mirna_mutation","bf.rna_mirna_cnv",
            "bf.rna_methy_mutation","bf.rna_methy_cnv", "bf.rna_mutation_cnv",
            "bf.mirna_methy_mutation", "bf.methy_mutation_cnv",
            
            "bf.rna_mirna_methy_mutation","bf.rna_mirna_methy_cnv", "bf.rna_mirna_mutation_cnv", 
            "bf.rna_methy_mutation_cnv", "bf.mirna_methy_mutation_cnv",
            
            "bf.rna_mirna_methy_mutation_cnv")
mnames_rf <- c("rf.rna","rf.mirna","rf.methy","rf.mutation","rf.cnv",
               
               "rf.rna_mirna","rf.rna_cnv","rf.rna_methy","rf.rna_mutation","rf.mirna_methy", 
               "rf.mirna_mutation", "rf.miran_cnv", "rf.methy_mutation", "rf.methy_cnv", "rf.mutation_cnv",
               
               "rf.miran_methy_cnv", "rf.miran_mutation_cnv", "rf.rna_mirna_methy","rf.rna_mirna_mutation","rf.rna_mirna_cnv",
               "rf.rna_methy_mutation","rf.rna_methy_cnv", "rf.rna_mutation_cnv",
               "rf.mirna_methy_mutation", "rf.methy_mutation_cnv",
               
               "rf.rna_mirna_methy_mutation","rf.rna_mirna_methy_cnv", "rf.rna_mirna_mutation_cnv", 
               "rf.rna_methy_mutation_cnv", "rf.mirna_methy_mutation_cnv",
               
               "rf.rna_mirna_methy_mutation_cnv")
mnames_lasso <- c("lasso.rna","lasso.mirna","lasso.methy","lasso.mutation","lasso.cnv",
               
               "lasso.rna_mirna","lasso.rna_cnv","lasso.rna_methy","lasso.rna_mutation","lasso.mirna_methy", 
               "lasso.mirna_mutation", "lasso.miran_cnv", "lasso.methy_mutation", "lasso.methy_cnv", "lasso.mutation_cnv",
               
               "lasso.miran_methy_cnv", "lasso.miran_mutation_cnv", "lasso.rna_mirna_methy","lasso.rna_mirna_mutation","lasso.rna_mirna_cnv",
               "lasso.rna_methy_mutation","lasso.rna_methy_cnv", "lasso.rna_mutation_cnv",
               "lasso.mirna_methy_mutation", "lasso.methy_mutation_cnv",
               
               "lasso.rna_mirna_methy_mutation","lasso.rna_mirna_methy_cnv", "lasso.rna_mirna_mutation_cnv", 
               "lasso.rna_methy_mutation_cnv", "lasso.mirna_methy_mutation_cnv",
               
               "lasso.rna_mirna_methy_mutation_cnv")
mnames_ipflasso <- c("ipflasso.rna","ipflasso.mirna","ipflasso.methy","ipflasso.mutation","ipflasso.cnv",
               
               "ipflasso.rna_mirna","ipflasso.rna_cnv","ipflasso.rna_methy","ipflasso.rna_mutation","ipflasso.mirna_methy", 
               "ipflasso.mirna_mutation", "ipflasso.miran_cnv", "ipflasso.methy_mutation", "ipflasso.methy_cnv", "ipflasso.mutation_cnv",
               
               "ipflasso.miran_methy_cnv", "ipflasso.miran_mutation_cnv", "ipflasso.rna_mirna_methy","ipflasso.rna_mirna_mutation","ipflasso.rna_mirna_cnv",
               "ipflasso.rna_methy_mutation","ipflasso.rna_methy_cnv", "ipflasso.rna_mutation_cnv",
               "ipflasso.mirna_methy_mutation", "ipflasso.methy_mutation_cnv",
               
               "ipflasso.rna_mirna_methy_mutation","ipflasso.rna_mirna_methy_cnv", "ipflasso.rna_mirna_mutation_cnv", 
               "ipflasso.rna_methy_mutation_cnv", "ipflasso.mirna_methy_mutation_cnv",
               
               "ipflasso.rna_mirna_methy_mutation_cnv")
mnames_prioritylasso <- c("prioritylasso.rna","prioritylasso.mirna","prioritylasso.methy","prioritylasso.mutation","prioritylasso.cnv",
               
               "prioritylasso.rna_mirna","prioritylasso.rna_cnv","prioritylasso.rna_methy","prioritylasso.rna_mutation","prioritylasso.mirna_methy", 
               "prioritylasso.mirna_mutation", "prioritylasso.miran_cnv", "prioritylasso.methy_mutation", "prioritylasso.methy_cnv", "prioritylasso.mutation_cnv",
               
               "prioritylasso.miran_methy_cnv", "prioritylasso.miran_mutation_cnv", "prioritylasso.rna_mirna_methy","prioritylasso.rna_mirna_mutation","prioritylasso.rna_mirna_cnv",
               "prioritylasso.rna_methy_mutation","prioritylasso.rna_methy_cnv", "prioritylasso.rna_mutation_cnv",
               "prioritylasso.mirna_methy_mutation", "prioritylasso.methy_mutation_cnv",
               
               "prioritylasso.rna_mirna_methy_mutation","prioritylasso.rna_mirna_methy_cnv", "prioritylasso.rna_mirna_mutation_cnv", 
               "prioritylasso.rna_methy_mutation_cnv", "prioritylasso.mirna_methy_mutation_cnv",
               
               "prioritylasso.rna_mirna_methy_mutation_cnv")
vals <- c("1","0","0","0","0",
       "1","1","1","1","0",
       "0","0","0","0","0",
       "0","0","1","1","1",
       "1","1","1","0","0",
       "1","1","1","1","0","1",
       
       "0","1","0","0","0",
       "1","0","0","0","1",
       "1","1","0","0","0",
       "1","1","1","1","1",
       "0","0","0","1","0",
       "1","1","1","0","1","1",
       
       "0","0","1","0","0",
       "0","0","1","0","1",
       "0","0","1","1","0",
       "1","0","1","0","0",
       "1","1","0","1","1",
       "1","1","0","1","1","1",
       
       "0","0","0","1","0",
       "0","0","0","1","0",
       "1","0","1","0","1",
       "0","1","0","1","0",
       "1","0","1","1","1",
       "1","0","1","1","1","1",
       
       "0","0","0","0","1",
       "0","1","0","0","0",
       "0","1","0","1","1",
       "1","1","0","0","1",
       "0","1","1","0","1",
       "0","1","1","1","1","1"
       
)
a <- c(rep(mnames_bf,5),rep(mnames_rf,5),rep(mnames_lasso,5),rep(mnames_ipflasso,5),rep(mnames_prioritylasso,5))

b <- c(rep("bf.rna", 31),rep("bf.mirna", 31),rep("bf.methy", 31),rep("bf.mutation", 31),rep("bf.cnv", 31),
       rep("rf.rna", 31),rep("rf.mirna", 31),rep("rf.methy", 31),rep("rf.mutation", 31),rep("rf.cnv", 31),
       rep("lasso.rna", 31),rep("lasso.mirna", 31),rep("lasso.methy", 31),rep("lasso.mutation", 31),rep("lasso.cnv", 31),
       rep("ipflasso.rna", 31),rep("ipflasso.mirna", 31),rep("ipflasso.methy", 31),rep("ipflasso.mutation", 31),rep("ipflasso.cnv", 31),
       rep("prioritylasso.rna", 31),rep("prioritylasso.mirna", 31),rep("prioritylasso.methy", 31),rep("prioritylasso.mutation", 31),rep("prioritylasso.cnv", 31))

c <- c(rep(vals,5))

h1 <- data.frame(as.factor(a),as.factor(b),as.character(c))
colnames(h1) <- c("a","b","c")

###### cindex #####

#block forest
resultscindex_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_bf)
colnames(resultscindex_bf) <- c("comb", "dat","cindex_bf")
resultswide_bf <- reshape(resultscindex_bf , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#random forest
resultscindex_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_rf)
colnames(resultscindex_rf) <- c("comb", "dat","cindex_rf")
resultswide_rf <- reshape(resultscindex_rf , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#lasso
resultscindex_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_lasso)
colnames(resultscindex_lasso) <- c("comb", "dat","cindex_lasso")
resultswide_lasso <- reshape(resultscindex_lasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#ipflasso
resultscindex_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_ipflasso)
colnames(resultscindex_ipflasso) <- c("comb", "dat","cindex_ipflasso")
resultswide_ipflasso <- reshape(resultscindex_ipflasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#priority lasso
resultscindex_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_prioritylasso)
colnames(resultscindex_prioritylasso) <- c("comb", "dat","cindex_prioritylasso")
resultswide_prioritylasso <- reshape(resultscindex_prioritylasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#combination
resultswide1 <- cbind(resultswide_bf,resultswide_rf,resultswide_lasso,resultswide_ipflasso,resultswide_prioritylasso)

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))
colnames(resultranks) <- gsub("cindex.", "", colnames(resultranks))
resultrankstemp <- data.frame(resultranks)
resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
#resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)
means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$sort <- seq(1,length(means$rank),1)

#means$combin
top30 <- merge(means[1:30,], resultranks2, by="combin")

# calculate CI
temp <- apply(resultrankstemp, 2, function(x) quantile(x, c(0.025, 0.975))) 
temp <- rbind(apply(resultrankstemp, 2, mean), temp)
CI_top30 <- temp[,match(means[1:30,]$combin,colnames(temp))]


p1 <- ggplot(data=top30, aes(x=combin, y=rank.y)) + 
  theme_bw() + 
  geom_boxplot() +
  ylim(0,155)+
  scale_x_discrete(limits=means[1:30,]$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Data set specific ranks")


h1 <- data.frame(as.factor(a),as.factor(b),as.character(c))
colnames(h1) <- c("a","b","c")
h1[,c('classifiers','b')] <- str_split_fixed(h1$b, "[.]", 2)
top30h1 <- merge(means[1:30,], h1, by.x ="combin", by.y ="a")

p2<-ggplot(top30h1, aes(x=combin, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"), 
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means[1:30,]$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')


##### prepare heatmap 2 ####
h2 <- data.frame(means$combin)
h2[,c('classifiers','combines')] <- str_split_fixed(h2$means.combin, "[.]", 2)
h2 <- cbind(h2,h2$classifiers)
colnames(h2) <- c("combines",  "classifiers" ,  "combine" ,     "classifier")

h2$w <- as.factor(rep(1, length(h2$combines)))
top30h2 <- merge(means[1:30,], h2, by.x ="combin", by.y ="combines")
top30h2 <- top30h2[order(top30h2$sort),]


p3<-ggplot(top30h2, aes(x=combin, y=w, fill=classifiers)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("rf"="blue","bf"="hotpink","lasso"="red", 
                               "ipflasso"="orange","prioritylasso"="green"),
                    labels = c("rsf","bf","lasso","ipflasso","prioritylasso"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_x_discrete(limits=top30h2$combin)+
  scale_y_discrete(labels=c("predictors")) +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank(),
        axis.ticks.x=element_blank(),
        legend.margin=margin(t=-15),
        legend.key.size = unit(8, "pt"))


## draw figure
top30_pall_cindex <- ggarrange(p1, NULL, p2, NULL, p3, nrow = 5, align="v",heights = c(1, -0.23, 1,-0.70, 1) )
top30_pall_cindex
ggsave(file="./figures/figureS6.png", 
       top30_pall_cindex, width=8, height=5.5)







###### ibrier #####
#### prepare heatmap 1 ####

a <- c(rep(mnames_bf,5),rep(mnames_rf,5),rep(mnames_lasso,5),rep(mnames_ipflasso,5),rep(mnames_prioritylasso,5))

b <- c(rep("bf.rna", 31),rep("bf.mirna", 31),rep("bf.methy", 31),rep("bf.mutation", 31),rep("bf.cnv", 31),
       rep("rf.rna", 31),rep("rf.mirna", 31),rep("rf.methy", 31),rep("rf.mutation", 31),rep("rf.cnv", 31),
       rep("lasso.rna", 31),rep("lasso.mirna", 31),rep("lasso.methy", 31),rep("lasso.mutation", 31),rep("lasso.cnv", 31),
       rep("ipflasso.rna", 31),rep("ipflasso.mirna", 31),rep("ipflasso.methy", 31),rep("ipflasso.mutation", 31),rep("ipflasso.cnv", 31),
       rep("prioritylasso.rna", 31),rep("prioritylasso.mirna", 31),rep("prioritylasso.methy", 31),rep("prioritylasso.mutation", 31),rep("prioritylasso.cnv", 31))

c <- c(rep(vals,5))

h1 <- data.frame(as.factor(a),as.factor(b),as.character(c))
colnames(h1) <- c("a","b","c")

###### ibrier #####

#block forest
resultsibrier_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_bf)
colnames(resultsibrier_bf) <- c("comb", "dat","ibrier_bf")
resultswide_bf <- reshape(resultsibrier_bf , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#random forest
resultsibrier_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_rf)
colnames(resultsibrier_rf) <- c("comb", "dat","ibrier_rf")
resultswide_rf <- reshape(resultsibrier_rf , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#lasso
resultsibrier_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_lasso)
colnames(resultsibrier_lasso) <- c("comb", "dat","ibrier_lasso")
resultswide_lasso <- reshape(resultsibrier_lasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#ipflasso
resultsibrier_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_ipflasso)
colnames(resultsibrier_ipflasso) <- c("comb", "dat","ibrier_ipflasso")
resultswide_ipflasso <- reshape(resultsibrier_ipflasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#priority lasso
resultsibrier_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_prioritylasso)
colnames(resultsibrier_prioritylasso) <- c("comb", "dat","ibrier_prioritylasso")
resultswide_prioritylasso <- reshape(resultsibrier_prioritylasso , idvar = c("dat"), timevar = "comb", direction = "wide")[,-1]
#combination
resultswide1 <- cbind(resultswide_bf,resultswide_rf,resultswide_lasso,resultswide_ipflasso,resultswide_prioritylasso)

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))
colnames(resultranks) <- gsub("ibrier.", "", colnames(resultranks))
resultrankstemp <- data.frame(resultranks)
resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
#resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)
means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$sort <- seq(1,length(means$rank),1)

#means$combin
top30 <- merge(means[1:30,], resultranks2, by="combin")

# calculate CI
temp <- apply(resultrankstemp, 2, function(x) quantile(x, c(0.025, 0.975))) 
temp <- rbind(apply(resultrankstemp, 2, mean), temp)
CI_top30 <- temp[,match(means[1:30,]$combin,colnames(temp))]


p1 <- ggplot(data=top30, aes(x=combin, y=rank.y)) + 
  theme_bw() + 
  geom_boxplot() +
  ylim(0,155)+
  scale_x_discrete(limits=means[1:30,]$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Data set specific ranks")
p1

h1 <- data.frame(as.factor(a),as.factor(b),as.character(c))
colnames(h1) <- c("a","b","c")
h1[,c('classifiers','b')] <- str_split_fixed(h1$b, "[.]", 2)
top30h1 <- merge(means[1:30,], h1, by.x ="combin", by.y ="a")

p2<-ggplot(top30h1, aes(x=combin, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"), 
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means[1:30,]$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')
p2

##### prepare heatmap 2 ####
h2 <- data.frame(means$combin)
h2[,c('classifiers','combines')] <- str_split_fixed(h2$means.combin, "[.]", 2)
h2 <- cbind(h2,h2$classifiers)
colnames(h2) <- c("combines",  "classifiers" ,  "combine" ,     "classifier")

h2$w <- as.factor(rep(1, length(h2$combines)))
top30h2 <- merge(means[1:30,], h2, by.x ="combin", by.y ="combines")
top30h2 <- top30h2[order(top30h2$sort),]

p3<-ggplot(top30h2, aes(x=combin, y=w, fill=classifiers)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("rf"="blue","bf"="hotpink","lasso"="red", "ipflasso"="orange",
                               "prioritylasso"="green"),
                    labels = c("rsf","bf","lasso","ipflasso","prioritylasso"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_x_discrete(limits=top30h2$combin)+
  scale_y_discrete(labels=c("predictors")) +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank(),
        axis.ticks.x=element_blank(),
        legend.margin=margin(t=-15),
        legend.key.size = unit(8, "pt"))
p3

## draw figure
top30_pall_ibrier <- ggarrange(p1, NULL, p2, NULL, p3, nrow = 5, align="v",heights = c(1, -0.23, 1,-0.70, 1) )
top30_pall_ibrier
ggsave(file="./figures/figure2.png", 
       top30_pall_ibrier, width=8, height=5.5)
