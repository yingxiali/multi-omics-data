##### load data 
rm(list = ls())
library("grid")
library("gridExtra")
library("ggplot2")
library("ggpubr")
library("patchwork")
setwd("C:/Users/yingxiali/Desktop/paper3/LRZ_Jul_Results")
load("./resultsumsum.RData")


mnames <- c("rna","mirna","methy","mutation","cnv",
            
            "rna_mirna","rna_cnv","rna_methy","rna_mutation","mirna_methy", 
            "mirna_mutation", "miran_cnv", "methy_mutation", "methy_cnv", "mutation_cnv",
            
            "miran_methy_cnv", "miran_mutation_cnv", "rna_mirna_methy","rna_mirna_mutation","rna_mirna_cnv",
            "rna_methy_mutation","rna_methy_cnv", "rna_mutation_cnv",
            "mirna_methy_mutation", "methy_mutation_cnv",
            
            "rna_mirna_methy_mutation","rna_mirna_methy_cnv", "rna_mirna_mutation_cnv", 
            "rna_methy_mutation_cnv", "mirna_methy_mutation_cnv",
            
            "rna_mirna_methy_mutation_cnv")

###### prepare the data of heatmap 
a <- rep(mnames,5)

b <- c(rep("rna", 31),rep("mirna", 31),rep("methy", 31),rep("mutation", 31),rep("cnv", 31))

c <- c("1","0","0","0","0",
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

d <- data.frame(as.factor(a),as.factor(b),as.character(c))
colnames(d) <- c("a","b","c")

###### non parameter methods #####
#### the cindex of block forest
resultscindex_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_bf)
colnames(resultscindex_bf) <- c("comb", "dat","cindex_bf")

resultswide <- reshape(resultscindex_bf , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_bf.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
#means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  ggtitle('bf')+
  geom_boxplot() + 
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')

## draw figure
pall_cindex_bf <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.25, 1) )


#### the ibrier of block forest
resultsibrier_bf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_bf)
colnames(resultsibrier_bf) <- c("comb", "dat","ibrier_bf")

resultswide <- reshape(resultsibrier_bf , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_bf.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  ggtitle('bf')+
  geom_boxplot() + 
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_ibrier_bf <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.25, 1) )




#### the cindex of random forest 
resultscindex_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_rf)
colnames(resultscindex_rf) <- c("comb", "dat","cindex_rf")

resultswide <- reshape(resultscindex_rf , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_rf.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
#means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('rsf')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_cindex_rf <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.25, 1) )


#### the ibrier of random forest 
resultsibrier_rf <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_rf)
colnames(resultsibrier_rf) <- c("comb", "dat","ibrier_rf")

resultswide <- reshape(resultsibrier_rf , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_rf.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('rsf')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_ibrier_rf <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1,-0.25, 1) )



###### parameter methods #####
#### the cindex of lasso
resultscindex_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_lasso)
colnames(resultscindex_lasso) <- c("comb", "dat","cindex_lasso")

resultswide <- reshape(resultscindex_lasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_lasso.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('lasso')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_cindex_lasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.32, 1) )


#### the ibrier of lasso
resultsibrier_lasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_lasso)
colnames(resultsibrier_lasso) <- c("comb", "dat","ibrier_lasso")

resultswide <- reshape(resultsibrier_lasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_lasso.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
#means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('lasso')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_ibrier_lasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.32, 1) )



#### the cindex of prioritylasso
resultscindex_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_prioritylasso)
colnames(resultscindex_prioritylasso) <- c("comb", "dat","cindex_prioritylasso")

resultswide <- reshape(resultscindex_prioritylasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_prioritylasso.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() +
  ggtitle('prioritylasso')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_cindex_prioritylasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.32, 1) )



#### the ibrier of prioritylasso
resultsibrier_prioritylasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_prioritylasso)
colnames(resultsibrier_prioritylasso) <- c("comb", "dat","ibrier_prioritylasso")

resultswide <- reshape(resultsibrier_prioritylasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_prioritylasso.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('prioritylasso')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_ibrier_prioritylasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.32, 1) )




#### the cindex of ipflasso
resultscindex_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$cindex_ipflasso)
colnames(resultscindex_ipflasso) <- c("comb", "dat","cindex_ipflasso")

resultswide <- reshape(resultscindex_ipflasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(-x)))

colnames(resultranks) <- gsub("cindex_ipflasso.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('ipflasso')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_cindex_ipflasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.25, 1) )



#### the ibrier of ipflasso
resultsibrier_ipflasso <- cbind(resultsumsum[,c(1,2)],resultsumsum$ibrier_ipflasso)
colnames(resultsibrier_ipflasso) <- c("comb", "dat","ibrier_ipflasso")

resultswide <- reshape(resultsibrier_ipflasso , idvar = c("dat"), timevar = "comb", direction = "wide")
resultswide1 <- resultswide[,-1]

# Calculate ranks of the methods:
resultranks <- t(apply(resultswide1, 1, function(x) rank(x)))

colnames(resultranks) <- gsub("ibrier_ipflasso.", "", colnames(resultranks))

resultrankstemp <- data.frame(resultranks)


resultranks2 <- reshape(resultrankstemp, varying=colnames(resultranks), 
                        v.names="rank", 
                        timevar="combin", times=colnames(resultranks),
                        direction="long")
resultranks2$combin <- factor(resultranks2$combin, levels=mnames)
set.seed(1234)

means <- aggregate(rank ~  combin, resultranks2, mean)
means <- means[order(means$rank),]
means$combin

p1 <- ggplot(data=resultranks2, aes(x=combin, y=rank)) + 
  theme_bw() + 
  geom_boxplot() + 
  ggtitle('ipflasso')+
  scale_x_discrete(limits=means$combin) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = " ", y = "Data set specific ranks")

p2<-ggplot(d, aes(x=a, y=b, fill=c)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_manual(values = c("grey70","darkviolet"))+
  labs(x = " ", y = "") +
  coord_fixed() + 
  scale_y_discrete(limits=c("mutation", "methy","cnv","mirna","rna"),
                   labels= c("mut", "met","cnv","mirna","rna")) +
  scale_x_discrete(limits=means$combin) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')

## draw figure
pall_ibrier_ipflasso <- ggarrange(p1, NULL, p2,  nrow = 3, align="v",heights = c(1, -0.25, 1) )


#### combine figures ####
p3_cindex <- ggarrange(pall_cindex_rf,  NULL, pall_cindex_bf,  NULL, pall_cindex_ipflasso,
                             nrow = 5, align="v", heights = c(1, -0.14, 1,-0.14,1))
ggsave(file="C:/Users/yingxiali/Desktop/paper3/3_rcode/Figures/tables and figures/p3_cindex_nov.png", 
       p3_cindex, width=8, height=14)


p2_cindex <- ggarrange(pall_cindex_lasso,  NULL, pall_cindex_prioritylasso,
                       nrow = 3, align="v", heights = c(1, -0.2, 1))
ggsave(file="C:/Users/yingxiali/Desktop/paper3/3_rcode/Figures/tables and figures/p2_cindex_nov.png", 
       p2_cindex, width=8, height=12)


p3_ibrier <- ggarrange(pall_ibrier_rf,  NULL, pall_ibrier_bf,  NULL, pall_ibrier_ipflasso,
                       nrow = 5, align="v", heights = c(1, -0.14, 1,-0.14,1))
ggsave(file="C:/Users/yingxiali/Desktop/paper3/3_rcode/Figures/tables and figures/p3_ibrier_nov.png", 
       p3_ibrier, width=8, height=14)


p2_ibrier <- ggarrange(pall_ibrier_lasso,  NULL, pall_ibrier_prioritylasso,
                       nrow = 3, align="v", heights = c(1, -0.2, 1))
ggsave(file="C:/Users/yingxiali/Desktop/paper3/3_rcode/Figures/tables and figures/p2_ibrier_nov.png", 
       p2_ibrier, width=8, height=12)







