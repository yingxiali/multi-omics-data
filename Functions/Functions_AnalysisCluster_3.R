########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data")

########################################################


evaluatesetting <- function(iter) {
  
  library("dplyr")
  
  library("survival")
  
  library("blockForest")
  
  library("glmnet")
  library("prioritylasso")
  library("ipflasso")
  
  #library("carSurv") # for feature selection
  library("ranger")# for feature selection
  library("prodlim")
  library("survcomp") # for cindex
  library("pec") # for ibrier
  
  # Initiate lists in which the results will be stored:
  
  ibrier_rf <- list()
  cindex_rf <- list()
  
  ibrier_ipflasso <- list()
  cindex_ipflasso <- list()
  
  ibrier_lasso <- list()
  cindex_lasso <- list()
  
  ibrier_bf <- list()
  cindex_bf <- list()
   paramvalues_bf <- list()
  
   ibrier_prioritylasso <- list()
   cindex_prioritylasso <- list() 
   pf_prioritylasso <- list()
  
  # Obtain information for the iter-th setting:
  
  dat <- scenariogrid$dat[iter]
  seed <- scenariogrid$seed[iter]
  comb <- scenariogrid$comb[iter]
  
  cvind <- scenariogrid$cvind[iter]
  cvfoldind <- scenariogrid$cvfoldind[iter]
  
  
  # Load data set:
  
  load(paste("./Data/", dat, sep=""))
  #load(paste("./data/5_omics_data/", dat, sep=""))
  
  # Make covariate matrix and target variable:
  
  if(comb=="rna_mirna_methy") {
    
    X <- cbind(clindata, rnadata, mirnadata, methydata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(rnadata),ncol(mirnadata),ncol(methydata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="rna_mirna_mutation"){
    
    X <- cbind(clindata, rnadata, mirnadata, mutationdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(rnadata),ncol(mirnadata),ncol(mutationdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="rna_mirna_cnv"){
    X <- cbind(clindata, rnadata, mirnadata, cnvdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(rnadata),ncol(mirnadata),ncol(cnvdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="rna_methy_mutation"){
    
    X <- cbind(clindata, rnadata, methydata, mutationdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(rnadata),ncol(methydata),ncol(mutationdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="rna_methy_cnv"){
    
    X <- cbind(clindata, rnadata, methydata, cnvdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(rnadata),ncol(methydata),ncol(cnvdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="rna_mutation_cnv"){
    
    X <- cbind(clindata,rnadata, mutationdata, cnvdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(rnadata), ncol(mutationdata),ncol(cnvdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="mirna_methy_mutation"){     
    
    X <- cbind(clindata, mirnadata, methydata,mutationdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(mirnadata),ncol(methydata),ncol(mutationdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="miran_methy_cnv"){
    
    X <- cbind(clindata, mirnadata,methydata, cnvdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(mirnadata),ncol(methydata),ncol(cnvdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="miran_mutation_cnv"){
    
    X <- cbind(clindata, mirnadata, mutationdata, cnvdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(mirnadata), ncol(mutationdata), ncol(cnvdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } else if(comb=="methy_mutation_cnv"){
    
    X <- cbind(clindata, methydata, mutationdata, cnvdata)
    
    block <- rep(1:4, times=c(ncol(clindata), ncol(methydata), ncol(mutationdata), ncol(cnvdata)))
    
    block <- lapply(1:4, function(x) which(block==x))
    
    y <- cbind(surdata$time, surdata$status)
    
  } 
  
  
  # Divide data set into training and test data:
  
  ncv <- 5
  set.seed(seed)
  
  cvdiv <- makeCVdiv(n=nrow(X), ncv=ncv)
  
  Xtrain <- X[cvdiv!=cvfoldind,]
  ytrain <- y[cvdiv!=cvfoldind,]
  
  Xtest <- X[cvdiv==cvfoldind,]
  ytest <- y[cvdiv==cvfoldind,]
  
  # Feature selection per block using ranger
  indsel <- c()
  blocksub <- c()
  
  for(i in seq(along=block)) {
    if(length(block[[i]]) <= 2500) {
      indsel <- c(indsel, block[[i]])
      blocksub <- c(blocksub, rep(i, length(block[[i]])))
    }
    else {
      varinds <- block[[i]]
      Xtrain_s = Xtrain[,varinds]
      Xtrain_s <- Xtrain_s[, apply(Xtrain_s, 2, sum) != 0]
      Xtrain_s$survivaltime <- ytrain[,1]
      Xtrain_s$statusvariable <- ytrain[,2]
      varimp <- ranger(data=Xtrain_s, dependent.variable.name = "survivaltime", status.variable.name = "statusvariable", 
                       importance="permutation", splitrule="extratrees", num.trees=10000, num.threads=1)$variable.importance
      inds <- which(colnames(Xtrain) %in% colnames( Xtrain_s[,order(varimp,decreasing = T)[1:2500]]))
      indsel <- c(indsel, inds)
      blocksub <- c(blocksub, rep(i, 2500))
      
    }
  }
  
  blocksub <- blocksub[order(indsel)]
  indsel <- indsel[order(indsel)]
  
  blocksub <- lapply(1:length(block), function(x) which(blocksub==x))
  
  Xtrain <- Xtrain[,indsel]
  Xtest <- Xtest[,indsel]
  
  ##reduce memory
  rm(methydata)
  gc()
 
  
  # blockforest
  
  bf <- blockforestwrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, blocksub=blocksub)
  ibrier_bf[[1]] <- bf$ibrier
  cindex_bf[[1]] <- bf$cindex
  paramvalues_bf[[1]] <- bf$paramvalues
  
  # prioritylasso:
  
  prioritylasso <- prioritylassowrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, blocksub=blocksub)
  ibrier_prioritylasso[[1]] <- prioritylasso$ibrier
  cindex_prioritylasso[[1]] <- prioritylasso$cindex
  pf_prioritylasso[[1]] <- prioritylasso$pf
  
  # random forest
  
  rf <- rangerwrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, blocksub=blocksub)
  ibrier_rf[[1]] <- rf$ibrier
  cindex_rf[[1]] <- rf$cindex
  
  # ipflasso
  
  ipflasso <- ipflassowrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, blocksub=blocksub)
  ibrier_ipflasso[[1]] <- ipflasso$ibrier
  cindex_ipflasso[[1]] <- ipflasso$cindex
  
  # lasso
  
  lasso <- lassowrap(Xtrain=Xtrain, ytrain=ytrain, Xtest=Xtest, ytest=ytest, blocksub=blocksub)
  ibrier_lasso[[1]] <- lasso$ibrier
  cindex_lasso[[1]] <- lasso$cindex
  
  
  # Combine results in list:
  
  res <- list( 
    
    ibrier_rf = ibrier_rf,
    cindex_rf = cindex_rf,
    
    ibrier_ipflasso = ibrier_ipflasso,
    cindex_ipflasso = cindex_ipflasso,
    
    ibrier_lasso = ibrier_lasso,
    cindex_lasso = cindex_lasso,
    
    
    ibrier_bf = ibrier_bf,
    cindex_bf = cindex_bf,
    paramvalues_bf = paramvalues_bf,
    
    ibrier_prioritylasso = ibrier_prioritylasso,
    cindex_prioritylasso = cindex_prioritylasso,
    pf_prioritylasso = pf_prioritylasso,
    
    settingind=iter)
  
  
  # Save results in folder:
  save(res, file=paste("./Results/rda_files/res", iter, ".Rda", sep=""))
  
}

# Function to generate the splittings for cross-validation:

# Input parameters:

# n    - number of observations in the data set
# ncv  - number of folds to use

makeCVdiv <- function(n, ncv) {
  
  nperfold <- ceiling(n/ncv)
  partition <- sample(rep(1:ncv, nperfold), n)
  
  partition
  
}
# pec_blockforest
predictSurvProb.blockForest <- function (object, newdata, times, ...) {
  ptemp <- blockForest:::predict.blockForest(object, data = newdata[,-c(1:2)], block.method=block.method, num.threads=1)$survival
  pos <- prodlim::sindex(jump.times = object$unique.death.times, 
                         eval.times = times)
  p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  p
}

# pec_randomforest
predictSurvProb.ranger <- function (object, newdata, times, ...) {
  ptemp <- ranger:::predict.ranger(object, data = newdata[,-c(1:2)], importance = "none")$survival
  pos <- prodlim::sindex(jump.times = object$unique.death.times, 
                         eval.times = times)
  p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  p
}

#pec_ipflasso, the code comes from https://github.com/HerrMo/multi-omics_benchmark_study/blob/master/R/ancillary_code_bench.R 
compute_coxph_mod <- function(coeffs, index_nz, data, ...) {
  index_nz_y <- c(TRUE, TRUE, index_nz)
  sub_data <- data[, index_nz_y, drop = FALSE]
  coxph(
    Surv(time, status) ~ .,
    data = sub_data,
    init = coeffs[index_nz],
    iter.max = 0,
    x = TRUE
  )
}

check_coeffs <- function(coeffs) {
  # check for NA coeffs and set to zero 
  if (any(is.na(coeffs))) {
    coeffs[is.na(coeffs)] = 0
  }
  coeffs
}

make_nz <- function(coeffs, check_all_zero = TRUE) {
  index_nz <- coeffs != 0
  
  if (check_all_zero) {
    if (!any(index_nz)) {
      stop("Null model fitted.")
    }
  }
  
  index_nz
}

get_nz <- function(coeffs, check_all_zero = TRUE) {
  coeffs %>%
    check_coeffs() %>%
    make_nz()
}


# Input parameters:

# Xtrain        - covariate matrix of the training data set
# ytrain        - target variable of the training data set. Matrix with two 
#                 columns, where the first column contains the vector of 
#                 survival/censoring times (one for each observation) and the 
#                 second column contains the status variable, that has the 
#                 value '1' if the corresponding time is a survival time and 
#                 '0' if that time is a censoring time.
# Xtest         - covariate matrix of the test data set
# ytest         - target variable in the test data set. Matrix with two
#                 columns, where the first column contains the vector of 
#                 survival/censoring times (one for each observation) and the 
#                 second column contains the status variable, that has the 
#                 value '1' if the corresponding time is a survival time and 
#                 '0' if that time is a censoring time.
# blocksub         - A list of length equal to the number of blocksub considered. 
#                 Each entry contains the vector of column indices in 'Xtrain' 
#                 of the covariates in one of the blocks.


#ipflasso
ipflassowrap <- function(Xtrain, ytrain, Xtest, ytest, blocksub) {
  
  colnames(ytrain) <- c("time","status")
  colnames(ytest) <- c("time","status")
  
  traindata <- cbind(ytrain, Xtrain)
  testdata <- cbind(ytest, Xtest)
  
  
  #ytrain[,1][ytrain[,1] <= 0] <- min(ytrain[,1][ytrain[,1]>0])
  
  
  # Train IPF-LASSO 
  # the parameters were set from Herrmann's code
  cvres <- cvr.adaptive.ipflasso(X=as.matrix(traindata[,-c(1:2)]),Y=as.matrix(traindata[,c(1:2)]),
                                 family="cox", type.measure = "deviance", standardize = TRUE,
                                 alpha = 0, type.step1 = "sep", nfolds=5,ncv=10,
                                 blocks= blocksub)
  # calculate cindex
  # Obtain predictions on test data:
  riskpred <- ipflasso.predict(cvres, Xtest = as.matrix(testdata[,-c(1:2)]))$linpredtest[,1]
  cindex <- concordance.index(x=riskpred, surv.time=testdata[,1], surv.event=testdata[,2])$c.index
  
  # calculate ibreir
  mod <- cvres
  coeffs <- mod$coeff[-1, mod$ind.bestlambda] # -1 to remove intercept (which is NA for Cox!)
  index_nz = get_nz(coeffs)
  cph <- compute_coxph_mod(coeffs, index_nz, traindata)
  mod$cph <- cph
  mod$index_nz <- index_nz
  
  index_nz_y1 <- c(F, F, mod$index_nz)
  probs = pec::predictSurvProb(mod$cph,
                               newdata = testdata[, index_nz_y1, drop = FALSE],
                               times = unique(traindata[,"time"]))
  
  # A formula to be inputted into the pec command
  frm <- as.formula(paste("Surv(time, status)~",
                          paste(colnames(traindata[,-c(1:2)]), collapse="+")))
  
  PredError <- pec(object=mod$cph,
                   formula = frm, cens.model="marginal",
                   data=testdata, verbose=F, exact = T)
  ibrier <- crps(PredError,times=max(testdata[,1]),start=0)[2,]
  
  # Return results:
  
  return(list(cindex=cindex, ibrier=ibrier))
  
}

#####
#lasso
lassowrap <- function(Xtrain, ytrain, Xtest, ytest, blocksub) {
  
  
  colnames(ytrain) <- c("time","status")
  colnames(ytest) <- c("time","status")
  
  traindata <- cbind(ytrain, Xtrain)
  testdata <- cbind(ytest, Xtest)
  
  
  #ytrain[,1][ytrain[,1] <= 0] <- min(ytrain[,1][ytrain[,1]>0])
  
  
  # TrainLasso
  
  
  lambdaopt <- cv.glmnet(x=as.matrix(Xtrain), y=Surv(ytrain[,1], ytrain[,2]), family="cox", alpha=1,
                         penalty.factor =c(rep(0,length(unlist(blocksub[1]))),rep(1,length(unlist(blocksub[-1])))))$lambda.min
  fittedmodel <- glmnet(x=as.matrix(Xtrain), y=Surv(ytrain[,1], ytrain[,2]), lambda = lambdaopt, family="cox", alpha=1,
                        penalty.factor =c(rep(0,length(unlist(blocksub[1]))),rep(1,length(unlist(blocksub[-1]))))) 
  
  
  # calculate cindex
  # Obtain predictions on test data:
  riskpred <- predict(fittedmodel, newx = as.matrix(testdata[,-c(1:2)]))[,1]
  cindex <- concordance.index(x=riskpred, surv.time=testdata[,1], surv.event=testdata[,2])$c.index
  
  # calculate ibreir
  mod <- fittedmodel
  coeffs <- as.vector(coef(mod, s = mod$lambda.min))
  index_nz = get_nz(coeffs)
  cph <- compute_coxph_mod(coeffs, index_nz, traindata)
  mod$cph <- cph
  mod$index_nz <- index_nz
  
  index_nz_y1 <- c(F, F, mod$index_nz)
  probs = pec::predictSurvProb(mod$cph,
                               newdata = testdata[, index_nz_y1, drop = FALSE],
                               times = unique(traindata[,"time"]))
  
  # A formula to be inputted into the pec command
  frm <- as.formula(paste("Surv(time, status)~",
                          paste(colnames(traindata[,-c(1:2)]), collapse="+")))
  
  PredError <- pec(object=mod$cph,
                   formula = frm, cens.model="marginal",
                   data=testdata, verbose=F, exact = T)
  ibrier <- crps(PredError,times=max(testdata[,1]),start=0)[2,]
  
  # Return results:
  
  return(list(cindex=cindex, ibrier=ibrier, pf=pf))
  
}

#random forest
rangerwrap <- function(Xtrain, ytrain, Xtest, ytest, blocksub) {
  
  
  colnames(ytrain) <- c("time","status")
  colnames(ytest) <- c("time","status")
  
  traindata <- cbind(ytrain, Xtrain)
  testdata <- cbind(ytest, Xtest)
  
  
  # Train random forest:
  
  rangerobj <- ranger(Surv(time, status) ~ ., data = traindata,  splitrule="extratrees",replace = FALSE,
                      probability = FALSE,num.threads=1, 
                      always.split.variables = c(colnames(Xtrain[,unlist( blocksub[1])])))
  
  # calculate cindex
  
  # Obtain predictions on test data
  riskpred <- rowSums(predict(rangerobj, data=testdata[,-c(1:2)], num.threads=1)$chf)
  cindex <- concordance.index(x=riskpred, surv.time=testdata[,1], surv.event=testdata[,2])$c.index
  
  # calculate ibrier
  
  # A formula to be inputted into the pec command
  frm <- as.formula(paste("Surv(time, status)~",
                          paste(colnames(traindata[,-c(1:2)]), collapse="+")))
  #note: colnames of testdata[,1:2] must be "time"and"status"
  PredError <- pec(object=rangerobj,
                   formula = frm, cens.model="marginal",
                   data=testdata, verbose=F, exact = T)
  ibrier <- crps(PredError,times=max(testdata[,1]),start=0)[2,]
  
  # Return results:
  
  return(list(cindex=cindex, ibrier=ibrier))
  
}

blockforestwrap <- function(Xtrain, ytrain, Xtest, ytest, blocksub) {
  
  
  colnames(ytrain) <- c("time","status")
  colnames(ytest) <- c("time","status")
  
  traindata <- cbind(ytrain, Xtrain)
  testdata <- cbind(ytest, Xtest)
  
  
  # Train blockforest:
  
  blockforobj <- blockfor(traindata[,-c(1:2)], as.matrix(traindata[,c(1:2)]), num.trees = 2000, replace = FALSE, probability = FALSE, blocks=blocksub,
                          nsets = 300, num.trees.pre = 1500, splitrule="extratrees", 
                          block.method = "BlockForest", num.threads=1, always.select.block=1)
  
  paramvalues <- blockforobj$paramvalues
  
  # calculate cindex
  riskpred <- rowSums(predict(blockforobj$forest, data=testdata[,-c(1:2)], block.method= "BlockForest", num.threads=1)$chf)
  cindex <- concordance.index(x=riskpred, surv.time=testdata[,1], surv.event=testdata[,2])$c.index
  
  # calculate ibrier
  # A formula to be inputted into the pec command
  frm <- as.formula(paste("Surv(time, status)~",
                          paste(colnames(traindata[,-c(1:2)]), collapse="+")))
  #note: colnames of testdata[,1:2] must be "time"and"status"
  PredError <- pec(object=blockforobj$forest,
                   formula = frm, cens.model="marginal",
                   data=testdata, verbose=F, exact = T)
  ibrier <- crps(PredError,times=max(testdata[,1]),start=0)[2,]
  
  # Return results:
  
  return(list(cindex=cindex, ibrier=ibrier, paramvalues=paramvalues))
  
}

prioritylassowrap <- function(Xtrain, ytrain, Xtest, ytest, blocksub) {
  
  
  colnames(ytrain) <- c("time","status")
  colnames(ytest) <- c("time","status")
  
  traindata <- cbind(ytrain, Xtrain)
  testdata <- cbind(ytest, Xtest)
  
  # Train priority-Lasso (including the training of the penalty factors):
  means <- c(rep(0,length(blocksub)))
  
  for(i in seq(along=blocksub))
    means[i] <- mean(abs(coef(cv.glmnet(x=as.matrix(Xtrain[,blocksub[[i]]]), y=Surv(ytrain[,1], ytrain[,2]), 
                                        family="cox", alpha=0), s="lambda.min", standardize = TRUE)))
  
  means_check <- 1/means
  
  badblocks <- which(is.infinite(means_check))
  if(length(badblocks) > 0) {
    blocksub <- blocksub[-badblocks]
    Xtrain <- Xtrain[,unlist(blocksub)]
    Xtest <- Xtest[,unlist(blocksub)]
    means_check <- means_check[-badblocks]
    
    cumsums <- cumsum(sapply(blocksub, length))
    blocksub <- mapply(function(x, y) x:y, c(1, cumsums[-length(cumsums)]+1), cumsums, SIMPLIFY=FALSE)
  }
  
  pf <- means_check 
  
  # the parameters were set from Herrmann's code
  cvres <- prioritylasso(X=as.matrix(Xtrain), Y=Surv(ytrain[,1], ytrain[,2]), family="cox",
                         type.measure="deviance", block1.penalization=FALSE, blocks=blocksub[c(1, order(pf[-1])+1)])
  # calculate cindex
  riskpred <- predict(cvres, newdata = as.matrix(testdata[,-c(1:2)]), type="link")[,1]
  cindex <- concordance.index(x=riskpred, surv.time=testdata[,1], surv.event=testdata[,2])$c.index
  
  # calculate ibreir
  mod <- cvres
  coeffs <- mod$coeff 
  index_nz = get_nz(coeffs)
  cph <- compute_coxph_mod(coeffs, index_nz, traindata)
  mod$cph <- cph
  mod$index_nz <- index_nz
  
  index_nz_y1 <- c(F, F, mod$index_nz)
  probs = pec::predictSurvProb(mod$cph,
                               newdata = testdata[, index_nz_y1, drop = FALSE],
                               times = unique(traindata[,"time"]))
  
  # A formula to be inputted into the pec command
  frm <- as.formula(paste("Surv(time, status)~",
                          paste(colnames(traindata[,-c(1:2)]), collapse="+")))
  
  PredError <- pec(object=mod$cph,
                   formula = frm, cens.model="marginal",
                   data=testdata, verbose=F, exact = T)
  ibrier <- crps(PredError,times=max(testdata[,1]),start=0)[2,]
  
  # Return results:
  
  return(list(cindex=cindex, ibrier=ibrier, pf=pf))
  
}

