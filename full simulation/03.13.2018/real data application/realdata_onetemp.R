# dim(fmridata)=943, 26
Whole_fmridata <- read.table("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Given/FDA_data.txt",
                             sep = "", header = FALSE, fill = T)
colnames(Whole_fmridata)=c('ID','temp','rating',paste("score",1:23,sep=""))
center_whole_fmri <- cbind(Whole_fmridata[,1:3], Whole_fmridata[,4:26]-as.vector(rep(1, 19803)) %*% t(colMeans( Whole_fmridata[,4:26]))) #center each column 

pvalue_21 <- c()  # pvalues from equal-variance test
test_21 <- list() # test details from equal-variance test for each region
pvalue_21.bonf <- list() # save each bonf adjusted p-values for each region
result_21.bonf <- c() # save reject(True)/accept(False) result from bonf test

indexfortest=1
library(refund)
library(lme4)
library(nlme)
library(arm)
library(RLRsim)
library(MASS)

for (i in 1:21){
  fmridata <- (Whole_fmridata)[(1+943*(i-1)):(943*i),]
  fmridata <- cbind(fmridata[,1:3,], fmridata[,4:26]-as.vector(rep(1, 943)) %*% t(colMeans( fmridata[,4:26]))) #center
  
  M <- fmridata[,4:26]
  D <- 23
  t <- (1:D)/D
  knots <- 7 # previous setting 10
  p <- 5  # previous setting p <- 7, the number of degree for B-splines we use, p=order-1
  percent <- .9 # npc=8, percent=0.9 then npc=7; percent=0.99, npc=9
  
  
  # FPCA on M
  results <- fpca.face(data.matrix(M), center = TRUE, argvals = t, knots = knots, pve = percent, p = p, lambda = 0) # pve need to be chosen!
  npc <- results$npc
  lambda <- results$evalues/D
  score <- results$scores
  ascore <- score[, 1:npc]/sqrt(D)
  efunctions <- results$efunctions*sqrt(D)
  #par(mfrow=c(4,2), mar=c(2.9,3,1,1), mgp = c(2,1,0), cex.main = 1, cex.axis = 0.8, xaxs = 'r')
  #myat=seq(0, 1, by=0.2)
  #for (i in 1:npc){
  #  plot(t, efunctions[,i],type="l",
  #       xlab="t",
  #       ylab="",
  #       xaxt="n",
  #       xlim=c(0.0,1.0),
  #       main=paste("eigenfunction ",i,sep=""),lwd=2, ylim=c(min(efunctions)-0.5,max(efunctions)+0.5))
  #  axis(1, at = myat, labels = round(myat,digits=1))
  #}
  
  nSubj=length(unique(fmridata$ID))
  nRep = rep(1, nSubj) #nRep in c(39, 48)
  for (i in 1:nSubj){
    nRep[i] = length(which(fmridata$ID==i))
  }
  
  ###########################################################################
  ID <- fmridata$ID
  dummyX <- fmridata$temp
  Y <- fmridata$rating
  
  designMatrix.reg <- data.frame(rating = Y, 
                             temp.1 = dummyX,
                             ID = as.factor(ID),
                             ascore = ascore)
  
  ###########################################################################
  ID.uni <- c()
  index <- matrix(1:(nSubj*npc), nrow = npc, ncol = nSubj)
  for (i in 1:length(nRep)){
    ID.uni = c(ID.uni, c(index[,i], rep(0, nRep[i] - npc)))
  }
  
  z.sim.uni = c()
  for(k in 1:nSubj){
    if(k==1){
      svd <- svd(ascore[1:nRep[1], ] %*% t(ascore[1:nRep[1], ])) #SVD on A_i
    }else{
      svd <- svd(ascore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ] %*% t(ascore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ])) #SVD on A_i
    }
    u.tra <- t(svd$v)
    u <- svd$u
    d <- (svd$d)[1:npc]
    if(k==1){
      Y[1:nRep[k]] <- u.tra %*% Y[1:nRep[k]]
      dummyX[1:nRep[k]] <- u.tra %*% dummyX[1:nRep[k]]
      ascore[1:nRep[k], ] <- rbind(u.tra[1:npc, ] %*% ascore[1:nRep[k], ], 
                                   matrix(0, nrow = nRep[k] - npc, 
                                          ncol = npc))
    }else{
      Y[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])] <- u.tra %*% Y[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])]
      dummyX[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])] <- u.tra %*% dummyX[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])]
      ascore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ] <- rbind(u.tra[1:npc, ] %*% ascore[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ], 
                                                               matrix(0, 
                                                                      nrow = nRep[k] - npc, 
                                                                      ncol = npc))
    }
    z.sim.uni <- c(z.sim.uni, sqrt(d), rep(0, nRep[k] - npc))
  }
  
  ###########################################################################
  ID <- fmridata$ID
  designMatrix <- data.frame(rating = Y, 
                             temp.1 = dummyX,
                             ID = as.factor(ID),
                             ID.uni = as.factor(ID.uni),
                             ascore = ascore,
                             z.sim.uni = z.sim.uni)
  
  
  # 'lmer' model
  designMatrix.lmm <- designMatrix
  
  additive0.sim <- paste(1:npc, collapse = " + ascore.")
  additive.sim <- paste(1:npc, collapse = " | ID) + (0 + ascore.")
  
  model.sim <- as.formula(paste("rating ~ 1 + temp.1 + ascore.", 
                                additive0.sim, 
                                " + (0 + temp.1 | ID) + (0 + z.sim.uni | ID.uni)", 
                                sep = ""))
  fullReml <- lmer(model.sim, data = designMatrix.lmm)
  f.slope <- as.formula(paste("rating ~ 1 + temp.1 + ascore.", 
                              additive0.sim, 
                              " + (0 + z.sim.uni | ID.uni)", 
                              sep = ""))
  m.slope <- lmer(f.slope, data = designMatrix.lmm)
  f0 <- as.formula(" . ~ . - (0 + z.sim.uni | ID.uni)")
  m0 <- update(fullReml, f0)
  
  test <- exactRLRT(m.slope, fullReml, m0)
  pvalue <- test$p[1] 
  pvalue_21 <- c(pvalue_21, pvalue)
  test_21[[indexfortest]] <- test
  
  
  
  
  ###################################################################
  ## bonferroni test ##
  additive.heter <- paste0(" + (0 + ascore.", 1:npc, " | ID)", collapse = "")
  bonf.test <- list()
  for(i in 1:npc){
    ii <- paste("ascore.", i, sep = "")
    # f.slope only contains the random effect to be tested
    f.slope <- as.formula(paste("Y ~ 1 + temp.1 + ascore.", 
                                additive0.sim, " + (0 +", ii, " | ID)", 
                                sep = ""))
    m.slope <- lmer(f.slope, data = designMatrix.reg)
    # mA is the model under alternative
    mA <- as.formula(paste("Y ~ 1 + temp.1 + ascore.", 
                           additive0.sim, " + (0 +", ii, " | ID)", 
                           "+ (0 + temp.1 | ID)",
                           sep = ""))
    fullReml <- lmer(mA, data = designMatrix.reg)
    #m0 is model under the null
    f0 <- as.formula(paste(" . ~ . - (0 + ", ii, "| ID)"))
    m0 <- update(fullReml, f0)
    bonf.test[[i]] <- exactRLRT(m.slope, fullReml, m0)
  }
  multiTest <- sapply(bonf.test, function(x) {
    c(statistic = x$statistic[1],
      "p-value" = x$p[1])})
  # use bonferroni correctiong method to adjust p-value
  bonf.pvalue <- p.adjust(multiTest[2,], "bonferroni")
  
  result_21.bonf <- c(result_21.bonf, sum(bonf.pvalue<0.05)>=1) # result=False denotes accepting the null; True denote rejecting the null; 
  pvalue_21.bonf[[indexfortest]] <- bonf.pvalue
  
  ###################################################################
  indexfortest = indexfortest + 1
}

# using lambda=0, we have equal-variance test all signifiant results, but ROI2 is relatively not highly significant(p-value=0.0045);
# using lambda=0, we have bonf test insignificant at all ROIs
# using lambda=NULL, we have equal-variance test 20 out of 21 significant results, only insignificant at ROI 2
# using lambda=NULL, we have bonf test only significant at 11th ROI

save(pvalue_21, test_21, result_21.bonf, pvalue_21.bonf, file="/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/03.13.2018/real data application/realdata_onetemp_bonf_EV.Rdata")

# find the most significant 3 p-values
#small <- order(pvalue_21)[1:3] #1 10 19





