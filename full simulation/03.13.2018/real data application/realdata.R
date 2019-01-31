# dim(fmridata)=943, 26
Whole_fmridata <- read.table("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Given/FDA_data.txt",
                      sep = "", header = FALSE, fill = T)
colnames(Whole_fmridata)=c('ID','temp','rating',paste("score",1:23,sep=""))
#center_whole_fmri <- cbind(Whole_fmridata[,1:3], Whole_fmridata[,4:26]-rowMeans(Whole_fmridata[,4:26])) #center each row 
center_whole_fmri <- cbind(Whole_fmridata[,1:3], Whole_fmridata[,4:26]-as.vector(rep(1, 19803)) %*% t(colMeans( Whole_fmridata[,4:26]))) #center each column 

pvalue_21 <- c()
test_21 <- list()
indexfortest=1
for (i in 1:21){
  fmridata <- (center_whole_fmri)[(1+943*(i-1)):(943*i),]
  M <- fmridata[,4:26]
  D <- 23
  t <- (1:D)/D
  knots <- 7 # previous setting 10
  p <- 5  # previous setting p <- 7, the number of degree for B-splines we use, p=order-1
  percent <- .9 # npc=8, percent=0.9 then npc=7; percent=0.99, npc=9
  
  
  # FPCA on M
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  library(MASS)
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
  
  ID.uni <- c()
  index <- matrix(1:(nSubj*npc), nrow = npc, ncol = nSubj)
  for (i in 1:length(nRep)){
    ID.uni = c(ID.uni, c(index[,i], rep(0, nRep[i] - npc)))
  }
  
  dummyX <- fmridata$temp
  dummyX <- cbind(dummyX, -dummyX + 1) #hot, warm
  Y <- fmridata$rating
  
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
      dummyX[1:nRep[k], ] <- u.tra %*% dummyX[1:nRep[k], ]
      ascore[1:nRep[k], ] <- rbind(u.tra[1:npc, ] %*% ascore[1:nRep[k], ], 
                                   matrix(0, nrow = nRep[k] - npc, 
                                          ncol = npc))
    }else{
      Y[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])] <- u.tra %*% Y[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])]
      dummyX[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ] <- u.tra %*% dummyX[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ]
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
                             temp.1 = dummyX[, 1],
                             temp.2 = dummyX[, 2],
                             ID = as.factor(ID),
                             ID.uni = as.factor(ID.uni),
                             ascore = ascore,
                             z.sim.uni = z.sim.uni)
  
  
  # 'lmer' model
  designMatrix.lmm <- designMatrix
  
  additive0.sim <- paste(1:npc, collapse = " + ascore.")
  additive.sim <- paste(1:npc, collapse = " | ID) + (0 + ascore.")
  
  model.sim <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                additive0.sim, 
                                " + (0 + temp.1 | ID) + (0 + temp.2 | ID) + (0 + z.sim.uni | ID.uni)", 
                                sep = ""))
  fullReml <- lmer(model.sim, data = designMatrix.lmm)
  f.slope <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
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
  indexfortest = indexfortest + 1
}

save(pvalue_21, test_21, center_whole_fmri, file="/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/03.13.2018/real data application/realdata.Rdata")

# find the most significant 3 p-values
small <- order(pvalue_21)[1:3] #1 10 19


  
  
  
  