#load("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/03.13.2018/real data application/realdata.Rdata")
# 1th voxel 
Whole_fmridata <- read.table("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Given/FDA_data.txt",
                             sep = "", header = FALSE, fill = T)
colnames(Whole_fmridata)=c('ID','temp','rating',paste("score",1:23,sep=""))
#center_whole_fmri <- cbind(Whole_fmridata[,1:3], Whole_fmridata[,4:26]-as.vector(rep(1, 19803)) %*% t(colMeans( Whole_fmridata[,4:26]))) #center each column 


realanalysis_outofsample <- function(voxel, Whole_fmridata){
  fmridata <- (Whole_fmridata)[(1+943*(voxel-1)):(943*voxel),]
  nSubj=length(unique(fmridata$ID))
  nRep = rep(1, nSubj) #nRep in c(39, 48)
  for (i in 1:nSubj){
    nRep[i] = length(which(fmridata$ID==i))
  }
  
  traindata <-  do.call(rbind, lapply(1:nSubj, function(x) fmridata[which(fmridata$ID==x),][1:(nRep[x]-10),]))
  traindata <- cbind(traindata[,1:3], traindata[,4:26]-as.vector(rep(1, dim(traindata)[1])) %*% t(colMeans(traindata[,4:26]))) #center
  
  testdata <- do.call(rbind, lapply(1:nSubj, function(x) fmridata[which(fmridata$ID==x),][(nRep[x]-10+1):nRep[x],]))
  
  #################################################################################
  
  M <- traindata[,4:26]
  D <- 23
  t <- (1:D)/D
  knots <- 7 # previous setting 10
  p <- 5  # previous setting p <- 7, the number of degree for B-splines we use, p=order-1
  m <- 3
  percent <- .95 # npc=8, percent=0.9 then npc=7; percent=0.99, npc=9
  
  #################################################################################
  
  # FPCA on M
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  library(MASS)
  results <- fpca.face(data.matrix(M), center = TRUE, argvals = t, knots = knots, pve = percent, p = p, m = m, lambda = 0) # pve need to be chosen!
  npc <- results$npc # voxel=1 then npc=9
  lambda <- results$evalues/D
  score <- results$scores
  ascore <- score[, 1:npc]/sqrt(D)
  efunctions <- results$efunctions*sqrt(D)
  
  #################################################################################
  
  # project testdata to get ascores
  M.test <- data.matrix(testdata[,4:26])
  xi <- (M.test - as.vector(rep(1, dim(M.test)[1])) %*% t(colMeans(M.test))) %*% efunctions/D
  # data for test
  designMatrix0 <- data.frame(rating = testdata$rating, 
                             temp.1 = testdata$temp,
                             ID = as.factor(testdata$ID),
                             ascore = xi)
  
  #################################################################################
  
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
  
  #################################################################################
  
  # data for train
  designMatrix1 <- data.frame(rating = traindata$rating, 
                              temp.1 = traindata$temp,
                              ID = as.factor(traindata$ID),
                              ascore = ascore)
  
  #################################################################################
  
  noRandom.simpart <- paste(1:npc, collapse = " + ascore.")
  noRandom.sim1 <- as.formula(paste("rating ~ 1 + temp.1 + ascore.", 
                                    noRandom.simpart, 
                                    " + (0 + temp.1 | ID) ", 
                                    sep = ""))
  NoRandomReml <- lmer(noRandom.sim1, data = designMatrix1)
  mseNoRandom <- mean((traindata$rating - predict(NoRandomReml))^2) #7365.821
  # npc=9
  Random.simpart <- paste0(" + (0 + ascore.", 1:npc,"|ID) ", collapse="")
  Random.sim1 <- as.formula(paste("rating ~ 1 + temp.1 + ascore.", 
                                  noRandom.simpart, 
                                  " + (0 + temp.1 | ID) ", 
                                  Random.simpart,
                                  sep = ""))
  RandomReml <- lmer(Random.sim1, data = designMatrix1)
  mseRandom <- mean((traindata$rating - predict(RandomReml))^2) #5913.915
  
  calculateRMSE_i <- function(Y, EST, designMatrix){
    out <- (Y- EST)^2
    mse=rep(0,length(unique(nSubj)))
    for(i in c(1:nSubj)){
      mse[i] <- mean(out[which(designMatrix$ID==i)])
    }
    return(list(rmse=sqrt(mean(mse)),mse=mse))
  }
  
  #################################################################################
  
  ## in sample: in the traindata
  RMSE_RandomReml <- calculateRMSE_i(traindata$rating, predict(RandomReml), designMatrix1) 
  RMSE_NoRandomReml <- calculateRMSE_i(traindata$rating, predict(NoRandomReml), designMatrix1)
  rsquare <- c(MuMIn::r.squaredGLMM(RandomReml)[,2], MuMIn::r.squaredGLMM(NoRandomReml)[,2])
  names(rsquare) <- c("FMM", "FLM")
  
  ## out of sample: in the testdata
  yhat.test.fmm <- predict(RandomReml, newdata=designMatrix0)
  yhat.test.flm <- predict(NoRandomReml, newdata=designMatrix0)
  RMSE.test.fmm <- calculateRMSE_i(testdata$rating, yhat.test.fmm, designMatrix0) 
  RMSE.test.flm <- calculateRMSE_i(testdata$rating, yhat.test.flm, designMatrix0) 
  
  #################################################################################
  
  fixeffWithTemp1 <- fixef(NoRandomReml)
  betaWithTemp1 <- efunctions %*% as.vector(fixeffWithTemp1[3:(npc+2)])
  
  ### RandomReml
  fixeffWithTemp2 <- fixef(RandomReml)
  randeffWithTemp_IDuni2 <- as.matrix(ranef(RandomReml)$ID[,2:(npc+1)])
  betaWithTemp2 <- efunctions %*% as.vector(fixeffWithTemp2[3:(npc+2)])
  
  theta_ik2 <- t(randeffWithTemp_IDuni2) #defalt is by column
  betai2 <- efunctions %*% theta_ik2  #23*8 x 8*20 = 23*20
  beta2 <- cbind(betaWithTemp2, apply(betai2, 2, function(x) x+betaWithTemp2))
  colnames(beta2) <- c("Fixed_effect", paste("Subject", 1:nSubj, sep = ""))
  
  ######  violin plot ######
  violin1 <- data.frame(cbind(c(RMSE.test.flm$mse, RMSE.test.fmm$mse),
                              rep(1:2, each=nSubj)))
  colnames(violin1) <- c("MSE","Model")
  violin1$Model <- as.factor(violin1$Model)
  
  library(ggplot2)
  MSEviolin <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/outsample_violin_MSE",voxel,".png")
  p1 <- ggplot(violin1, aes(x=Model, y=MSE)) +
    geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.alpha = 0.1) +
    geom_jitter( width = 0.2, alpha = 0.6) +
    geom_hline(aes(yintercept = median(RMSE.test.fmm$mse)), color="red", linetype="dashed", size = 0.8)+
    scale_x_discrete(breaks=c("1", "2"),
                     labels=c("FLM", "FMM")) +
    labs(title=paste("ROI ", voxel, "(out of sample)")) + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = rel(0.9), face="bold")) +
    theme(axis.title.y = element_text(size = rel(0.8))) +
    theme(axis.title.x = element_text(size = rel(0.8))) 
  cowplot::ggsave(filename=MSEviolin, p1)
  
  betaplot <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/outsample_betaplot_",voxel,".png")
  png(filename=betaplot)
  plot(t, beta2[,2], col="darkgrey",type="l", ylim=c(min(beta2)-1, max(beta2)+1),
       xaxt='n', ylab=expression(beta(t)+beta[i](t)), xlab="Time (Seconds)", lty=2, 
       cex.lab = 1, cex.axis=1,cex.main=1, cex.sub=1,
       main=paste0("ROI ",voxel))
  lines(t, beta2[,3], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,4], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,5], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,6], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,7], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,8], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,9], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,10], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,11], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,12], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,13], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,14], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,15], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,16], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,17], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,18], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,19], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,20], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,21], col="darkgrey",type="l", lty=2)
  lines(t, beta2[,1], col="red", type="l", lwd=1.5)
  abline(h=0, col="black", lty=2, lwd=1.2)
  axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)
  dev.off()
  
  
  return(list(RMSE_FMM.train=RMSE_RandomReml, 
              RMSE_FLM.train=RMSE_NoRandomReml,
              RMSE_FMM.test=RMSE.test.fmm,
              RMSE_FLM.test=RMSE.test.flm,
              rsquare=rsquare,
              npc=npc))
}


ROI_4 <- realanalysis_outofsample(4, Whole_fmridata)
ROI_5 <- realanalysis_outofsample(5, Whole_fmridata)
ROI_6 <- realanalysis_outofsample(6, Whole_fmridata)
ROI_16 <- realanalysis_outofsample(16, Whole_fmridata)
ROI_19 <- realanalysis_outofsample(19, Whole_fmridata)
ROI_10 <- realanalysis_outofsample(10, Whole_fmridata)

