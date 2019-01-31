load("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/03.13.2018/real data application/realdata.Rdata")
# 1th voxel 
voxel=10
fmridata <- (center_whole_fmri)[(1+943*(voxel-1)):(943*voxel),]
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
npc <- results$npc # voxel=1 then npc=9
lambda <- results$evalues/D
score <- results$scores
ascore <- score[, 1:npc]/sqrt(D)
efunctions <- results$efunctions*sqrt(D)
par(mfrow=c(4,2), mar=c(2.9,3,1,1), mgp = c(2,1,0), cex.main = 1, cex.axis = 0.8, xaxs = 'r')
myat=seq(0, 1, by=0.2)
for (i in 1:npc){
  plot(t, efunctions[,i],type="l",
       xlab="t",
       ylab="",
       xaxt="n",
       xlim=c(0.0,1.0),
       main=paste("eigenfunction ",i,sep=""),lwd=2, ylim=c(min(efunctions)-0.5,max(efunctions)+0.5))
  axis(1, at = myat, labels = round(myat,digits=1))
}

nSubj=length(unique(fmridata$ID))
nRep = rep(1, nSubj) #nRep in c(39, 48)
for (i in 1:nSubj){
  nRep[i] = length(which(fmridata$ID==i))
}

dummyX <- fmridata$temp
dummyX <- cbind(dummyX, -dummyX + 1) #hot, warm
Y <- fmridata$rating

designMatrix1 <- data.frame(rating = fmridata$rating, 
                            temp.1 = dummyX[,1],
                            temp.2 = dummyX[,2],
                            temp = factor(fmridata$temp, labels=c("warm","hot")),
                            ID = as.factor(fmridata$ID),
                            ascore = ascore)
noRandom.simpart <- paste(1:npc, collapse = " + ascore.")
noRandom.sim1 <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                  noRandom.simpart, 
                                  " + (0 + temp.1 | ID) + (0 + temp.2 | ID) ", 
                                  sep = ""))
NoRandomReml <- lmer(noRandom.sim1, data = designMatrix1)
mseNoRandom <- mean((fmridata$rating - predict(NoRandomReml))^2) #6504.301 
# npc=9
Random.sim1 <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                noRandom.simpart, 
                                " + (0 + temp.1 | ID) + (0 + temp.2 | ID)", 
                                "+ (0 + ascore.1 | ID)",
                                "+ (0 + ascore.2 | ID) + (0 + ascore.3 | ID)",
                                "+ (0 + ascore.4 | ID) + (0 + ascore.5 | ID)", 
                                "+ (0 + ascore.6 | ID) + (0 + ascore.7 | ID)",
                                "+ (0 + ascore.8 | ID)", 
                                "+ (0 + ascore.9 | ID)",
                                sep = ""))
RandomReml <- lmer(Random.sim1, data = designMatrix1)
mseRandom <- mean((fmridata$rating - predict(RandomReml))^2) #5572.025 

# svd decomposition model
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

designMatrix <- data.frame(rating = Y, 
                           temp.1 = dummyX[, 1],
                           temp.2 = dummyX[, 2],
                           ID = as.factor(fmridata$ID),
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
msefullReml<- mean(residuals(fullReml)^2) #5564.165

yhat <- predict(fullReml)
yhat_fullRemlBack <- rep(1, length(yhat)) #transform yhat back

ascore_tran <- score[, 1:npc]/sqrt(D)
for(k in 1:nSubj){
  if(k==1){
    svd <- svd(ascore_tran[1:nRep[1], ] %*% t(ascore_tran[1:nRep[1], ])) #SVD on A_i
  }else{
    svd <- svd(ascore_tran[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ] %*% t(ascore_tran[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ])) #SVD on A_i
  }
  u <- svd$v
  
  if(k==1){
    yhat_fullRemlBack[1:nRep[k]] <- u %*% yhat[1:nRep[k]]
  }else{
    yhat_fullRemlBack[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])] <- u %*% yhat[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k])]
  }
}
mean((yhat_fullRemlBack-fmridata$rating)^2) #5564.165



f.slope <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                            additive0.sim, 
                            " + (0 + z.sim.uni | ID.uni)", 
                            sep = ""))
m.slope <- lmer(f.slope, data = designMatrix.lmm)
f0 <- as.formula(" . ~ . - (0 + z.sim.uni | ID.uni)")
m0 <- update(fullReml, f0)

test <- exactRLRT(m.slope, fullReml, m0)
pvalue <- test$p[1] 


library(ggplot2)
library(reshape2)

fixeffWithTemp <- fixef(fullReml)
randeffWithTemp_IDuni <- as.matrix(ranef(fullReml)$ID.uni[2:(npc*nSubj+1),])
betaWithTemp <- efunctions %*% as.vector(fixeffWithTemp[4:(npc+3)])

theta_ik <- matrix(randeffWithTemp_IDuni,nrow=npc) #defalt is by column
score_thetaik <- score/sqrt(D)
for(k in 1:nSubj){
  if(k==1){
    svd <- svd(score_thetaik[1:nRep[1], ]) #SVD on A_i
  }else{
    svd <- svd(score_thetaik[(sum(nRep[1:(k-1)])+1):sum(nRep[1:k]), ]) #SVD on A_i
  }
  v <- svd$v
  theta_ik[,k] = v %*% theta_ik[,k]
}


betai <- efunctions %*% theta_ik  #23*9 x 9*20 = 23*20
beta <- cbind(betaWithTemp, apply(betai, 2, function(x) x+betaWithTemp))
colnames(beta) <- c("Fixed_effect", paste("Subject", 1:nSubj, sep = ""))

#dfBetaWithTemp <- data.frame(time = t, beta)
#matplot(t, beta[,2:21], col=colors()[c(3, 9, 26, 31, 568, 48, 509, 70,
#                                411, 89, 340, 117, 142, 151, 503,
#                                125, 468, 652, 542, 544)],
#        type="l", lty=2, ylab = expression(beta(t)), xaxt='n')
#lines(t, beta[,1], type="l", lwd=1.5, col="red")
#abline(h=0, lty=2, lwd=1.2)

par(mgp=c(2,1,0), mar=c(3, 3, 1, 1) + 0.1, mfrow=c(1,1))
plot(t, beta[,2], col=colors()[c(3)],type="l", ylim=c(min(beta)-1, max(beta)+1),
     xaxt='n', ylab=expression(beta(t)+beta[i](t)), xlab="Time (Seconds)", lty=2, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="ROI 10")
lines(t, beta[,3], col=colors()[c(9)],type="l", lty=2)
lines(t, beta[,4], col=colors()[c(26)],type="l", lty=2)
lines(t, beta[,5], col=colors()[c(31)],type="l", lty=2)
lines(t, beta[,6], col=colors()[c(568)],type="l", lty=2)
lines(t, beta[,7], col=colors()[c(48)],type="l", lty=2)
lines(t, beta[,8], col=colors()[c(509)],type="l", lty=2)
lines(t, beta[,9], col=colors()[c(70)],type="l", lty=2)
lines(t, beta[,10], col=colors()[c(411)],type="l", lty=2)
lines(t, beta[,11], col=colors()[c(89)],type="l", lty=2)
lines(t, beta[,12], col=colors()[c(340)],type="l", lty=2)
lines(t, beta[,13], col=colors()[c(117)],type="l", lty=2)
lines(t, beta[,14], col=colors()[c(142)],type="l", lty=2)
lines(t, beta[,15], col=colors()[c(151)],type="l", lty=2)
lines(t, beta[,16], col=colors()[c(503)],type="l", lty=2)
lines(t, beta[,17], col=colors()[c(125)],type="l", lty=2)
lines(t, beta[,18], col=colors()[c(468)],type="l", lty=2)
lines(t, beta[,19], col=colors()[c(652)],type="l", lty=2)
lines(t, beta[,20], col=colors()[c(542)],type="l", lty=2)
lines(t, beta[,21], col=colors()[c(544)],type="l", lty=2)
lines(t, beta[,1], col="red", type="l", lwd=1.5)
abline(h=0, col="black", lty=2, lwd=1.2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)
legend(x = "bottom",inset = 0,
       legend = c(paste("Subject", 1:nSubj, sep = ""), "Population effect"), 
       col=colors()[c(3, 9, 26, 31, 568, 48, 509, 70,
                      411, 89, 340, 117, 142, 151, 503,
                      125, 468, 652, 542, 544, 552)], 
       lwd=c(rep(1,20),1.5),lty=c(rep(2,20),1), 
       cex=.5, ncol = 3, bty='n')

plot(t, betai[,1], col=colors()[c(3)],type="l", ylim=c(min(betai)-1, max(betai)+1),
     xaxt='n', ylab=expression(beta[i](t)), xlab="Time (Seconds)", lty=2, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="ROI 10")
lines(t, betai[,2], col=colors()[c(9)],type="l", lty=2)
lines(t, betai[,3], col=colors()[c(26)],type="l", lty=2)
lines(t, betai[,4], col=colors()[c(31)],type="l", lty=2)
lines(t, betai[,5], col=colors()[c(568)],type="l", lty=2)
lines(t, betai[,6], col=colors()[c(48)],type="l", lty=2)
lines(t, betai[,7], col=colors()[c(509)],type="l", lty=2)
lines(t, betai[,8], col=colors()[c(70)],type="l", lty=2)
lines(t, betai[,9], col=colors()[c(411)],type="l", lty=2)
lines(t, betai[,10], col=colors()[c(89)],type="l", lty=2)
lines(t, betai[,11], col=colors()[c(340)],type="l", lty=2)
lines(t, betai[,12], col=colors()[c(117)],type="l", lty=2)
lines(t, betai[,13], col=colors()[c(142)],type="l", lty=2)
lines(t, betai[,14], col=colors()[c(151)],type="l", lty=2)
lines(t, betai[,15], col=colors()[c(503)],type="l", lty=2)
lines(t, betai[,16], col=colors()[c(125)],type="l", lty=2)
lines(t, betai[,17], col=colors()[c(468)],type="l", lty=2)
lines(t, betai[,18], col=colors()[c(652)],type="l", lty=2)
lines(t, betai[,19], col=colors()[c(542)],type="l", lty=2)
lines(t, betai[,20], col=colors()[c(544)],type="l", lty=2)
abline(h=0, col="black", lty=2, lwd=1.2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)
legend(x = "bottom",inset = 0,
       legend = paste("Subject", 1:nSubj, sep = ""), 
       col=colors()[c(3, 9, 26, 31, 568, 48, 509, 70,
                      411, 89, 340, 117, 142, 151, 503,
                      125, 468, 652, 542, 544)], 
       lwd=rep(1,20),lty=rep(2,20), 
       cex=.45, ncol = 3, bty='n')

# use ggplot2
#par(mfrow=c(1,1))
#dfBetaWithTemp <- melt(dfBetaWithTemp, id = c("time"))
#ggplot(data = dfBetaWithTemp, aes(x=time, y=value)) +
#  geom_line(data=dfBetaWithTemp[dfBetaWithTemp$variable!="Fixed_effect",],
#            aes(colour=variable, group = variable),
#            size = .5,
#            alpha = 1) +
#  geom_line(data=dfBetaWithTemp[dfBetaWithTemp$variable=="Fixed_effect",],
#            aes(colour=variable, group = variable),
#            size = 0.8,
#            alpha = 1) +
#  ggtitle(expression(beta~'(t)')) 


#dfbetai<- data.frame(time = t, betai)#23*20
#colnames(dfbetai) <- c("time", paste("Subject", 1:nSubj, sep = ""))
#dfbetai <- melt(dfbetai, id = c("time"))

#ggplot(data = dfbetai, aes(x=time, y=value)) +
#  geom_point(aes(color = variable)) +
#  geom_line(data=dfbetai[dfbetai$variable!="Fixed_effect",],
#            aes(colour=variable, group = variable),
#            size = .5,
#            alpha = .5) +
#  ggtitle(expression(beta[i]~'(t)'))

### pick 4 subjects, plot the beta_i(t) and compare ###
### pick 4, 8, 12, 17 ###
### NoRandomReml
fixeffWithTemp1 <- fixef(NoRandomReml)
betaWithTemp1 <- efunctions %*% as.vector(fixeffWithTemp1[3:(npc+2)])

### RandomReml
fixeffWithTemp2 <- fixef(RandomReml)
randeffWithTemp_IDuni2 <- as.matrix(ranef(RandomReml)$ID[,3:(npc+2)])
betaWithTemp2 <- efunctions %*% as.vector(fixeffWithTemp2[3:(npc+2)])

theta_ik2 <- t(randeffWithTemp_IDuni2) #defalt is by column
betai2 <- efunctions %*% theta_ik2  #23*8 x 8*20 = 23*20
beta2 <- cbind(betaWithTemp2, apply(betai2, 2, function(x) x+betaWithTemp2))
colnames(beta2) <- c("Fixed_effect", paste("Subject", 1:nSubj, sep = ""))

######  violin plot ######
resi2_fullReml <- (yhat_fullRemlBack-fmridata$rating)^2
resi2_RandomReml <- residuals(RandomReml)^2
resi2_NoRandomReml <- residuals(NoRandomReml)^2
violin_fullReml <- violin_RandomReml <- violin_NoRandomReml <- c()
for (i in 1:length(nRep)){
  if(i==1){
    violin_fullReml <- c(violin_fullReml, mean(resi2_fullReml[1:nRep[i]]))
    violin_RandomReml <- c(violin_RandomReml, mean(resi2_RandomReml[1:nRep[i]]))
    violin_NoRandomReml <- c(violin_NoRandomReml, mean(resi2_NoRandomReml[1:nRep[i]]))
  }else{
    violin_fullReml <- c(violin_fullReml, mean(resi2_fullReml[(sum(nRep[1:(i-1)])+1):sum(nRep[1:i])]))
    violin_RandomReml <- c(violin_RandomReml, mean(resi2_RandomReml[(sum(nRep[1:(i-1)])+1):sum(nRep[1:i])]))
    violin_NoRandomReml <- c(violin_NoRandomReml, mean(resi2_NoRandomReml[(sum(nRep[1:(i-1)])+1):sum(nRep[1:i])]))
  }
}

violin1 <- data.frame(cbind(c(violin_NoRandomReml, violin_fullReml, violin_RandomReml),
                            rep(1:3, each=nSubj)))
colnames(violin1) <- c("MSE","Model")
violin1$Model <- as.factor(violin1$Model)

library(ggplot2)
ggplot(violin1, aes(x=Model, y=MSE)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.alpha = 0.1) +
  geom_jitter( width = 0.2, alpha = 0.6) +
  geom_hline(aes(yintercept = median(violin_RandomReml)), color="red", linetype="dashed", size = 0.8)+
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("No ranef","Homogenous", "Heteroscedastic")) +
  labs(title="ROI 10") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = rel(0.9), face="bold")) +
  theme(axis.title.y = element_text(size = rel(0.8))) +
  theme(axis.title.x = element_text(size = rel(0.8))) 

resi_fullReml <- yhat_fullRemlBack-fmridata$rating
resi_RandomReml <- residuals(RandomReml)
resi_NoRandomReml <- residuals(NoRandomReml)
violin1_resi <- data.frame(cbind(c(resi_NoRandomReml, resi_fullReml, resi_RandomReml),
                                 rep(1:3, each=dim(fmridata)[1])))
colnames(violin1_resi) <- c("Residual","Model")
violin1_resi$Model <- as.factor(violin1_resi$Model)

ggplot(violin1_resi, aes(x=Model, y=Residual)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.alpha = 0.1) +
  geom_jitter( width = 0.3, alpha = 0.3) +
  #geom_hline(aes(yintercept = median(resi_RandomReml)), color="red", linetype="dashed", size = 0.8)+
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("No ranef","Homogenous", "Heteroscedastic")) +
  labs(title="ROI 10") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = rel(0.9), face="bold")) +
  theme(axis.title.y = element_text(size = rel(0.8))) +
  theme(axis.title.x = element_text(size = rel(0.8)))

# compare suject 4 with beta(t)

f_ymin = min(c(beta2[,c(2,5,10,20)], beta[,c(2,5,10,20)], betaWithTemp1))-1
f_ymax = max(c(beta2[,c(2,5,10,20)], beta[,c(2,5,10,20)], betaWithTemp1))+1
fi_ymin = min(c(betai2[,c(1,4,9,19)], betai[,c(1,4,9,19)]))-1
fi_ymax = max(c(betai2[,c(1,4,9,19)], betai[,c(1,4,9,19)]))+1
par(mfrow=c(2,2), mar=c(3, 3, 1, 1)+0.1)

# sub 1, beta+betai
plot(t, beta2[,2],type="l", ylim=c(f_ymin,f_ymax),
     xaxt='n', ylab=expression(beta(t)+beta[i](t)), xlab="Time (Seconds)", lty=1, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="Suject 1 (ROI 10)", col="blue", lwd=2)
lines(t, beta[,2],type="l", col="green", lwd=2, lty=3)
#lines(t, betaWithTemp1, type="l",col="red", lty=2, lwd=2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)

# sub 4, beta+betai
plot(t, beta2[,5],type="l", ylim=c(f_ymin,f_ymax),
     xaxt='n', ylab=expression(beta(t)+beta[i](t)), xlab="Time (Seconds)", lty=1, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="Suject 4 (ROI 10)", col="blue", lwd=2)
lines(t, beta[,5],type="l", col="green", lwd=2, lty=3)
#lines(t, betaWithTemp1, type="l",col="red", lty=2, lwd=2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)

# sub 9, beta+betai
plot(t, beta2[,10],type="l", ylim=c(f_ymin,f_ymax),
     xaxt='n', ylab=expression(beta(t)+beta[i](t)), xlab="Time (Seconds)", lty=1, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="Suject 9 (ROI 10)", col="blue", lwd=2)
lines(t, beta[,10],type="l", col="green", lwd=2, lty=3)
#lines(t, betaWithTemp1, type="l",col="red", lty=2, lwd=2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)

# sub 19, beta+betai
plot(t, beta2[,20],type="l", ylim=c(f_ymin,f_ymax),
     xaxt='n', ylab=expression(beta(t)+beta[i](t)), xlab="Time (Seconds)", lty=1, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="Suject 19 (ROI 10)", col="blue", lwd=2)
lines(t, beta[,20],type="l", col="green", lwd=2, lty=3)
#lines(t, betaWithTemp1, type="l",col="red", lty=2, lwd=2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)

# sub 1, betai
plot(t, betai2[,1],type="l", ylim=c(fi_ymin,fi_ymax),
     xaxt='n', ylab=expression(beta[i](t)), xlab="Time (Seconds)", lty=1, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="Suject 1 (ROI 10)", col="blue", lwd=2)
lines(t, betai[,1],type="l", col="green", lwd=2, lty=3)
abline(h=0, col="red", lwd=1.2, lty=2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)

# sub 4, betai
plot(t, betai2[,4],type="l", ylim=c(fi_ymin,fi_ymax),
     xaxt='n', ylab=expression(beta[i](t)), xlab="Time (Seconds)", lty=1, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="Suject 4 (ROI 10)", col="blue", lwd=2)
lines(t, betai[,4],type="l", col="green", lwd=2, lty=3)
abline(h=0, col="red", lwd=1.2, lty=2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)
# sub 9, betai
plot(t, betai2[,9],type="l", ylim=c(fi_ymin,fi_ymax),
     xaxt='n', ylab=expression(beta[i](t)), xlab="Time (Seconds)", lty=1, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="Suject 9 (ROI 10)", col="blue", lwd=2)
lines(t, betai[,9],type="l", col="green", lwd=2, lty=3)
abline(h=0, col="red", lwd=1.2, lty=2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)
# sub 19, betai
plot(t, betai2[,19],type="l", ylim=c(fi_ymin,fi_ymax),
     xaxt='n', ylab=expression(beta[i](t)), xlab="Time (Seconds)", lty=1, 
     cex.lab = 0.8, cex.axis=0.8,cex.main=0.8, cex.sub=0.8,
     main="Suject 19 (ROI 10)", col="blue", lwd=2)
lines(t, betai[,19],type="l", col="green", lwd=2, lty=3)
abline(h=0, col="red", lwd=1.2, lty=2)
axis(1, at=t[c(5,10,15,20)],labels=c(10,20,30,40), cex.lab=0.8, cex.axis=0.8)
