### RandomReml
designMatrix1 <- data.frame(rating = fmridata$rating, 
                            temp.1 = dummyX[,1],
                            temp.2 = dummyX[,2],
                            temp = factor(fmridata$temp, labels=c("warm","hot")),
                            ID = as.factor(fmridata$ID),
                            ascore = ascore)
### fitting process
# npc=8
noRandom.simpart <- paste(1:npc, collapse = " + ascore.")
Random.sim1 <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                noRandom.simpart, 
                                " + (0 + temp.1 | ID) + (0 + temp.2 | ID)", 
                                "+ (0 + ascore.1 | ID)",
                                "+ (0 + ascore.2 | ID) + (0 + ascore.3 | ID)",
                                "+ (0 + ascore.4 | ID) + (0 + ascore.5 | ID)", 
                                "+ (0 + ascore.6 | ID) + (0 + ascore.7 | ID)",
                                "+ (0 + ascore.8 | ID)", 
                                sep = ""))
RandomReml <- lmer(Random.sim1, data = designMatrix1)
mseRandom <- mean((fmridata$rating - predict(RandomReml))^2) #5428.054 

#estimation 
fixeffWithTemp2 <- fixef(RandomReml)
randeffWithTemp_IDuni2 <- as.matrix(ranef(RandomReml)$ID[,3:(npc+2)]) #first two column is for temp.1 and temp.2
betaWithTemp2 <- efunctions %*% as.vector(fixeffWithTemp2[3:(npc+2)])

theta_ik2 <- t(randeffWithTemp_IDuni2) #defalt is by column
betai2 <- efunctions %*% theta_ik2  #23*8 x 8*20 = 23*20, prediction for betai(t), (beta1(t),...beta20(t))
beta2 <- cbind(betaWithTemp2, apply(betai2, 2, function(x) x+betaWithTemp2)) #(beta(t), beta(t)+beta1(t),...,beta(t)+beta20(t)))
colnames(beta2) <- c("Fixed_effect", paste("Subject", 1:nSubj, sep = ""))
