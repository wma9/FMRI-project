library(parallel)
simRep <- 1000 # Replication times in one simulation
#pvalue.true <- .05 # Testing type I error 
#b.var <- c(0.04) # The set of varaince of random covariates b as random slope
#smooth <- 1 # measurement error is added to M if smooth = 0; no measurement error is added if sooth = 1
cores <- 8
#r.sim <- b.var
nSubj <- 20
nRep <- 50
seediter <- 1  

run_one_sample <- function(iter){
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  library(MASS)
  
  #set.seed(iter)
  nNewvisit <- 10
  set.seed(iter+(seediter-1)*5000)
  D <- 80 # grid number total
  #nSubj <- 20 # 200 # I the number of curves
  #nRep <-  20 # 20 # datasets for each covariance function 
  totalN <- nSubj * nRep
  thetaK.true <- 2
  timeGrid <- (1:D)/D
  
  npc.true <- 3
  percent <- 0.95
  SNR <- 3 # 5, signal noise ratio'
  
  sd.epsilon <- 1 # or 0.5
  delta.true <- 0.5
  a.mean <- 0
  
  gamma.true <- 2
  gammaVar.true <- 1
  # hot
  gammaI.true.i <- mapply(rnorm, nSubj, 0, rep(sqrt(gammaVar.true), 1))
  gammaI.true <- gammaI.true.i[rep(1:nrow(gammaI.true.i), each = nRep), ]
  # warm
  #gammaI2.true.i <- mapply(rnorm, nSubj, gamma.true, rep(sqrt(gammaVar.true), 1))
  #gammaI2.true <- gammaI2.true.i[rep(1:nrow(gammaI2.true.i), each = nRep), ]
  
  dummyX <- rbinom(n = totalN, size = 1, prob = 0.5) # dummyX
  
  #generate functional covariates
  lambda.sim <- function(degree) {
    return(0.5^(degree - 1))
  }
  psi.fourier <- function(t, degree) {
    result <- NA
    if(degree == 1){
      result <- sqrt(2) * sinpi(2*t)
    }else if(degree == 2){
      result <- sqrt(2) * cospi(4*t)
    }else if(degree == 3){
      result <- sqrt(2) * sinpi(4*t)
    }
    return(result)
  }
  
  lambdaVec.true <- mapply(lambda.sim, 1: npc.true) 
  
  psi.true <- matrix(data = mapply(psi.fourier, rep(timeGrid, npc.true), rep(1:npc.true, each=D)),
                     nrow = npc.true, 
                     ncol = D, 
                     byrow = TRUE)
  ascore.true <- mvrnorm(totalN, rep(a.mean, npc.true), diag(lambdaVec.true))   
  Mt.true <-  ascore.true %*% psi.true     
  error <- rnorm(totalN, mean = 0, sd = sd.epsilon)
  
  #thetaIK.true <- mvrnorm(nSubj, rep(thetaK.true, npc.true), diag(rep(r.sim, npc.true)))
  thetaIK.true1 <- mvrnorm(nSubj, rep(0, npc.true), diag(c(r.sim, r.sim/2, r.sim/4)))
  thetaIK.true <- thetaIK.true1[rep(1:nrow(thetaIK.true1), each = nRep), ]
  betaM.true <- (thetaK.true+thetaIK.true) * ascore.true
  betaM.true <- rowSums(betaM.true)
  
  
  #Y <- delta.true + dummyX * gammaI.true + (dummyX - 1) * gammaI2.true + betaM.true + error
  Y <- delta.true + dummyX * (gammaI.true+gamma.true) + betaM.true + error
  
  
  ##########################################################################
  
  ID <- rep(1:nSubj, each = nRep)
  if(smooth == 0){
    Merror.Var <- sum(lambdaVec.true) / SNR  #SNR = sum(lambdaVec.true)/Merror.Var
    Mt.hat <- Mt.true + matrix(rnorm(totalN*D, mean=0, sd = sqrt(Merror.Var)), totalN, D)
  }
  if(smooth == 1){
    Merror.Var <- 0  #SNR = sum(lambdaVec.true)/Merror.Var
    Mt.hat <- Mt.true 
  }
  
  M <- Mt.hat
  # M <- M - matrix(rep(colMeans(M), each = totalN), totalN, D) # center:column-means are 0
  
  t <- (1:D)/D
  knots <- 5 #5 previous setting 10
  p <- 5  # previous setting p <- 7, the number of degree for B-splines we use
  
  results <- fpca.face(M, center = TRUE, argvals = t, knots = knots, pve = percent, p = p, lambda = 0) # pve need to be chosen!
  npc <- results$npc
  score <- results$scores
  ascore1 <- score[, 1:npc]/sqrt(D)
  efunctions1 <- results$efunctions*sqrt(D)
  sign_correct <- diag(sign(colSums(efunctions1*t(psi.true))))
  efunctions <- efunctions1 %*% sign_correct
  ascore <- ascore1 %*% sign_correct
##############################################################
#dummyX <- cbind(dummyX, -dummyX + 1)

##############################################################
designMatrix <- data.frame(rating = Y, 
                           temp.1 = dummyX,
                           #temp.2 = dummyX[, 2],
                           ID = as.factor(ID),
                           ascore = ascore)
  
noRandom.simpart <- paste(1:npc, collapse = " + ascore.")
noRandom.sim1 <- as.formula(paste("rating ~ 1 + temp.1  + ascore.", 
                                  noRandom.simpart, 
                                  " + (0 + temp.1 | ID) ", 
                                  sep = ""))

NoRandomReml <- lmer(noRandom.sim1, data = designMatrix)

# npc=3
additive.heter <- paste0(" + (0 + ascore.", 1:npc, " | ID)", collapse = "")

Random.sim1 <- as.formula(paste("rating ~ 1 + temp.1  + ascore.", 
                                noRandom.simpart, 
                                " + (0 + temp.1 | ID) ", 
                                additive.heter, 
                                sep = ""))
RandomReml <- lmer(Random.sim1, data = designMatrix)

###true beta_t
beta_t <- t(psi.true)%*%as.vector(rep(thetaK.true, npc.true))
beta_it <- t(psi.true) %*% t(thetaIK.true1)

### NoRandomReml
fixeffWithTemp1 <- fixef(NoRandomReml)
betaWithTemp1 <- efunctions %*% as.vector(fixeffWithTemp1[3:(npc+2)])
ISE_beta_1 <- sum((betaWithTemp1-beta_t)^2/D)      #ISE for beta(t)
e1 <- Y - predict(NoRandomReml)
mseNoRandom <- mean(e1^2) #mse for response
NoRandom_vcov <- summary(NoRandomReml)$vcov

### RandomReml
fixeffWithTemp2 <- fixef(RandomReml)
betaWithTemp2 <- efunctions %*% as.vector(fixeffWithTemp2[3:(npc+2)])
ISE_beta_2 <- sum((betaWithTemp2-beta_t)^2/D)  #ISE for beta(t)
e2 <- Y - predict(RandomReml)
mseRandom <- mean(e2^2) #mse for response
betai_2 <- efunctions %*% t(as.matrix(ranef(RandomReml)$ID[,2:(npc+1)])) # beta_i(t)
MISE_2 <- mean(colMeans((betai_2-beta_it)^2)) #MISE for beta_i(t)
Random_vcov <- summary(RandomReml)$vcov

gene_newdata <- function(newvisitN){
  #generate new dummyX
  dummyX <- rbinom(n = nSubj*newvisitN, size = 1, prob = 0.5) # new dummyX
  #generate new x_ij(t) for new visit
  ascore.true <- mvrnorm(nSubj*newvisitN, rep(a.mean, npc.true), diag(lambdaVec.true)) 
  #sum new x_ij(t)(beta_t + betai_t)
  Mt.true <- ascore.true %*% psi.true    
  #generate new error
  error <- rnorm(nSubj*newvisitN, mean = 0, sd = sd.epsilon)
  
  thetaIK.true <- thetaIK.true1[rep(1:nrow(thetaIK.true1), each = newvisitN), ]
  betaM.true <- (thetaIK.true+thetaK.true) * ascore.true
  betaM.true <- rowSums(betaM.true)
  
  # hot
  gammaI.true <- gammaI.true.i[rep(1:nrow(gammaI.true.i), each = newvisitN), ]
  # warm
  #gammaI2.true <- gammaI2.true.i[rep(1:nrow(gammaI2.true.i), each = newvisitN), ]
  # generate new Y
  Y <- delta.true + dummyX * (gammaI.true+gamma.true) + betaM.true + error
  totalN <- nSubj*newvisitN
  ID <- rep(1:nSubj, each = newvisitN)
  
  if(smooth == 0){
    Merror.Var <- sum(lambdaVec.true) / SNR  #SNR = sum(lambdaVec.true)/Merror.Var
    Mt.hat <- Mt.true + matrix(rnorm(totalN*D, mean=0, sd = sqrt(Merror.Var)), totalN, D)
  }
  if(smooth == 1){
    Merror.Var <- 0  #SNR = sum(lambdaVec.true)/Merror.Var
    Mt.hat <- Mt.true 
  }
  
  M <- Mt.hat
  #projection to get new xi_{ijk}
  xi <- M %*% efunctions/D
  
  designMatrix <- data.frame(rating = Y, 
                             temp.1 = dummyX,
                             #temp.2 = -dummyX + 1,
                             ID = as.factor(ID),
                             ascore = xi)
  return(list(newY=Y,
              newdata=designMatrix,
              Xt = M))
}

#predict for the new dataset
test1 <- gene_newdata(nNewvisit)
#noRandom beta_t
#betaIT1.i <- betaWithTemp1
#betaIT1 <- betaIT1.i[,rep(1,nSubj*nNewvisit)]
#Random beta_t
#betaIT2.i <- betaWithTemp2 %*% rep(1,nSubj) + betai_2
#betaIT2 <- betaIT2.i[, rep(1:nSubj, each=nNewvisit)]
#calculate new MSE
#de<-model.matrix(~0+test1$newdata$ID)
#test_est1 <- rowMeans(test1$Xt * t(betaIT1)) + as.vector(test1$newdata$temp.1)*(de%*% as.vector(ranef(NoRandomReml)$ID$temp.1)) + as.matrix(cbind(1,as.matrix(test1$newdata$temp.1))) %*% as.matrix(fixeffWithTemp1[1:2])
#test_est2 <- rowMeans(test1$Xt * t(betaIT2)) + as.vector(test1$newdata$temp.1)*(de%*% as.vector(ranef(RandomReml)$ID$temp.1)) + cbind(1,as.matrix(test1$newdata$temp.1)) %*% as.matrix(fixeffWithTemp2[1:2])
#test_est_MSE1 <- mean((test_est1 - test1$newY)^2)
#test_est_MSE2 <- mean((test_est2 - test1$newY)^2)

test_est1 <- predict(NoRandomReml, test1$newdata)
test_est2 <- predict(RandomReml, test1$newdata)
test_est_MSE1 <- mean((test_est1 - test1$newY)^2)
test_est_MSE2 <- mean((test_est2 - test1$newY)^2)

#test_est1 <- as.matrix(test1$newdata[,4:(3+npc)])%*%as.matrix(fixeffWithTemp1[3:(2+npc)]) + as.vector(test1$newdata$temp.1)*(de%*% as.vector(ranef(NoRandomReml)$ID$temp.1)) + as.matrix(cbind(1,as.matrix(test1$newdata$temp.1))) %*% as.matrix(fixeffWithTemp1[1:2])
#test_est2 <- as.matrix(test1$newdata[,4:(3+npc)])%*%as.matrix(fixeffWithTemp2[3:(2+npc)]) + rowSums(as.matrix(test1$newdata[,4:(3+npc)]) * (de%*% as.matrix(ranef(RandomReml)$ID[,2:(1+npc)]))) + as.vector(test1$newdata$temp.1)*(de%*% as.vector(ranef(RandomReml)$ID$temp.1)) + cbind(1,as.matrix(test1$newdata$temp.1)) %*% as.matrix(fixeffWithTemp2[1:2])
#test_est_MSE1 <- mean((test_est1 - test1$newY)^2)
#test_est_MSE2 <- mean((test_est2 - test1$newY)^2)

return(list(realTau = r.sim,  
            ISE_beta.norandom = ISE_beta_1,
            ISE_beta.random = ISE_beta_2,
            mse_Y.NoRandom = mseNoRandom,
            mse_Y.Random = mseRandom,
            MISE_betai.Random = MISE_2,
            beta.norandom = betaWithTemp1,
            beta.random = betaWithTemp2,
            betait.random = betai_2,
            NoRandom_vcov = NoRandom_vcov,
            Random_vcov = Random_vcov,
            test_est_MSE.noran = test_est_MSE1,
            test_est_MSE.ran = test_est_MSE2))
}

# Setup parallel
#cores <- detectCores()
cluster <- makeCluster(cores)
clusterExport(cluster, c("nSubj", "nRep", "seediter"))

for(r.sim in c(.02, .04, .08)){
  for (smooth in c(1,0)){
    clusterExport(cluster, c("r.sim", "smooth")) # casting the coefficient parameter on the random effects' covariance function
    b.var <- r.sim
    fileName <- paste("heter_pred_", smooth, "_",b.var,"_seed",
                      seediter,"_grp", nSubj, "-rep", nRep, 
                      ".RData", sep = "") # Saving file's name
    
    # run the simulation
    loopIndex <- 1
    est.sim <- list()
    
    node_results <- parLapply(cluster, 1:simRep, run_one_sample)
    
    # result1.sim <- lapply(node_results, function(x) {list(realTau = x$realTau, 
    #                                                      pvalue = x$pvalue)})
    #result2.sim <- lapply(node_results, function(x) {list(realTau = x$realTau, 
    #                                                      pvalues.bonf = x$pvalues.bonf,
    #												      smooth = x$smooth,
    #												      npc = x$npc)})
    #resultDoubleList.sim[[loopIndex]] <- node_results
    
    #save.image(file=fileName) # Auto Save
    
    #table1.sim <- sapply(result1.sim, function(x) {       
    #  c(sens = (sum(x$pvalue <= pvalue.true) > 0))})
    #Power1 <- mean(table1.sim)
    #cat("nRandCovariate: ", nRandCovariate, fill = TRUE)
    #cat("Power1: ", Power1, fill = TRUE)
    #power1.sim[[loopIndex]] <- list(Power = Power1, realTau = r.sim) 
    
    
    #ISE_beta.norandom,ISE_beta.random, mse_Y.NoRandom, mse_Y.Random, MISE_betai.Random
    est <- sapply(node_results, function(x){       
      return(c(x$ISE_beta.norandom, x$ISE_beta.random, x$mse_Y.NoRandom, x$mse_Y.Random, x$MISE_betai.Random,
               x$test_est_MSE.noran, x$test_est_MSE.ran))})
    est.mean <- rowMeans(est)
    
    beta.norandom <- sapply(node_results, function(x){       
      x$beta.norandom})
    beta.random <- sapply(node_results, function(x){       
      x$beta.random})
    betait.random <- lapply(node_results, function(x){       
      x$betait.random})
    NoRandom_vcov <- lapply(node_results, function(x){       
      x$NoRandom_vcov})
    Random_vcov <- lapply(node_results, function(x){       
      x$Random_vcov})
    
    
    est.sim <- list(est.mean = list(test_est_MSE.noran = est.mean[6],
                                    test_est_MSE.ran = est.mean[7],
                                    ISE_beta.norandom = est.mean[1],
                                    ISE_beta.random = est.mean[2],
                                    mse_Y.NoRandom = est.mean[3], 
                                    mse_Y.Random = est.mean[4], 
                                    MISE_betai.Random = est.mean[5]), 
                                realTau = c(r.sim,r.sim/2,r.sim/4), 
                                smooth = smooth)
    
    save(est.sim, beta.norandom,beta.random,betait.random, NoRandom_vcov,Random_vcov, file=fileName) # Auto Save
  }
}

stopCluster(cluster)

