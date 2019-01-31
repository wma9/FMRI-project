# data will be generated from tau1=tau^2, tau2=tau^2/2, and tau3=tau^2/4
library(parallel)
simRep <- 5000 # Replication times in one simulation
pvalue.true <- .05 # Testing type I error 
#b.var <- c(0.04) # The set of varaince of random covariates b as random slope
smooth <- 0 # measurement error is added to M if smooth = 0; no measurement error is added if sooth = 1
cores <- 4
#r.sim <- b.var


run_one_sample <- function(iter){
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  library(MASS)
  
  set.seed(iter+10000)
  D <- 80 # grid number total
  nSubj <- 20 # 200 # I the number of curves
  nRep <-  20 # 20 # datasets for each covariance function 
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
  gammaI.true <- mapply(rnorm, nSubj, gamma.true, rep(sqrt(gammaVar.true), 1))
  gammaI.true <- gammaI.true[rep(1:nrow(gammaI.true), each = nRep), ]
  # warm
  gammaI2.true <- mapply(rnorm, nSubj, gamma.true, rep(sqrt(gammaVar.true), 1))
  gammaI2.true <- gammaI2.true[rep(1:nrow(gammaI2.true), each = nRep), ]
  
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
  thetaIK.true <- mvrnorm(nSubj, rep(thetaK.true, npc.true), diag(c(r.sim, r.sim/2, r.sim/4)))
  thetaIK.true <- thetaIK.true[rep(1:nrow(thetaIK.true), each = nRep), ]
  betaM.true <- thetaIK.true * ascore.true
  betaM.true <- rowSums(betaM.true)
  
  
  Y <- delta.true + dummyX * gammaI.true + (dummyX - 1) * gammaI2.true + betaM.true + error
  
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
  knots <- 5 # previous setting 10
  p <- 5  # previous setting p <- 7, the number of degree for B-splines we use
  
  results <- fpca.face(M, center = TRUE, argvals = t, knots = knots, pve = percent, p = p, lambda = 0) # pve need to be chosen!
  npc <- results$npc
  score <- results$scores
  ascore <- score[, 1:npc]/sqrt(D)
  
  # plot(results$efunctions[,2]*sqrt(D))
  # lines(1:80, psi.fourier(timeGrid, 2)) #match very well 
  # to compare lambda: results$evalues/(D)
  # to compare estimated M, Mt.hat, Mt.true
  # a<-results$scores %*% t(results$efunctions)
  # plot(M[300,]) #Mt.hat
  # lines(a[300,]+results$mu,col="red") # estimated M
  # lines(Mt.true[300,], col="blue") #true Mt
  
  ###########################################################################
  dummyX <- cbind(dummyX, -dummyX + 1)
  
  ###########################################################################
  designMatrix <- data.frame(rating = Y, 
                             temp.1 = dummyX[, 1],
                             temp.2 = dummyX[, 2],
                             ID = as.factor(ID),
                             ascore = ascore)
  
  
  # 'lmer' model
  designMatrix.lmm <- designMatrix
  
  additive0.sim <- paste(1:npc, collapse = " + ascore.")
  additive.sim <- paste(1:npc, collapse = " | ID) + (0 + ascore.")
  additive.heter <- paste0(" + (0 + ascore.", 1:npc, " | ID)", collapse = "")
  
  # Confusion of modifying
  #model.sim <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
  #                              additive0.sim, 
  #                              " + (0 + temp.1 | ID) + (0 + temp.2 | ID) + (0 + z.sim.uni | ID.uni)", 
  #                              sep = ""))
  
  model.sim <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                additive0.sim, 
                                " + (0 + temp.1 | ID) + (0 + temp.2 | ID) ", 
                                additive.heter,
                                sep = ""))
  test_time <-system.time({
    
    tests2 <- list()
    for(i in 1:npc){
      ii <- paste("ascore.", i, sep = "")
      #model only contains the random effect that want to be tested
      f.slope <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                  additive0.sim, " + (0 +", ii, " | ID)", 
                                  sep = ""))
      m.slope <- lmer(f.slope, data = designMatrix.lmm)
      #full model under the alternative
      mA <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                             additive0.sim, " + (0 +", ii, " | ID)", 
                             "+ (0 + temp.1 | ID) + (0 + temp.2 | ID)",
                             sep = ""))
      fullReml <- lmer(mA, data = designMatrix.lmm)
      #model under the null
      f0 <- as.formula(paste(" . ~ . - (0 + ", ii, "| ID)"))
      m0 <- update(fullReml, f0)
      tests2[[i]] <- exactRLRT(m.slope, fullReml, m0)
    }
    multiTest1 <- sapply(tests2, function(x) {
      c(statistic = x$statistic[1],
        "p-value" = x$p[1])})
    pvalues.bonf <- p.adjust(multiTest1[2,], "bonferroni")
  })  
  
  #f.slope <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
  #                            additive0.sim, 
  #                            " + (0 + z.sim.uni | ID.uni)", 
  #                            sep = ""))
  #m.slope <- lmer(f.slope, data = designMatrix.lmm)
  #f0 <- as.formula(" . ~ . - (0 + z.sim.uni | ID.uni)")
  #m0 <- update(fullReml, f0)
  
  #tests2 <- exactRLRT(m.slope, fullReml, m0)
  #pvalues.bonf <- tests2$p[1] 
  
  ###################################################################################
  
  
  return(list(realTau = r.sim,  
              pvalues.bonf = pvalues.bonf, 
              Merror.Var = Merror.Var,
              smooth = smooth,
              npc = npc, 
              tests2 = tests2,
              test_time = test_time[3]))
}

# Setup parallel
#cores <- detectCores()
cluster <- makeCluster(cores)
for(r.sim in c(.1^2, .12^2, .15^2, .2^2, .25^2, .29^2, .31^2, .35^2, .4^2)){
  clusterExport(cluster, c("r.sim", "smooth")) # casting the coefficient parameter on the random effects' covariance function
  b.var <- r.sim
  fileName <- paste("test_power_", smooth, "_",b.var,"_seed3_grp20-rep20.RData", sep = "") # Saving file's name
  
  # run the simulation
  loopIndex <- 1
  power2.sim <- list()
  
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
  
  table2.sim <- sapply(node_results, function(x){       
    c(overall.sens = (sum(x$pvalues.bonf <= pvalue.true) > 0))})
  Power2 <- mean(table2.sim)
  
  table2.time <- sapply(node_results, function(x){       
    x$test_time })
  single_testtime <- mean(table2.time)
  
  power2.sim[[loopIndex]] <- list(Power = Power2, realTau = c(r.sim,r.sim/2,r.sim/4), smooth = smooth, single_test_sec = single_testtime)
  
  save(power2.sim, file=fileName) # Auto Save
  
}

stopCluster(cluster)