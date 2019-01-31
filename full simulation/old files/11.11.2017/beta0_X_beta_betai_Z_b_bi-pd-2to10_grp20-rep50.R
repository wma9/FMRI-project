# For your information, briefly, the model is that, as below,
# Y_ij = intercept + dummyX_ij * (beta + beta_i) + Z_ij * (b + b_i) + error_ij, 
# where 'beta' and 'b' are fixed effects while 'b and 'b_i' are random slope effects.
# We conduct the test that H_0: variance of b_i = 0 v.s. H_a: variance of b_i > 0 through exactRLR.
# In this simulation, we use parallel computing in order to reduce the performing time.

# Setups
library(parallel)
simRep <- 1000 # Replication times in one simulation
pvalue.true <- .05 # Testing type I error 
b.var <- c(0) # The set of varaince of random covariates b as random slope

# Below is the function defined for each node in parallel
run_one_sample <- function(iter){
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  
  set.seed(iter)
  nGroup <- 20 # Group number
  nRep.sim <- 50 # Duplication in each group n_i 
  epsilon.sd <- 1 # or 0.5
  intercept.true <- 0.5
  # nRandCovariate <- 
  b.sim <- 2
  beta.sim <- 2
  betaVar.sim <- 1
  z.mean <- 1
  totalN <- nGroup * nRep.sim
  
  z.var <- c(1)
  ID.sim <- rep(1:nGroup, each = nRep.sim)
  error.sim <- rnorm(n = totalN, mean = 0, sd = epsilon.sd)
  z.sim <- mapply(rnorm, totalN, z.mean, rep(sqrt(z.var), nRandCovariate))
  bV.sim <- mapply(rnorm, nGroup, b.sim, rep(sqrt(r.sim), nRandCovariate))
  bV.sim <- bV.sim[rep(1:nrow(bV.sim), each = nRep.sim), ]
  bV.sim <- bV.sim * z.sim
  bV.sim <- rowSums(bV.sim)
  
  betaV.sim <- mapply(rnorm, nGroup, beta.sim, rep(sqrt(betaVar.sim), 1))
  betaV.sim <- betaV.sim[rep(1:nrow(betaV.sim), each = nRep.sim), ]

  dummyX <- 2 * rbinom(n = totalN, size = 1, prob = 0.5) - 1 # Shift dummyX
  Y.sim <- intercept.true + bV.sim + dummyX * beta.sim + error.sim # NEW add 'dummyX'
  
  
  ID = ID.sim
  Y = Y.sim
  npc = nRandCovariate

  designMatrix.pdIdent <- data.frame(rating = Y, 
                             temp = factor(dummyX, labels=c("warm", "hot")),
                             ID = as.factor(ID), 
                             a.score = z.sim)
  # 'lme' model with 'pdIdent'
  fullReml.pdIdent <- NA
  noAScoreReml.pdIdent <- NA
  notempReml.pdIdent <- NA
  

  # if(npc == 1){
  #   fullReml.pdIdent <- lme(fixed = Y ~ 1 + dummyX + z.sim, 
  #                           random = list(ID = pdIdent(~ 0 + dummyX), 
  #                                         ID = pdIdent(~ 0 + z.sim)), 
  #                           data = designMatrix, 
  #                           control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
  #   noZReml.pdIdent <- lme(fixed = Y ~ 1 + dummyX + z.sim,
  #                          random = list(ID = pdIdent(~ 0 + dummyX)),
  #                          data = designMatrix,
  #                          control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
  #   noDummyXReml.pdIdent <- lme(fixed = Y ~ 1 + dummyX + z.sim,
  #                               random = list(ID = pdIdent(~ 0 + z.sim)),
  #                               data = designMatrix,
  #                               control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
  # }else
    if(npc == 2){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1 + a.score.2)), 
                            data = designMatrix.pdIdent,
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2, 
                              random = list(ID = pdIdent(~ 0 + a.score.1 + a.score.2)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else if(npc == 3){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3)), 
                            data = designMatrix.pdIdent,
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3, 
                              random = list(ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else if(npc == 4){    
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + a.score.4, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + a.score.4)), 
                            data = designMatrix.pdIdent,
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + a.score.4, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + a.score.4, 
                              random = list(ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + a.score.4)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else if(npc == 5){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                              a.score.4 + a.score.5, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                         a.score.4 + a.score.5)), 
                            data = designMatrix.pdIdent,
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                  a.score.4 + a.score.5, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                a.score.4 + a.score.5, 
                              random = list(ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                           a.score.4 + a.score.5)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else if(npc == 6){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                              a.score.4 + a.score.5 + a.score.6, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                         a.score.4 + a.score.5 + a.score.6)), 
                            data = designMatrix.pdIdent,
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                  a.score.4 + a.score.5 + a.score.6, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                a.score.4 + a.score.5 + a.score.6, 
                              random = list(ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                           a.score.4 + a.score.5 + a.score.6)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else if(npc == 7){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                              a.score.4 + a.score.5 + a.score.6 + a.score.7, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                         a.score.4 + a.score.5 + a.score.6 + a.score.7)), 
                            data = designMatrix.pdIdent,
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                  a.score.4 + a.score.5 + a.score.6 + a.score.7, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                a.score.4 + a.score.5 + a.score.6 + a.score.7, 
                              random = list(ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                           a.score.4 + a.score.5 + a.score.6 + a.score.7)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else if(npc == 8){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                              a.score.4 + a.score.5 + a.score.6 + a.score.7 + a.score.8, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                         a.score.4 + a.score.5 + a.score.6 + 
                                                         a.score.7 + a.score.8)), 
                            data = designMatrix.pdIdent,
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                  a.score.4 + a.score.5 + a.score.6 + a.score.7 + a.score.8, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                a.score.4 + a.score.5 + a.score.6 + a.score.7 + a.score.8, 
                              random = list(ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                           a.score.4 + a.score.5 + a.score.6 + 
                                                           a.score.7 + a.score.8)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else if(npc == 9){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                              a.score.4 + a.score.5 + a.score.6 + a.score.7 + a.score.8 + a.score.9, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                         a.score.4 + a.score.5 + a.score.6 + 
                                                         a.score.7 + a.score.8 + a.score.9)), 
                            data = designMatrix.pdIdent,
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                  a.score.4 + a.score.5 + a.score.6 + a.score.7 + a.score.8 + a.score.9, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1 + a.score.2 + a.score.3 + 
                                a.score.4 + a.score.5 + a.score.6 + a.score.7 + a.score.8 + a.score.9, 
                              random = list(ID = pdIdent(~ 0 + a.score.1 + a.score.2 + a.score.3 + 
                                                           a.score.4 + a.score.5 + a.score.6 + 
                                                           a.score.7 + a.score.8 + a.score.9)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else{
    #    additive0.sim <- paste(1:pca.npc, collapse = " + a.score.")
    #    modelFix.sim <- as.formula(paste("rating ~ 1 + temp + a.score.", 
    #                                     additive0.sim, 
    #                                     sep = ""))
    #    modelRan.sim <- as.formula(paste("~ 0 + a.score.", 
    #                                     additive0.sim, 
    #                                     sep = ""))
    # fullReml.pdIdent <- lme(fixed = modelFix.sim, 
    #                         random = list(ID = pdIdent(~ 0 + temp), 
    #                                       ID = pdIdent(modelRan.sim)), 
    #                         data = designMatrix.pdIdent)
    # noAScoreReml.pdIdent <- lme(fixed = modelFix.sim, 
    #                             random = list(ID = pdIdent(~ 0 + temp)), 
    #                             data = designMatrix.pdIdent)
    # notempReml.pdIdent <- lme(fixed = modelFix.sim, 
    #                           random = list(ID = pdIdent(modelRan.sim)), 
    #                           data = designMatrix.pdIdent)
  }
  
  tests1 <- exactRLRT(notempReml.pdIdent) # , nsim = 100000

  
  # 'lmer' model
  designMatrix.lmm <- designMatrix.pdIdent
  
  additive0.sim <- paste(1:npc, collapse = " + a.score.")
  additive.sim <- paste(1:npc, collapse = " | ID) + (0 + a.score.")
  # Confusion of modifying
  model.sim <- as.formula(paste("rating ~ 1 + temp + a.score.", 
                                additive0.sim, 
                                " + (0 + a.score.", 
                                additive.sim, 
                                " | ID)", 
                                sep = ""))
  fullReml <- lmer(model.sim, data = designMatrix.lmm)
  
  tests2 <- list()
  for(i in 1:npc){
    ii <- paste("a.score.", i, sep = "")
    f0 <- as.formula(paste(" . ~ . - (0 + ", ii, "| ID)"))
    m0 <- update(fullReml, f0)
    f.slope <- as.formula(paste("rating ~ 1 + temp + a.score.", additive0.sim, " + (0 +", ii, " | ID)", 
                                sep = ""))
    m.slope <- lmer(f.slope, data = designMatrix.lmm)
    tests2[[i]] <- exactRLRT(m.slope, fullReml, m0)
  }
  multiTest1 <- sapply(tests2, function(x) {
    c(statistic = x$statistic[1],
      "p-value" = x$p[1])})
  pvalues.bonf <- p.adjust(multiTest1[2,], "bonferroni")
  
  
  # fullReml <- lmer(Y ~ 1 + dummyX + z.sim + (0 + dummyX | ID) + (0 + z.sim | ID), data = designMatrix) ## + (1 | temp : ID)
  # m0 <- update(fullReml,  . ~ . - (0 + z.sim | ID))
  # m.slope <- update(fullReml, . ~ . -  (0 + dummyX | ID))
  # tests2 <- list()
  # tests2[[1]] <- exactRLRT(m.slope, fullReml, m0) # , nsim = 100000
  # 
  # 
  # tests2 <- sapply(tests2, function(x){c(statistic = x$statistic[1], "p-value" = x$p[1])})
  # pvalues.bonf <- p.adjust(tests2[2,], "bonferroni")
  
  return(list(realTau = r.sim, 
              pvalue = tests1$p[1], 
              pvalues.bonf = pvalues.bonf, 
              tests1 = tests1, 
              tests2 = multiTest1))
}

# Setup parallel
cores <- detectCores()
cluster <- makeCluster(cores)
clusterSetRNGStream(cluster, 20170822)

for(nRandCovariate in 2:9){  # START out-outer loop
  clusterExport(cluster, c("nRandCovariate")) # casting the coefficient parameter on the random effects' covariance function
  fileName <- paste("beta0_X_beta_betai_Z_b_bi-pd-", nRandCovariate,".RData", sep = "") # Saving file's name
  
  
# run the simulation
loopIndex <- 1
resultDoubleList.sim <- list()
power1.sim <- list()
power2.sim <- list()
for(r.sim in b.var){  # START outer loop
  clusterExport(cluster, c("r.sim")) # casting the coefficient parameter on the random effects' covariance function
  
  node_results <- parLapply(cluster, 1:simRep, run_one_sample)
  
  result1.sim <- lapply(node_results, function(x) {list(realTau = x$realTau, 
                                                        pvalue = x$pvalue)})
  result2.sim <- lapply(node_results, function(x) {list(realTau = x$realTau, 
                                                        pvalues.bonf = x$pvalues.bonf)})
  resultDoubleList.sim[[loopIndex]] <- node_results
  
  save.image(file=fileName) # Auto Save
  
  table1.sim <- sapply(result1.sim, function(x) {       
    c(sens = (sum(x$pvalue <= pvalue.true) > 0))})
  Power1 <- mean(table1.sim)
  cat("nRandCovariate: ", nRandCovariate, fill = TRUE)
  cat("Power1: ", Power1, fill = TRUE)
  power1.sim[[loopIndex]] <- list(Power = Power1, realTau = r.sim) 
  
  table2.sim <- sapply(result2.sim, function(x) {       
    c(overall.sens = (sum(x$pvalues.bonf <= pvalue.true) > 0))})
  Power2 <- mean(table2.sim)
  cat("Power2: ", Power2, fill = TRUE)
  power2.sim[[loopIndex]] <- list(Power = Power2, realTau = r.sim)
  
  loopIndex <- loopIndex + 1
} # End outer loop

save.image(file=fileName) # Auto Save


par(mfrow=c(2,1))
# Histogram plots
hist(sapply(result1.sim, function(x) x$pvalue), 
     main = "Histogram of p-value for lme model", 
     xlab = "p-value")

hist(sapply(result2.sim, function(x) x$pvalues.bonf), 
     main = "Histogram of p-value for lmer model", 
     xlab = "p-value")

# hist(sapply(resultDoubleList.sim[[1]], function(x) (x$tests1)$statistic[1]),
#      breaks = (0:110)/10, 
#      main = "Histogram of test-statistic for lme model",
#      xlab = "Test Statistics")
# 
# hist(sapply(resultDoubleList.sim[[1]], function(x) (x$tests2)[1,1]), 
#      breaks = (0:100)/10, 
#      main = "Histogram of test-statistic for lmer model",
#      xlab = "Test Statistics")

} # End out-outer loop

stopCluster(cluster)


