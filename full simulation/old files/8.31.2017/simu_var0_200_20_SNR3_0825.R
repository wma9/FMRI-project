# We assign each sampling within each covariance to each core
#

library(refund)
library(lme4)
library(nlme)
library(arm)
library(RLRsim)
library(parallel)

sample.sim <- 500
pvalue.total <- .05
fileName <- "simu_var0_200_20_SNR3_0825.RData"
rRange.sim <- c(0)

run_one_sample <- function(numSim){
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  
  nTime.sim <- 80 # J
  nSubj.sim <- 200 # 50 # I the number of curves
  nRep.sim <- 20 # 50 # datasets for each covariance function 
  sd.epsilon <- 1 # or 0.5
  true.delta <- 0.5
  true.npc <- 3
  percent <- 0.95
  SNR.sim <- 3 # 5
  beta.sim <- 2
  a.mean <- 1
  n.sim <- nSubj.sim * nRep.sim
  
  lambda.sim <- function(degree) {
    return(0.5^(degree - 1))
  }
  
  psi.fourier <- function(t, degree) {
    result <- NA
    if(degree == 1){
      result <- sqrt(2) * sin(2*pi*t)
    }else if(degree == 2){
      result <- sqrt(2) * cos(4*pi*t)
    }else if(degree == 3){
      result <- sqrt(2) * sin(4*pi*t)
    }
    return(result)
  }
  
  psi.polyn <- function(t, degree) {
    result <- NA
    if(degree == 1){
      result <- sqrt(3) * (2*t - 1)
    }else if(degree == 2){
      result <- sqrt(5) * (6*t^2 - 6*t +1)
    }else if(degree == 3){
      result <- sqrt(7) * (20*t^3 - 30*t^2 + 12*t - 1)
    }
    return(result)
  }
  
  lambdaVec.sim <- c(lambda.sim(1), lambda.sim(2), lambda.sim(3))
  amply.sim <- sd.epsilon^2 * SNR.sim * 2 / sum(lambdaVec.sim)
  timeSlot.sim <- (1:nTime.sim)/nTime.sim
  ID.sim <- rep(1:nSubj.sim, each = nRep.sim)
  
  
  
  set.seed(numSim)
  
  
  
  ptm <- proc.time()
  
  eta.sim <- rnorm(n = n.sim, mean = 0, sd = sd.epsilon)
  ascore.sim <- mapply(rnorm, n.sim, a.mean, sqrt(lambdaVec.sim))
  
  b.sim <- mapply(rnorm, nSubj.sim, beta.sim, rep(sqrt(r.sim), true.npc))
  b.sim <- b.sim[rep(1:nrow(b.sim), each = nRep.sim), ]
  b.sim <- b.sim * ascore.sim
  b.sim <- rowSums(b.sim)
  
  Y.sim <- (true.delta + b.sim)*amply.sim + eta.sim
  
  z.sim <- rbinom(n = n.sim, size = 1, prob = 0.5)
  
  if(1 == 1){
    phi.func <- psi.fourier
  }else{
    phi.func <- psi.polyn
  }
  
  phi.sim <- matrix(data = c(phi.func(timeSlot.sim, 1), 
                             phi.func(timeSlot.sim, 2), 
                             phi.func(timeSlot.sim, 3)), 
                    nrow = true.npc, 
                    ncol = nTime.sim, 
                    byrow = TRUE)
  Mt.sim <- ascore.sim %*% phi.sim
  Mt.sim <- amply.sim* Mt.sim + matrix(rnorm((n.sim)*nTime.sim, mean=0, sd = sd.epsilon), n.sim, nTime.sim) 
  
  ##########################################################################
  
  Y = Mt.sim
  ID = ID.sim
  rating = Y.sim
  temp = z.sim
  Y <- t(scale(t(Y))) # standardized
  t <- 1:nTime.sim
  knots <- 5 # previous setting 10
  p <- 5  # previous setting p <- 7
  
  results <- fpca.face(Y, center = TRUE, argvals = t, knots = knots, pve = percent, p = p, lambda = 0) # pve need to be chosen!
  npc <- results$npc
  score <- results$scores
  pca.npc <- npc
  a.score <- score[, 1:pca.npc]
  
  # # Single test
  # a.score.single <- sqrt((rowSums(a.score^2))^(-1))
  # designMatrix.lmm <- data.frame(rating = rating, 
  #                                temp = factor(temp, labels=c("warm", "hot")),
  #                                ID = as.factor(ID), 
  #                                a.score = a.score,
  #                                a.score.single = a.score.single)
  # 
  # if(1 == 1){
  #   fullReml.single <- lmer(rating ~ 1 + temp + a.score 
  #                           # + (1 | ID) 
  #                           + (1 | temp : ID) 
  #                           #                   + (0 + temp | ID) 
  #                           + (0 + a.score.single | ID), data = designMatrix.lmm)
  #   
  #   noAScoreReml.single <- lmer(rating ~ 1 + temp + a.score 
  #                               # + (1 | ID) 
  #                               + (1 | temp : ID), data = designMatrix.lmm)
  #   notempReml.single <- lmer(rating ~ 1 + temp + a.score 
  #                             # + (1 | temp : ID) 
  #                             + (0 + a.score.single | ID), data = designMatrix.lmm)
  # }
  # 
  # #test for subject specific slopes in brain given temp:
  # tests1 <- exactRLRT(notempReml.single, fullReml.single, noAScoreReml.single)
  # # #test for subject specific slopes in temp No brain:
  # # exactRLRT(noAScoreReml.single)
  # # resultSingle.sim[[numSim]] <- list(realTau = r.sim, pvalue = tests1$p[1])
  
  
  
  # Single test
  designMatrix.pdIdent <- data.frame(rating = rating, 
                                     temp = factor(temp, labels=c("warm", "hot")),
                                     ID = as.factor(ID), 
                                     a.score = a.score)
  fullReml.pdIdent <- NA
  noAScoreReml.pdIdent <- NA
  notempReml.pdIdent <- NA
  
  if(pca.npc == 1){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score.1)), 
                            data = designMatrix.pdIdent, 
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1, 
                                random = list(ID = pdIdent(~ 0 + temp)), 
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score.1, 
                              random = list(ID = pdIdent(~ 0 + a.score.1)), 
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    
  }else if(pca.npc == 2){
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
    
  }else if(pca.npc == 3){
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
    
  }else if(pca.npc == 4){    
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
    
  }else if(pca.npc == 5){
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
    
  }else if(pca.npc == 6){
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
    
  }else if(pca.npc == 7){
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
    
  }else if(pca.npc == 8){
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
    
  }else if(pca.npc == 9){
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
  
  
  #test for subject specific slopes in brain given temp:
  tests1 <- exactRLRT(notempReml.pdIdent, fullReml.pdIdent, noAScoreReml.pdIdent) 
  # When testing, it recalls the model, but "lme" lost all components!
  # #test for subject specific slopes in temp No brain:
  # exactRLRT(noAScoreReml.single)
  # resultSingle.sim[[numSim]] <- list(realTau = r.sim, pvalue = tests1$p[1])
  
  
  
  
  
  # Multiple tests
  # designMatrix.lmm <- data.frame(rating = rating, 
  #                                temp = factor(temp, labels=c("warm", "hot")),
  #                                ID = as.factor(ID), 
  #                                a.score = a.score)
  designMatrix.lmm <- designMatrix.pdIdent
  
  additive0.sim <- paste(1:pca.npc, collapse = " + a.score.")
  additive.sim <- paste(1:pca.npc, collapse = " | ID) + (0 + a.score.")
  # Confusion of modifying
  model.sim <- as.formula(paste("rating ~ 1 + temp + a.score.", 
                                additive0.sim, 
                                " + (1 | temp : ID) + (0 + a.score.", 
                                additive.sim, 
                                " | ID)", 
                                sep = ""))
  fullReml <- lmer(model.sim, data = designMatrix.lmm)
  
  tests2 <- list()
  for(i in 1:pca.npc){
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
  
  # resultMulti.sim[[numSim]] <- list(realTau = r.sim, pvalues.bonf = pvalues.bonf)
  ###################################################################################
  
  print(paste("r.sim = ", r.sim, ", iteration: ", numSim, ", time: ", proc.time() - ptm))
  
  
  return(list(realTau = r.sim, 
              pvalue = tests1$p[1], 
              pvalues.bonf = pvalues.bonf, 
              npc = npc, 
              tests1 = tests1, 
              tests2 = multiTest1))
}








# cores <- detectCores()
# cluster <- makeCluster(cores)
cluster <- makeCluster(as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")))
clusterSetRNGStream(cluster, 20170822)


loopIndex <- 1
resultDoubleList.sim <- list()
npc.stack <- c()
powerSingle.sim <- list()
powerMulti.sim <- list()
for(r.sim in rRange.sim){  # START Primary Loop
  # resultSingle.sim <- list()
  # resultMulti.sim <- list()
  
  # for(numSim in 1:sample.sim){  # START Secondary Loop
  #   return(list(realTau = r.sim, pvalue = tests1$p[1], pvalues.bonf = pvalues.bonf))
  # } # END LOOP
  
  clusterExport(cluster, c("r.sim"))
  
  node_results <- parLapply(cluster, 1:sample.sim, run_one_sample)
  
  resultSingle.sim <- lapply(node_results, function(x) {list(realTau = x$realTau, 
                                                             pvalue = x$pvalue)})
  resultMulti.sim <- lapply(node_results, function(x) {list(realTau = x$realTau, 
                                                            pvalues.bonf = x$pvalues.bonf)})
  npc.stack <- c(npc.stack, sapply(node_results, function(x) x$npc))
  resultDoubleList.sim[[loopIndex]] <- node_results
  
  # test_centers <- 2:6
  # node_results <- parLapply(cluster, test_centers, function(n_centers) run_one_sample)
  # final_results <- by(node_results, test_centers, median)
  
  
  save.image(file=fileName) # Auto Save
  
  
  
  tableSingle.sim <- sapply(resultSingle.sim, function(x) {       
    c(sens = (sum(x$pvalue <= pvalue.total) > 0))})
  PowerSingle <- mean(tableSingle.sim)
  powerSingle.sim[[loopIndex]] <- list(Power = PowerSingle, realTau = r.sim) 
  
  tableMulti.sim <- sapply(resultMulti.sim, function(x) {       
    c(overall.sens = (sum(x$pvalues.bonf <= pvalue.total) > 0))})
  PowerMulti <- mean(tableMulti.sim)
  powerMulti.sim[[loopIndex]] <- list(Power = PowerMulti, realTau = r.sim)
  
  loopIndex <- loopIndex + 1
} # END LOOP


# powertableSingle.sim <- sapply(powerSingle.sim, function(x) {
#   c(Power = x$Power,
#     realTau = x$realTau)})
# plot(powertableSingle.sim[2, ],powertableSingle.sim[1, ],type="l",lwd = 3,
#      ylim = c(0,1),
#      ylab = "Power",
#      xlab = "Random Effect Covariance", 
#      main = "Plot of Power against Random Effect Covariance by Single test")
# 
# 
# 
# powertableMulti.sim <- sapply(powerMulti.sim, function(x) {
#   c(Power = x$Power,
#     realTau = x$realTau)})
# plot(powertableMulti.sim[2, ],powertableMulti.sim[1, ],type="l",lwd = 3,
#      ylim = c(0,1),
#      ylab = "Power",
#      xlab = "Random Effect Covariance", 
#      main = "Plot of Power against Random Effect Covariance by Multiple tests")


save.image(file=fileName)


stopCluster(cluster)



