library(parallel)

sample.sim <- 5000
pvalue.total <- .05
fileName <- "simu_var0_20_100_5000.RData"
rRange.sim <- c(0)

run_one_sample <- function(numSim){
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  
  set.seed(numSim)
  
  nTime.sim <- 80 # J
  nSubj.sim <- 20 # 200 # I the number of curves
  nRep.sim <- 100 # 50 # datasets for each covariance function 
  sd.epsilon <- 1 # or 0.5
  true.delta <- 0.5
  true.npc <- 1
  percent <- 0.99
  SNR.sim <- .5 # 3
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
  
  lambdaVec.sim <- c(lambda.sim(1))
  amply.sim <- sd.epsilon^2 * SNR.sim * 2 / sum(lambdaVec.sim)
  timeSlot.sim <- (1:nTime.sim)/nTime.sim
  ID.sim <- rep(1:nSubj.sim, each = nRep.sim)
  
  
  
  
  
  
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
  
  phi.sim <- matrix(data = c(phi.func(timeSlot.sim, 1)), 
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
  
  pca.npc = true.npc
  npc = true.npc
  a.score = ascore.sim
  
  

  # Single test
  designMatrix.pdIdent <- data.frame(rating = rating, 
                                     temp = factor(temp, labels=c("warm", "hot")),
                                     ID = as.factor(ID), 
                                     a.score = a.score)
  fullReml.pdIdent <- NA
  noAScoreReml.pdIdent <- NA
  notempReml.pdIdent <- NA
  
  if(pca.npc == 1){
    fullReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score, 
                            random = list(ID = pdIdent(~ 0 + temp), 
                                          ID = pdIdent(~ 0 + a.score)), 
                            data = designMatrix.pdIdent, 
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noAScoreReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score,
                                random = list(ID = pdIdent(~ 0 + temp)),
                                data = designMatrix.pdIdent,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    notempReml.pdIdent <- lme(fixed = rating ~ 1 + temp + a.score,
                              random = list(ID = pdIdent(~ 0 + a.score)),
                              data = designMatrix.pdIdent,
                              control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
  }else{

  }
  
  
  #test for subject specific slopes in brain given temp:
  tests1 <- exactRLRT(notempReml.pdIdent, fullReml.pdIdent, noAScoreReml.pdIdent) # , nsim = 100000
  

  
  
  
  
  # Multiple tests
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
  fullReml <- lmer(rating ~ 1 + temp + a.score + (0 + temp | ID) + (0 + a.score | ID), data = designMatrix.lmm)
  
  tests2 <- list()
  for(i in 1:pca.npc){
    ii <- paste("a.score.", i, sep = "")
    f0 <- as.formula(paste(" . ~ . - (0 + ", ii, "| ID)"))
    m0 <- update(fullReml,  . ~ . - (0 + a.score | ID))
    f.slope <- as.formula(paste("rating ~ 1 + temp + a.score.", additive0.sim, " + (0 +", ii, " | ID)", 
                                sep = ""))
    m.slope <- update(fullReml, . ~ . -  (0 + temp | ID))
    tests2[[i]] <- exactRLRT(m.slope, fullReml, m0) # , nsim = 100000
  }
  multiTest1 <- sapply(tests2, function(x) {
    c(statistic = x$statistic[1],
      "p-value" = x$p[1])})
  pvalues.bonf <- p.adjust(multiTest1[2,], "bonferroni")
  
  ###################################################################################
  
  print(paste("r.sim = ", r.sim, ", iteration: ", numSim, ", time: ", proc.time() - ptm))
  
  
  return(list(realTau = r.sim, 
              pvalue = tests1$p[1], 
              pvalues.bonf = pvalues.bonf, 
              npc = npc, 
              tests1 = tests1, 
              tests2 = multiTest1))
}








cores <- detectCores()
cluster <- makeCluster(cores)
# cluster <- makeCluster(as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")))
clusterSetRNGStream(cluster, 20170822)


loopIndex <- 1
resultDoubleList.sim <- list()
npc.stack <- c()
powerSingle.sim <- list()
powerMulti.sim <- list()
for(r.sim in rRange.sim){  # START Primary Loop

  
  clusterExport(cluster, c("r.sim"))
  
  node_results <- parLapply(cluster, 1:sample.sim, run_one_sample)
  
  resultSingle.sim <- lapply(node_results, function(x) {list(realTau = x$realTau, 
                                                             pvalue = x$pvalue)})
  resultMulti.sim <- lapply(node_results, function(x) {list(realTau = x$realTau, 
                                                            pvalues.bonf = x$pvalues.bonf)})
  npc.stack <- c(npc.stack, sapply(node_results, function(x) x$npc))
  resultDoubleList.sim[[loopIndex]] <- node_results
  

  
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





save.image(file=fileName)


stopCluster(cluster)





hist(sapply(resultSingle.sim, function(x) x$pvalue), main = "Plot of p-value for lme model")

hist(sapply(resultMulti.sim, function(x) x$pvalues.bonf), main = "Plot of p-value for lmer model")

hist(sapply(resultDoubleList.sim[[1]], function(x) (x$tests1)$statistic[1]),
     breaks = (0:100)/10, 
     main = "Plot of test-statistic for lme model")

hist(sapply(resultDoubleList.sim[[1]], function(x) (x$tests2)[1,1]), 
     breaks = (0:100)/10, 
     main = "Plot of test-statistic for lmer model")


