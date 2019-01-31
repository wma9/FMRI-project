# For your information, briefly, the model is that, as below,
# Y_ij = intercept + dummyX_ij * (beta + beta_i) + Z_ij * (b + b_i) + error_ij, 
# where 'beta' and 'b' are fixed effects while 'b and 'b_i' are random slope effects.
# We conduct the test that H_0: variance of b_i = 0 v.s. H_a: variance of b_i > 0 through exactRLR.
# In this simulation, we use parallel computing in order to reduce the performing time.

# Setups
library(parallel)
simRep <- 1000 # Replication times in one simulation
pvalue.true <- .05 # Testing type I error 
fileName <- "beta0_X_beta_betai_Z_b_bi-pd.RData" # Saving file's name
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
  nRandCovariate <- 1
  b.sim <- 2
  beta.sim <- 2
  betaVar.sim <- 1
  z.mean <- 1
  totalN <- nGroup * nRep.sim
  
  z.var <- c(1)
  ID.sim <- rep(1:nGroup, each = nRep.sim)
  error.sim <- rnorm(n = totalN, mean = 0, sd = epsilon.sd)
  z.sim <- mapply(rnorm, totalN, z.mean, sqrt(z.var))
  bV.sim <- mapply(rnorm, nGroup, b.sim, rep(sqrt(r.sim), nRandCovariate))
  bV.sim <- bV.sim[rep(1:nrow(bV.sim), each = nRep.sim), ]
  bV.sim <- bV.sim * z.sim
  bV.sim <- rowSums(bV.sim)
  
  betaV.sim <- mapply(rnorm, nGroup, beta.sim, rep(sqrt(betaVar.sim), 1))
  betaV.sim <- betaV.sim[rep(1:nrow(betaV.sim), each = nRep.sim), ]

  dummyX <- 2 * rbinom(n = totalN, size = 1, prob = 0.5) - 1 # Shift dummyX
  Y.sim <- (intercept.true + bV.sim) + dummyX * betaV.sim + error.sim # NEW add 'dummyX'
  
  
  ID = ID.sim
  Y = Y.sim
  dummyX = dummyX
  npc = nRandCovariate
  npc = nRandCovariate
  z.sim = z.sim
  
  designMatrix <- data.frame(Y = Y, 
                             dummyX = factor(dummyX, labels=c("warm", "hot")),
                             ID = as.factor(ID), 
                             z.sim = z.sim)
  # 'lme' model with 'pdIdent'
  fullReml.pdIdent <- NA
  noZReml.pdIdent <- NA
  noDummyXReml.pdIdent <- NA
  
  if(npc == 1){
    fullReml.pdIdent <- lme(fixed = Y ~ 1 + dummyX + z.sim, 
                            random = list(ID = pdIdent(~ 0 + dummyX), 
                                          ID = pdIdent(~ 0 + z.sim)), 
                            data = designMatrix, 
                            control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noZReml.pdIdent <- lme(fixed = Y ~ 1 + dummyX + z.sim,
                           random = list(ID = pdIdent(~ 0 + dummyX)),
                           data = designMatrix,
                           control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
    noDummyXReml.pdIdent <- lme(fixed = Y ~ 1 + dummyX + z.sim,
                                random = list(ID = pdIdent(~ 0 + z.sim)),
                                data = designMatrix,
                                control = lmeControl(msVerbose = TRUE, opt = 'optim', singular.ok=TRUE, returnObject=TRUE))
  }
  
  tests1 <- exactRLRT(noDummyXReml.pdIdent, fullReml.pdIdent, noZReml.pdIdent) # , nsim = 100000

  
  # 'lmer' model
  fullReml <- lmer(Y ~ 1 + dummyX + z.sim + (1 | dummyX : ID) + (0 + z.sim | ID), data = designMatrix) # (0 + dummyX | ID)
  m0 <- update(fullReml,  . ~ . - (0 + z.sim | ID))
  m.slope <- update(fullReml, . ~ . -  (1 | dummyX : ID))
  tests2 <- list()
  tests2[[1]] <- exactRLRT(m.slope, fullReml, m0) # , nsim = 100000

  
  tests2 <- sapply(tests2, function(x){c(statistic = x$statistic[1], "p-value" = x$p[1])})
  pvalues.bonf <- p.adjust(tests2[2,], "bonferroni")
  
  return(list(realTau = r.sim, 
              pvalue = tests1$p[1], 
              pvalues.bonf = pvalues.bonf, 
              tests1 = tests1, 
              tests2 = tests2))
}

# Setup parallel
cores <- detectCores()
cluster <- makeCluster(cores)
clusterSetRNGStream(cluster, 20170822)

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
  power1.sim[[loopIndex]] <- list(Power = Power1, realTau = r.sim) 
  
  table2.sim <- sapply(result2.sim, function(x) {       
    c(overall.sens = (sum(x$pvalues.bonf <= pvalue.true) > 0))})
  Power2 <- mean(table2.sim)
  power2.sim[[loopIndex]] <- list(Power = Power2, realTau = r.sim)
  
  loopIndex <- loopIndex + 1
} # End outer loop

save.image(file=fileName) # Auto Save
stopCluster(cluster)


par(mfrow=c(2,2))
# Histogram plots
hist(sapply(result1.sim, function(x) x$pvalue), 
     main = "Histogram of p-value for lme model", 
     xlab = "p-value")

hist(sapply(result2.sim, function(x) x$pvalues.bonf), 
     main = "Histogram of p-value for lmer model", 
     xlab = "p-value")

hist(sapply(resultDoubleList.sim[[1]], function(x) (x$tests1)$statistic[1]),
     breaks = (0:110)/10, 
     main = "Histogram of test-statistic for lme model",
     xlab = "Test Statistics")

hist(sapply(resultDoubleList.sim[[1]], function(x) (x$tests2)[1,1]), 
     breaks = (0:100)/10, 
     main = "Histogram of test-statistic for lmer model",
     xlab = "Test Statistics")


