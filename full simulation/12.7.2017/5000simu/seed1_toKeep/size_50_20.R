# This code is for size test(suedo code)
# For your information, briefly, the model is that, as below,
# Y_ij = intercept + dummyX_ij * (beta + beta_i) + Z_ij * (b + b_i) + error_ij, 
# where 'beta' and 'b' are fixed effects while 'b and 'b_i' are random slope effects.
# We conduct the test that H_0: variance of b_i = 0 v.s. H_a: variance of b_i > 0 through exactRLR.
# In this simulation, we use parallel computing in order to reduce the performing time.

# Setups
library(parallel)
simRep <- 20000 # Replication times in one simulation
pvalue.true <- .05 # Testing type I error 
b.var <- c(0) # The set of varaince of random covariates b as random slope
cores <- 4

# Below is the function defined for each node in parallel
run_one_sample <- function(iter){
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  # library(MASS)
  
  set.seed(iter)
  nGroup <- 50 # Group number
  nRep.sim <- 20 # Duplication in each group n_i 
  epsilon.sd <- 1 # or 0.5
  intercept.true <- 0.5
  # nRandCovariate <- 2
  # r.sim <- 0
  b.sim <- 2
  beta.sim <- 2
  betaVar.sim <- 1
  z.mean <- 0
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
  betaV2.sim <- mapply(rnorm, nGroup, beta.sim, rep(sqrt(betaVar.sim), 1))
  betaV2.sim <- betaV2.sim[rep(1:nrow(betaV2.sim), each = nRep.sim), ]
  
  dummyX <- rbinom(n = totalN, size = 1, prob = 0.5) # Shift dummyX
  Y.sim <- (intercept.true + bV.sim) + dummyX * betaV.sim + (dummyX - 1) * betaV2.sim + error.sim # NEW add 'dummyX'
  
  
  ID = ID.sim
  Y = Y.sim
  npc = nRandCovariate
  
  ##12.7.2017
  dummyX <- cbind(dummyX, -dummyX + 1)
  z.sim.uni = c()
  ID.uni <- c(rbind(matrix(1:(nGroup*nRandCovariate), 
                           nrow = nRandCovariate, 
                           ncol = nGroup), 
                    matrix(0, nrow = nRep.sim - nRandCovariate, ncol = nGroup)))
  
  for(k in 1:nGroup){
    svd <- svd(z.sim[((k-1)*nRep.sim+1):(k*nRep.sim), ] %*% t(z.sim[((k-1)*nRep.sim+1):(k*nRep.sim), ])) #SVD on A_i
    u.tra <- t(svd$v)
    u <- svd$u
    d <- (svd$d)[1:nRandCovariate]
    # u <- cbind(u, Null(u))
    Y[((k-1)*nRep.sim+1):(k*nRep.sim)] <- u.tra %*% Y[((k-1)*nRep.sim+1):(k*nRep.sim)]
    dummyX[((k-1)*nRep.sim+1):(k*nRep.sim), ] <- u.tra %*% dummyX[((k-1)*nRep.sim+1):(k*nRep.sim), ]
    z.sim[((k-1)*nRep.sim+1):(k*nRep.sim), ] <- rbind(u.tra[1:nRandCovariate, ] %*% z.sim[((k-1)*nRep.sim+1):(k*nRep.sim), ], 
                                                      matrix(0, 
                                                             nrow = nRep.sim - nRandCovariate, 
                                                             ncol = nRandCovariate))
    z.sim.uni <- c(z.sim.uni, sqrt(d), rep(0, nRep.sim - nRandCovariate))
    
  }
  ##12.7.2017
  
  
  designMatrix.pdIdent <- data.frame(rating = Y, 
                                     temp.1 = dummyX[, 1],
                                     temp.2 = dummyX[, 2],
                                     ID = as.factor(ID),
                                     ID.uni = as.factor(ID.uni),
                                     a.score = z.sim,
                                     z.sim.uni = z.sim.uni)
  
  
  # 'lmer' model
  designMatrix.lmm <- designMatrix.pdIdent
  
  additive0.sim <- paste(1:npc, collapse = " + a.score.")
  additive.sim <- paste(1:npc, collapse = " | ID) + (0 + a.score.")
  # Confusion of modifying
  model.sim <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + a.score.", 
                                additive0.sim, 
                                " + (0 + temp.1 | ID) + (0 + temp.2 | ID) + (0 + z.sim.uni | ID.uni)", 
                                sep = ""))
  fullReml <- lmer(model.sim, data = designMatrix.lmm)
  
  # tests2 <- list()
  # for(i in 1:npc){
  #   ii <- paste("a.score.", i, sep = "")
  #   f0 <- as.formula(paste(" . ~ . - (0 + ", ii, "| ID)"))
  #   m0 <- update(fullReml, f0)
  #   f.slope <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + a.score.", additive0.sim, " + (0 +", ii, " | ID)", 
  #                               sep = ""))
  #   m.slope <- lmer(f.slope, data = designMatrix.lmm)
  #   tests2[[i]] <- exactRLRT(m.slope, fullReml, m0)
  # }
  # multiTest1 <- sapply(tests2, function(x) {
  #   c(statistic = x$statistic[1],
  #     "p-value" = x$p[1])})
  # pvalues.bonf <- p.adjust(multiTest1[2,], "bonferroni")
  f.slope <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + a.score.", 
                              additive0.sim, 
                              " + (0 + z.sim.uni | ID.uni)", 
                              sep = ""))
  m.slope <- lmer(f.slope, data = designMatrix.lmm)
  f0 <- as.formula(" . ~ . - (0 + z.sim.uni | ID.uni)")
  m0 <- update(fullReml, f0)
  
  tests2 <- exactRLRT(m.slope, fullReml, m0)
  pvalues.bonf <- tests2$p[1]
  
  # fullReml <- lmer(Y ~ 1 + dummyX + z.sim + (0 + dummyX | ID) + (0 + z.sim | ID), data = designMatrix) ## + (1 | temp : ID)
  # m0 <- update(fullReml,  . ~ . - (0 + z.sim | ID))
  # m.slope <- update(fullReml, . ~ . -  (0 + dummyX | ID))
  # tests2 <- list()
  # tests2[[1]] <- exactRLRT(m.slope, fullReml, m0) # , nsim = 100000
  # 
  # 
  # tests2 <- sapply(tests2, function(x){c(statistic = x$statistic[1], "p-value" = x$p[1])})
  # pvalues.bonf <- p.adjust(tests2[2,], "bonferroni")
  
  ######################################################################
  
  return(list(realTau = r.sim, 
              #pvalue = tests1$p[1], 
              pvalues.bonf = pvalues.bonf)) 
  #tests1 = tests1, 
  #tests2 = tests2))
}

# Setup parallel
#cores <- detectCores()
cluster <- makeCluster(cores)
nRandCovariate <- 2 
r.sim <- b.var
clusterExport(cluster, c("nRandCovariate", "r.sim")) # casting the coefficient parameter on the random effects' covariance function
fileName <- paste("248_grp50-rep20-", nRandCovariate,".RData", sep = "") # Saving file's name

# run the simulation
loopIndex <- 1
power2.sim <- list()

node_results <- parLapply(cluster, 1:simRep, run_one_sample)
Power2 <- mean(sapply(node_results, function(x) {       
  c(overall.sens = (sum(x$pvalues.bonf < pvalue.true) > 0))}))
power2.sim[[loopIndex]] <- list(Power = Power2, realTau = r.sim)

save(power2.sim,file=fileName) # Auto Save
stopCluster(cluster)





