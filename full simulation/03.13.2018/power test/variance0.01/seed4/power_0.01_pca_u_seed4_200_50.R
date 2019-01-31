library(parallel)
simRep <- 5000 # Replication times in one simulation
pvalue.true <- .05 # Testing type I error 
b.var <- c(0.01) # The set of varaince of random covariates b as random slope
smooth <- 0 # measurement error is added to M if smooth = 0; no measurement error is added if sooth = 1
cores <- 4
r.sim <- b.var


run_one_sample <- function(iter){
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  library(MASS)
  
  set.seed(iter+15000)
  D <- 80 # grid number total
  nSubj <- 200 # 200 # I the number of curves
  nRep <-  50 # 20 # datasets for each covariance function 
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
  
  thetaIK.true <- mvrnorm(nSubj, rep(thetaK.true, npc.true), diag(rep(r.sim, npc.true)))
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
  # to compare lambda: results$evalues/(D)) 
  # to compare estimated M, Mt.hat, Mt.true
  # a<-results$scores %*% t(results$efunctions)
  # plot(M[300,]) #Mt.hat
  # lines(a[300,]+results$mu,col="red") # estimated M
  # lines(Mt.true[300,], col="blue") #true Mt
  
  ###########################################################################
  dummyX <- cbind(dummyX, -dummyX + 1)
  z.sim.uni = c()
  ID.uni <- c(rbind(matrix(1:(nSubj*npc), 
                           nrow = npc, 
                           ncol = nSubj), 
                    matrix(0, nrow = nRep - npc, ncol = nSubj)))
  
  for(k in 1:nSubj){
    svd <- svd(ascore[((k-1)*nRep+1):(k*nRep), ] %*% t(ascore[((k-1)*nRep+1):(k*nRep), ])) #SVD on A_i
    u.tra <- t(svd$v)
    u <- svd$u
    d <- (svd$d)[1:npc]
    # u <- cbind(u, Null(u))
    Y[((k-1)*nRep+1):(k*nRep)] <- u.tra %*% Y[((k-1)*nRep+1):(k*nRep)]
    dummyX[((k-1)*nRep+1):(k*nRep), ] <- u.tra %*% dummyX[((k-1)*nRep+1):(k*nRep), ]
    ascore[((k-1)*nRep+1):(k*nRep), ] <- rbind(u.tra[1:npc, ] %*% ascore[((k-1)*nRep+1):(k*nRep), ], 
                                               matrix(0, 
                                                      nrow = nRep - npc, 
                                                      ncol = npc))
    z.sim.uni <- c(z.sim.uni, sqrt(d), rep(0, nRep - npc))
  }
  
  ###########################################################################
  designMatrix <- data.frame(rating = Y, 
                             temp.1 = dummyX[, 1],
                             temp.2 = dummyX[, 2],
                             ID = as.factor(ID),
                             ID.uni = as.factor(ID.uni),
                             ascore = ascore,
                             z.sim.uni = z.sim.uni)
  
  
  # 'lmer' model
  designMatrix.lmm <- designMatrix
  
  additive0.sim <- paste(1:npc, collapse = " + ascore.")
  additive.sim <- paste(1:npc, collapse = " | ID) + (0 + ascore.")
  # Confusion of modifying
  model.sim <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                additive0.sim, 
                                " + (0 + temp.1 | ID) + (0 + temp.2 | ID) + (0 + z.sim.uni | ID.uni)", 
                                sep = ""))
  fullReml <- lmer(model.sim, data = designMatrix.lmm)
  f.slope <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                              additive0.sim, 
                              " + (0 + z.sim.uni | ID.uni)", 
                              sep = ""))
  m.slope <- lmer(f.slope, data = designMatrix.lmm)
  f0 <- as.formula(" . ~ . - (0 + z.sim.uni | ID.uni)")
  m0 <- update(fullReml, f0)
  
  tests2 <- exactRLRT(m.slope, fullReml, m0)
  pvalues.bonf <- tests2$p[1] 
  
  ###################################################################################
  
  
  return(list(realTau = r.sim,  
              pvalues.bonf = pvalues.bonf, 
              Merror.Var = Merror.Var,
              smooth = smooth,
              npc = npc, 
              tests2 = tests2))
}

# Setup parallel
#cores <- detectCores()
cluster <- makeCluster(cores)
clusterSetRNGStream(cluster, 20170822)

# for(nRandCovariate in 1 * 2){  # START out-outer loop
# clusterExport(cluster, c("nRandCovariate")) # casting the coefficient parameter on the random effects' covariance function
# fileName <- paste("power_", b.var, "_grp20-rep20-", nRandCovariate,".RData", sep = "") # Saving file's name

clusterExport(cluster, c("r.sim", "smooth")) # casting the coefficient parameter on the random effects' covariance function
fileName <- paste("f_power_", smooth, "_",b.var,"_seed4_grp200-rep50.RData", sep = "") # Saving file's name

# run the simulation
loopIndex <- 1
# resultDoubleList.sim <- list()
#power1.sim <- list()
power2.sim <- list()
#  for(r.sim in b.var){  # START outer loop

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

table2.sim <- sapply(node_results, function(x) {       
  c(overall.sens = (sum(x$pvalues.bonf <= pvalue.true) > 0))})
Power2 <- mean(table2.sim)
#cat("Power2: ", Power2, fill = TRUE)
power2.sim[[loopIndex]] <- list(Power = Power2, realTau = r.sim, smooth = smooth)

#  loopIndex <- loopIndex + 1
# } # End outer loop

save.image(file=fileName) # Auto Save


#  par(mfrow=c(2,1))
#  Histogram plots
#  hist(sapply(result1.sim, function(x) x$pvalue), 
#       main = "Histogram of p-value for lme model", 
#       xlab = "p-value")

#  hist(sapply(result2.sim, function(x) x$pvalues.bonf), 
#       main = "Histogram of p-value for lmer model", 
#       xlab = "p-value")

# hist(sapply(resultDoubleList.sim[[1]], function(x) (x$tests1)$statistic[1]),
#      breaks = (0:110)/10, 
#      main = "Histogram of test-statistic for lme model",
#      xlab = "Test Statistics")
# 
# hist(sapply(resultDoubleList.sim[[1]], function(x) (x$tests2)[1,1]), 
#      breaks = (0:100)/10, 
#      main = "Histogram of test-statistic for lmer model",
#      xlab = "Test Statistics")

#} # End out-outer loop

stopCluster(cluster)