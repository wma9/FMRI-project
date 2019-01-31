# data will be generated from tau1=tau^2, tau2=tau^2/2, and tau3=tau^2/4
# parameters are setted as estimations from real dataset ROI1
#load("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/03.13.2018/real data application/realdata.Rdata")

library(parallel)
simRep <- 5000 # Replication times in one simulation
pvalue.true <- .05 # Testing type I error 
b.var <- c(0) # The set of varaince of random covariates b as random slope
smooth <- 1 # measurement error is added to M if smooth = 0; no measurement error is added if sooth = 1
cores <- 4
r.sim <- b.var


run_one_sample <- function(iter){
  library(refund)
  library(lme4)
  library(nlme)
  library(arm)
  library(RLRsim)
  library(MASS)
  
  set.seed(iter)
  D <- 23 # grid number total
  nSubj <- 20 # 200 # I the number of curves
  nRep <-  20 # 20 # datasets for each covariance function 
  totalN <- nSubj * nRep
  thetaK.true <- c(2.88,-5.68,-13.27,12.56,-4.22,15.05,-16.23,-7.43)
  timeGrid <- (1:D)/D
  
  npc.true <- 8
  percent <- 0.9
  SNR <- 3 # 5, signal noise ratio'
  
  sd.epsilon <- 77.20 # or 0.5
  delta.true <- 216.960115
  a.mean <- 0
  
  gamma.true <- 155.773673
  #gammaVar.true <- 1
  # hot
  gammaI.true <- mapply(rnorm, nSubj, gamma.true, rep(47.4058, 1))
  gammaI.true <- gammaI.true[rep(1:nrow(gammaI.true), each = nRep), ]
  # warm
  gammaI2.true <- mapply(rnorm, nSubj, 0, rep(39.4067, 1))
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
    }else if(degree == 4){
      result <- sqrt(2) * cospi(2*t)
    }
    return(result)
  }
  
  #lambdaVec.true <- mapply(lambda.sim, 1: npc.true) 
  lambdaVec.true <- c(8.511901, 4.792678, 2.890541, 1.815417, 1.664546, 1.570437, 1.373613, 1.234652)
  psi.true <- matrix(NA, npc.true, D)
  psi.true[,1] <- c(-4.58229221,-0.09848766, -0.93119824, -0.13341140, -0.20641040,  0.02754230,  0.03061633,  0.17396017)
  psi.true[,2] <- c(-0.2160786, -0.3596236,  1.1626193, -0.4617492,1.3248845,-0.5763519,-1.4869616,-1.3572735)
  psi.true[,3] <- c(0.04127565, -0.70870226,  1.21536810, -1.30944320,  1.32049482, -0.28114728, -2.22240573, -0.60154694)
  psi.true[,4] <- c(-0.18511811, -0.88923674,  1.44468954,  0.40853337,  0.98619285, -0.05925722, -0.36151249, -0.22290398)
  psi.true[,5] <- c(-0.15076488, -0.99997205,  1.51649487,  1.83572942,  0.62929514,  0.05230834,  1.44807692,  0.16632868)
  psi.true[,6] <- c(0.0576725, -1.1156363,  1.0726455,  1.2783912,  0.2776352,  0.1594476,  1.6865626,  0.8444019)
  psi.true[,7] <- c(0.28911516, -1.20771405,  0.32835404, -0.47496477, -0.07678195,  0.27672427,  0.90311669,  1.34927307)
  psi.true[,8] <- c(0.4316456, -1.2200804, -0.3270418, -1.8199218, -0.4282464,  0.3814498,  0.2580890,  1.1535780)
  psi.true[,9] <- c(0.4627584, -1.1386244, -0.7118289, -1.8973583, -0.7294464,  0.4695396,  0.1942568,  0.2749822)
  psi.true[,10] <- c(0.4336590, -0.9935591, -0.8954491, -0.9046887, -0.9098747,  0.5441143,  0.3118642, -0.7435449)
  psi.true[,11] <- c(0.40912989, -0.81637640, -1.01407183,  0.43573809, -0.92540105,  0.55886466,  0.05500454, -1.27811810)
  psi.true[,12] <- c(0.4119738, -0.6007555, -1.0950459,  1.4340992, -0.8002524,  0.3769678, -0.6333096, -1.0515947)
  psi.true[,13] <- c(0.4180155, -0.3107292, -1.0522481,  1.7144792, -0.6069558, -0.1524760, -1.2625788, -0.2417027)
  psi.true[,14] <- c(0.39545261,  0.06970304, -0.83132382,1.25395810, -0.39188358, -1.01497471, -1.27828193,  0.69119483)
  psi.true[,15] <- c(0.3450467,  0.4891059, -0.5511043,  0.3761372, -0.1069349, -1.8698062, -0.5617180,  1.2133724)
  psi.true[,16] <- c(0.3029794,  0.8405075, -0.4631683, -0.3939246,  0.3530907, -2.1479496,  0.4652873,  0.9562520)
  psi.true[,17] <- c(0.30114471,  1.04097588, -0.70409060, -0.63264196,  0.99264356, -1.44497814,  1.18621789, -0.01188167)
  psi.true[,18] <- c(0.32406110,  1.11504735, -1.04059513, -0.32316152,  1.55546857,  0.06819829,  1.19480228, -0.95581343)
  psi.true[,19] <- c(0.3087857,  1.1918427, -0.9346253,  0.1192947,  1.5743835,  1.5982221,  0.5203247, -0.7996496)
  psi.true[,20] <- c(0.20864465,1.37411806, -0.08302559,  0.23731576,  0.73218161,  2.22173590, -0.41408531,  0.72213973)
  psi.true[,21] <- c(0.06565368,  1.58724820,  0.99215159, -0.01191975, -0.72750170,  1.49216694, -0.90581223,  1.86831939)
  psi.true[,22] <- c(-0.01247229,  1.57453886,  1.30346661, -0.24289162, -1.92397624, -0.05822514, -0.31123979,  0.07777254)
  psi.true[,23] <- c(-0.06028787,  1.17641014,  1.59902737, -0.48759945, -1.91260485, -0.62211566,  1.18368617, -2.22754523)
  
  #psi.true <- matrix(data = mapply(psi.fourier, rep(timeGrid, npc.true), rep(1:npc.true, each=D)),
  #                   nrow = npc.true, 
  #                   ncol = D, 
  #                   byrow = TRUE)
  ascore.true <- mvrnorm(totalN, rep(a.mean, npc.true), diag(lambdaVec.true))   
  Mt.true <-  ascore.true %*% psi.true     
  error <- rnorm(totalN, mean = 0, sd = sd.epsilon)
  
  #thetaIK.true <- mvrnorm(nSubj, rep(thetaK.true, npc.true), diag(rep(r.sim, npc.true)))
  #thetaIK.true <- mvrnorm(nSubj, rep(thetaK.true, npc.true), diag(c(r.sim, r.sim/2, r.sim/4, r.sim/8)))
  thetaIK.true <- mvrnorm(nSubj, thetaK.true, diag(rep(r.sim, npc.true)))
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
  knots <- 7 # previous setting 10
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
                             ascore = ascore,
                             ID.uni = as.factor(ID.uni),
                             z.sim.uni = z.sim.uni)
  
  
  # 'lmer' model
  designMatrix.lmm <- designMatrix
  
  additive0.sim <- paste(1:npc, collapse = " + ascore.")
  additive.sim <- paste(1:npc, collapse = " | ID) + (0 + ascore.")
  #additive.heter <- paste0(" + (0 + ascore.", 1:npc, " | ID)", collapse = "")
  #model.sim <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
  #                              additive0.sim, 
  #                              " + (0 + temp.1 | ID) + (0 + temp.2 | ID) ", 
  #                              additive.heter,
  #                              sep = ""))
  test_time1 <-system.time({
    
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
  
  # Confusion of modifying
  test_time2 <-system.time({
    
    model.sim <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                  additive0.sim, 
                                  " + (0 + temp.1 | ID) + (0 + temp.2 | ID) + (0 + z.sim.uni | ID.uni)", 
                                  sep = ""))
    
    f.slope <- as.formula(paste("rating ~ 1 + temp.1 + temp.2 + ascore.", 
                                additive0.sim, 
                                " + (0 + z.sim.uni | ID.uni)", 
                                sep = ""))
    m.slope <- lmer(f.slope, data = designMatrix.lmm)
    fullReml <- lmer(model.sim, data = designMatrix.lmm)
    f0 <- as.formula(" . ~ . - (0 + z.sim.uni | ID.uni)")
    m0 <- update(fullReml, f0)
    
    tests1 <- exactRLRT(m.slope, fullReml, m0)
    pvalues <- tests1$p[1] 
    
  }) 
  ###################################################################################
  
  
  return(list(realTau = r.sim,  
              pvalues.bonf = pvalues.bonf, 
              pvalues = pvalues,
              Merror.Var = Merror.Var,
              smooth = smooth,
              npc = npc, 
              tests.bonf = tests2,
              tests = tests1,
              test_time = c(test_time1[3],test_time2[3])))#test_time1 is for bonferroni; test_time2 is for our model
}

# Setup parallel
#cores <- detectCores()
cluster <- makeCluster(cores)

clusterExport(cluster, c("r.sim", "smooth")) # casting the coefficient parameter on the random effects' covariance function
fileName <- paste("test_8eigen_size_", smooth, "_",b.var,"_seed1_grp20-rep20.RData", sep = "") # Saving file's name

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

table1.sim <- sapply(node_results, function(x){       
  c(overall.sens = (sum(x$pvalues.bonf <= pvalue.true) > 0))})
Power1 <- mean(table1.sim)

table2.sim <- sapply(node_results, function(x){       
  c(overall.sens = (sum(x$pvalues <= pvalue.true) > 0))})
Power2 <- mean(table2.sim) 

npc.sim <- sapply(node_results, function(x){       
  x$npc})
npc.sum <- table(npc.sim) 

table2.time <- sapply(node_results, function(x){       
  x$test_time })
single_testtime <- rowMeans(table2.time) #row1 is testtime for bonferroni; row2 is time for our method

power2.sim[[loopIndex]] <- list(Power1.bonf = Power1, Power2 = Power2, 
                                realTau = c(r.sim,r.sim/2,r.sim/4), 
                                smooth = smooth, 
                                npc.table = npc.sum, 
                                single_test_sec = single_testtime)

save(power2.sim, file=fileName) # Auto Save

stopCluster(cluster)