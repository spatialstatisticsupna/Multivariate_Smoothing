################################################################################
################################################################################
##########                     Simulation Study                       ##########
################################################################################
################################################################################
rm(list=ls())
library(spdep)
library(MASS)
library(sf)

#################################################
## Create Folder for results  ##
#################################################
if(!dir.exists("SimulationStudy_iCAR")){
  dir.create("SimulationStudy_iCAR")
}

#################################################
## Load the cartography file  ##
#################################################
load("../../Data/Carto_Spain_300areas_3Years.Rdata")
carto <- Carto.areas

sf::sf_use_s2(FALSE)
carto.nb <- poly2nb(carto)
W.nb <- nb2mat(carto.nb, style="B")
nbInfo <- nb2WB(carto.nb)


hotspot <-  c("S1") ##replace `S1` with `S2`, `S3` or `S4` for Scenarios 2, 3 and 4, respectively

library(nimble)
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)



code <- nimbleCode({
  # likelihood:
  for(j in 1:nDisease){
    for(i in 1:nReg){
      O[j,i] ~ dpois(pop[i]*r[j,i])
      logit(r[j,i]) <- mu[j]+beta[j,i]
      
      #Metrics?
      MSS.r[j,i] <- (rate[j,i]-r[j,i])^2
      RMSS.r[j,i] <- (rate[j,i]-r[j,i])^2/r[j,i]
    
    }
    phi[j, 1:nReg] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nReg], 1, zero_mean = 1)
    
    
    mu[j] ~ dflat()
    
    
    A.d[j] ~ dchisq(nDisease+2-j+1)
  }
  for (j in 1:(nDisease-1)) {
    for(k in (j+1):nDisease){
      A.low[k,j] ~ dnorm(0, sd =1)
    }
  }
  A[1:nDisease,1:nDisease] <- A.low[1:nDisease,1:nDisease] + diag(sqrt(A.d[1:nDisease]))
  
  
  Sigma.b[1:nDisease,1:nDisease] <- A[1:nDisease,1:nDisease] %*% t(A[1:nDisease,1:nDisease])
  eig.val[1:nDisease] <- eigen(Sigma.b[1:nDisease,1:nDisease], symmetric = TRUE)$values
  eig.vec[1:nDisease,1:nDisease] <- eigen(Sigma.b[1:nDisease,1:nDisease], symmetric = TRUE)$vectors
  M[1:nDisease,1:nDisease] <- eig.vec[1:nDisease,1:nDisease] %*% diag(sqrt(eig.val[1:nDisease]))
  
  
  # priors:
  beta[1:nDisease, 1:nReg] <-  M[1:nDisease,1:nDisease] %*% phi[1:nDisease, 1:nReg]
})





for (h in 1:length(hotspot)) {
  load(paste0("../../Data/Data_SimulationStudy_",hotspot[h],"_300.Rdata"))
  S.area <- length(unique(DataSIM$code))
  N.dis <- length(unique(DataSIM$disease))#define
  
  n.sim <- 1000
  act.sim <- 1
  c.save <- 0
  
  
  iCAR.res <- list()
  repeat{
    print(act.sim)
    data <- DataSIM[which(DataSIM$sim==act.sim),]
    
    ##########################################################################
    ##########################################################################
    ####iCAR
    constants <- list(nDisease = N.dis,
                      nReg = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                      weights = nbInfo$weights, adj = nbInfo$adj, 
                      pop = data$population[1:S.area],
                      rate = t(matrix(data$crude.rate/10^5, ncol = 2)))
    
    
    inits <- function(){list(mu = rep(0, N.dis), phi = rbind(rnorm(S.area, sd = 0.1),
                                                             rnorm(S.area, sd = 0.1)),
                             A.d = c(rchisq(1, df=N.dis-1+3),rchisq(1, df=N.dis-2+3)),
                             A.low = matrix(c(0,rnorm(1, sd = 0.1),
                                              0, 0), ncol = N.dis))}

    data.nimble <- list(O = t(matrix(data$observed, ncol = 2)))
    
    
    
    mcmc.out <- nimbleMCMC(code = code,
                           constants = constants,
                           data = data.nimble,
                           inits = inits(),
                           nchains = 3,
                           niter = 30000,
                           nburnin = 5000,
                           thin = 75,
                           summary = TRUE,
                           samples = TRUE,
                           monitors = c('mu',  'Sigma.b', 'M',
                                        'r', 'MSS.r', 'beta',
                                        'RMSS.r'),
                           samplesAsCodaMCMC = TRUE,
                           setSeed = c(280320251, 29032025,30032025),
                           WAIC = TRUE)
    
    iCAR.res[act.sim-c.save] <- list(mcmc.out)
    
    if (act.sim%%25 == 0){
      save(list = c("iCAR.res"),
           file = paste0("./SimulationStudy_iCAR/Results_SimulationStudy_300_",hotspot[h],"_",act.sim,".Rdata"))
      
      c.save<-c.save+25
      iCAR.res <- list()
    }
    
    act.sim <- act.sim +1
    
    if(act.sim==n.sim+1){break}
  }
  
}
