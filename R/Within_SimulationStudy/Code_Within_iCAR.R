################################################################################
################################################################################
##########                     Simulation Study                       ##########
################################################################################
################################################################################
# rm(list=ls())
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
load("../../Data/Carto_Spain_100areas_3Years.Rdata")
carto <- Carto.areas

sf::sf_use_s2(FALSE)
carto.nb <- poly2nb(carto)
W.nb <- nb2mat(carto.nb, style="B")
nbInfo <- nb2WB(carto.nb)



hotspot <- c("S1") ##replace `S1` with `S2`, `S3` or `S4` for Scenarios 2, 3 and 4, respectively

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
  }
  
  # priors:
  beta[1:nDisease, 1:nReg] <-  M[1:nDisease,1:nDisease] %*% phi[1:nDisease, 1:nReg]
  
})


var1 <- c(0.0025, 0.04, 0.25)
var2 <- c(0.0025, 0.04, 0.25)
cor <- c(0, 0.7)
for (h in 1:length(hotspot)) {
  for (v1 in var1) {
    for (v2 in var2) {
      for(r in cor){
        load(paste0("../../Data/Data_SimulationStudy_",hotspot[h],"_100.Rdata"))
        DataSIM <- DataSIM
        S.area <- length(unique(DataSIM$code))
        N.dis <- length(unique(DataSIM$disease))
        
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
          Sigma.b <-  matrix(c(v1, sqrt(v1*v2)*r, sqrt(v1*v2)*r, v2), ncol = 2)
          eig.val<- eigen(Sigma.b, symmetric = TRUE)$values
          eig.vec <- eigen(Sigma.b, symmetric = TRUE)$vectors
          M <- eig.vec %*% diag(sqrt(eig.val))
          
          constants <- list(pop = data$population[1:S.area],
                            rate = t(matrix(data$crude.rate/10^5, ncol = 2)),
                            nDisease = N.dis,
                            nReg = S.area, L = length(nbInfo$adj), num = nbInfo$num,
                            weights = nbInfo$weights, adj = nbInfo$adj,
                            M = M)
          
          inits = function(){
            list(mu = rep(0, N.dis), phi = rbind(rnorm(S.area, sd = 0.1),
                                                 rnorm(S.area, sd = 0.1)))}
          
          data.nimble <- list(O = t(matrix(data$observed, ncol = 2)))
          
          
          mcmc.out <- nimbleMCMC(code = code,
                                 constants = constants,
                                 data = data.nimble,
                                 inits = inits,
                                 nchains = 3,
                                 niter = 15000,
                                 nburnin = 5000,
                                 thin = 25,
                                 summary = TRUE,
                                 samples = TRUE,
                                 monitors = c('mu', 'r', 'MSS.r', 'beta',
                                              'RMSS.r'),
                                 samplesAsCodaMCMC = TRUE,
                                 setSeed = c(280320251, 29032025,30032025),
                                 WAIC = TRUE)
          
          iCAR.res[act.sim-c.save] <- list(mcmc.out)
          
          if (act.sim%%25 == 0){
            save(list = c("iCAR.res"),
                 file = paste0("./SimulationStudy_iCAR/Results_SimulationStudy_100_",hotspot[h],"_",v1,"_",v2,"_",r,"_",act.sim,".Rdata"))
            
            c.save<-c.save+25
            iCAR.res <- list()
          }
          
          act.sim <- act.sim +1
          
          if(act.sim==n.sim+1){break}
        }
      }
    }
  }
}


