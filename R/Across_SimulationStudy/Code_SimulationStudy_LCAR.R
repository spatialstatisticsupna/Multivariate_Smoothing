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
if(!dir.exists("SimulationStudy_LCAR")){
  dir.create("SimulationStudy_LCAR")
}

#################################################
## Load the cartography file  ##
#################################################
load("../../Data/Carto_Spain_300areas_3Years.Rdata")
carto <- Carto.areas
S.area <- length(unique(carto$ID))

sf::sf_use_s2(FALSE)
cv.nb <- poly2nb(carto)
W.nb <- nb2mat(cv.nb, style="B")
nbInfo <- nb2WB(cv.nb)
# Number of neighbors of each municipality
nadj <- card(cv.nb)
# Neighbors of each municipality
map <- unlist(cv.nb)


hotspot <-  c("S1") ##replace `S1` with `S2`, `S3` or `S4` for Scenarios 2, 3 and 4, respectively

library(nimble)
nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)

##LCAR
# Diagonal matrix with the number of neighbors of each area
D <- diag(nadj)
# Adjacency matrix
W <- nb2mat(cv.nb, style = "B", zero.policy = TRUE)
# Eigenvalues of D-W
Lambda <- eigen(D - W)$values
# Identity matrix
I <- diag(rep(1, S.area))

# All the neighborhoods j ~ k where k < j
from.to <- cbind(rep(1:S.area, times = nadj), map); colnames(from.to) <- c("from", "to")
from.to <- from.to[which(from.to[, 1] < from.to[, 2]), ]
NDist <- nrow(from.to)


dcar_leroux <- nimbleFunction(
  name = 'dcar_leroux',
  run = function(x = double(1),        # Spatial random effect (vector)
                 rho = double(0),      # Amount of spatial dependence (scalar)
                 sd.theta = double(0), # Standard deviation (scalar)
                 Lambda = double(1),   # Eigenvalues of matrix D - W
                 from.to = double(2),  # Matrix of distinct pairs of neighbors from.to[, 1] < from.to[, 2]
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    # Number of small areas
    N <- dim(x)[1]
    # Number of distinct pairs of neighbors
    NDist <- dim(from.to)[1]
    # Required vectors
    x.from <- nimNumeric(NDist)
    x.to <- nimNumeric(NDist)
    for (Dist in 1:NDist) {
      x.from[Dist] <- x[from.to[Dist, 1]]
      x.to[Dist] <- x[from.to[Dist, 2]]
    }
    
    # Log-density
    logDens <- sum(dnorm(x[1:N], mean = 0, sd = sd.theta * pow(1 - rho, -1/2), log = TRUE)) -
      N/2 * log(1 - rho) + 1/2 * sum(log(rho * (Lambda[1:N] - 1) + 1)) -
      1/2 * pow(sd.theta, -2) * rho * sum(pow(x.from[1:NDist] - x.to[1:NDist], 2))
    if(log) return(logDens)
    else return(exp(logDens))
  }
)

# Another function of type rcar_leroux must be defined, otherwise it will cause an error
rcar_leroux <- nimbleFunction(
  name = 'rcar_leroux',
  run = function(n = integer(0),
                 rho = double(0),
                 sd.theta = double(0),
                 Lambda = double(1),
                 from.to = double(2)) {
    returnType(double(1))
    
    nimStop("user-defined distribution dcar_leroux provided without random generation function.")
    x <- nimNumeric(542)
    return(x)
  }
)

assign('dcar_leroux', dcar_leroux, envir = .GlobalEnv)
assign('rcar_leroux', rcar_leroux, envir = .GlobalEnv)


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
    phi[j, 1:nReg] ~ dcar_leroux(rho = rho,
                                 sd.theta = 1,
                                 Lambda = Lambda[1:nReg],
                                 from.to = from.to[1:NDist, 1:2])
    
    #Constraint
    # Zero-mean constraint for theta[1:NMuni]
    zero.theta[j] ~ dnorm(mean.thetas[j], 10000)
    mean.thetas[j] <- mean(phi[j, 1:nReg])
    
    
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
  
  # Hyperparameters of the spatial random effect
  rho ~ dunif(0, 1)

})





for (h in 1:length(hotspot)) {
  load(paste0("../../Data/Data_SimulationStudy_",hotspot[h],"_300.Rdata"))
  S.area <- length(unique(DataSIM$code))
  N.dis <- length(unique(DataSIM$disease))#define
  
  n.sim <- 1000
  act.sim <- 1
  c.save <- 0
  
  
  LCAR.res <- list()
  repeat{
    print(act.sim)
    data <- DataSIM[which(DataSIM$sim==act.sim),]
    
    ##########################################################################
    ##########################################################################
    ####L-CAR
    constants <- list(nDisease = N.dis,
                      nReg = S.area, NDist = NDist, Lambda = Lambda, 
                      from.to = from.to, 
                      pop = data$population[1:S.area],
                      rate = t(matrix(data$crude.rate/10^5, ncol = 2)))
    
    inits <- function(){list(mu = rep(0, N.dis), phi = rbind(rnorm(S.area, sd = 0.1),
                                                             rnorm(S.area, sd = 0.1)),
                             rho = runif(1),
                             A.d = c(rchisq(1, df=N.dis-1+3),rchisq(1, df=N.dis-2+3)),
                             A.low = matrix(c(0,rnorm(1, sd = 0.1),
                                              0, 0), ncol = N.dis))}
    
    data.nimble <- list(O = t(matrix(data$observed, ncol = 2)),
                        zero.theta = c(0,0))

    mcmc.out <- nimbleMCMC(code = code,
                           constants = constants,
                           data = data.nimble,
                           inits = inits(),
                           nchains = 3,
                           niter = 70000,
                           nburnin = 20000,
                           thin = 75,
                           summary = TRUE,
                           samples = TRUE,
                           monitors = c('mu',  'Sigma.b', 'M', 'rho',
                                        'r', 'MSS.r', 'beta',
                                        'RMSS.r'),
                           samplesAsCodaMCMC = TRUE,
                           setSeed = c(280320251, 29032025,30032025),
                           WAIC = TRUE)
    
    LCAR.res[act.sim-c.save] <- list(mcmc.out)
    
    if (act.sim%%25 == 0){
      save(list = c("LCAR.res"),
           file = paste0("./SimulationStudy_LCAR/Results_SimulationStudy_300_",hotspot[h],"_",act.sim,".Rdata"))
      
      c.save<-c.save+25
      LCAR.res <- list()
    }
    
    act.sim <- act.sim +1
    
    if(act.sim==n.sim+1){break}
  }
  
}
