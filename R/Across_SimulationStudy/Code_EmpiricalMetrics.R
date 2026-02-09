################################################################################
################################################################################
##########                   Empirical metrics                        ##########
################################################################################
################################################################################

################################################################################
###############                      Table 3                     ###############
################################################################################
S.area <- c(101,302)
c.save <- seq(25,50, by = 25)
nsim <- 25
scenarios <- c("S1", "S2")#, "S3", "S4")
prior <- c("iCAR", "LCAR")#, "LjCAR")


mat.res <- matrix(NA, ncol = 4*length(prior)+2, nrow = length(S.area)*length(scenarios)*3+16)
for (A in 1:length(S.area)) {
  for (pr in 1:length(prior)){
    for (s in 1:length(scenarios)) {
      meanMSS.dis1 <- 0
      meanRMSS.dis1 <- 0
      meanMSS.dis2 <- 0
      meanRMSS.dis2 <- 0
      
      maxMSS.dis1 <- 0
      maxRMSS.dis1 <- 0
      maxMSS.dis2 <- 0
      maxRMSS.dis2 <- 0
      
      meanRS.dis1 <- c()
      meanRS.dis2 <- c()
      meanSP <- c()
      
      
      for (cs in 1:length(c.save)) {
        if(S.area[A]==101){
          eval(parse(text = paste0("load('../../Data/Data_SimulationStudy_",scenarios[s],"_100.Rdata')")))
          
          if(prior[pr]=="iCAR"){
            eval(parse(text = paste0("load('./SimulationStudy_iCAR/Results_SimulationStudy_100_",scenarios[s],"_",c.save[cs],".Rdata')")))
            res.T <- iCAR.res
          }
          else if(prior[pr]=="LCAR"){
            eval(parse(text = paste0("load('./SimulationStudy_LCAR/Results_SimulationStudy_100_",scenarios[s],"_",c.save[cs],".Rdata')")))
            res.T <- LCAR.res
          }
          else{
            eval(parse(text = paste0("load('./SimulationStudy_LjCAR/Results_SimulationStudy_100_",scenarios[s],"_",c.save[cs],".Rdata')")))
            res.T <- LjCAR.res
          }
        }
        else{
          eval(parse(text = paste0("load('../../Data/Data_SimulationStudy_",scenarios[s],"_300.Rdata')")))
          
          if(prior[pr]=="iCAR"){
            eval(parse(text = paste0("load('./SimulationStudy_iCAR/Results_SimulationStudy_300_",scenarios[s],"_",c.save[cs],".Rdata')")))
            res.T <- iCAR.res
          }
          else if(prior[pr]=="LCAR"){
            eval(parse(text = paste0("load('./SimulationStudy_LCAR/Results_SimulationStudy_300_",scenarios[s],"_",c.save[cs],".Rdata')")))
            res.T <- LCAR.res
          }
          else{
            eval(parse(text = paste0("load('./SimulationStudy_LjCAR/Results_SimulationStudy_300_",scenarios[s],"_",c.save[cs],".Rdata')")))
            res.T <- LjCAR.res
          }
        }
        
        for (ns in 1:nsim) {
          res <- res.T[[ns]]
          sim.n <- (25*(cs-1)+ns)
          mean.r1 <- sum(DataSIM$observed[which(DataSIM$sim==sim.n & DataSIM$disease==1)])/sum(DataSIM$population[which(DataSIM$sim==sim.n & DataSIM$disease==1)])*10^5
          mean.r2 <- sum(DataSIM$observed[which(DataSIM$sim==sim.n & DataSIM$disease==2)])/sum(DataSIM$population[which(DataSIM$sim==sim.n & DataSIM$disease==2)])*10^5
          
          mean.r <- sum(DataSIM$observed[which(DataSIM$sim==sim.n)])/sum(DataSIM$population[which(DataSIM$sim==sim.n)])*10^5
          
          ##disease1
          meanRMSS.dis1 <- meanRMSS.dis1 + res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area[A], "]"), "Mean"]*10^5
          maxRMSS.dis1 <- maxRMSS.dis1 + max(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area[A], "]"), "Mean"]*10^5)
          meanRS.dis1[25*(cs-1)+ns] <- sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area[A], "]"), "Mean"]*10^10)/sum((mean.r1 - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==1)])^2)
          
          ##disease2
          meanRMSS.dis2 <- meanRMSS.dis2 + res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area[A], "]"), "Mean"]*10^5
          maxRMSS.dis2 <- maxRMSS.dis2 + max(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area[A], "]"), "Mean"]*10^5)
          meanRS.dis2[25*(cs-1)+ns] <- sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area[A], "]"), "Mean"]*10^10)/sum((mean.r2 - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==2)])^2)
          
          ###SP
          meanSP[25*(cs-1)+ns] <- (sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area[A], "]"), "Mean"]*10^10) + sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area[A], "]"), "Mean"]*10^10))/(sum((mean.r1 - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==1)])^2)+sum((mean.r2 - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==2)])^2))
        }
      }
      
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+1,1+5*(pr-1)] <- round(sum(meanRMSS.dis1)/(length(c.save)*nsim),2)
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+1,2+5*(pr-1)] <- round(maxRMSS.dis1/(length(c.save)*nsim),2)
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+1,3+5*(pr-1)] <- round(sum(meanRS.dis1)/(length(c.save)*nsim),2)
      
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+2,1+5*(pr-1)] <- round(sum(meanRMSS.dis2)/(length(c.save)*nsim),2)
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+2,2+5*(pr-1)] <- round(maxRMSS.dis2/(length(c.save)*nsim),2)
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+2,3+5*(pr-1)] <- round(sum(meanRS.dis2)/(length(c.save)*nsim),2)
      
      
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+3,1+5*(pr-1)] <- round(sum(meanRMSS.dis1 + meanRMSS.dis2)/(length(c.save)*nsim),2)
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+3,3+5*(pr-1)] <- round(sum(meanRS.dis1 + meanRS.dis2)/(length(c.save)*nsim),2)
      mat.res[5*(s-1)+(5*length(scenarios)+1)*(A-1)+3,4+5*(pr-1)] <- round(sum(meanSP)/(length(c.save)*nsim),2)
    }
  }
}

names <- rep(c(NA, "Disease 1 (j=1)", "Disease 2 (j=2)", "Total", NA),length(scenarios))
names[1] <- "Scenario 1"
names[6] <- "Scenario 2"
names[11] <- "Scenario 3"
names[16] <- "Scenario 4"

mat.res <- rbind(rep(NA, 4*length(prior)+length(scenarios)-1), mat.res)
mat.res <- cbind(c(rep(NA,5*length(scenarios)), "G=300", rep(NA,5*length(scenarios)-1)), 
                 c(names, NA, names[-length(names)]),
                 mat.res)

mat.res <- rbind(rep(NA, 4*length(prior)+length(scenarios)+1), mat.res)

mat.res[1,1] <- "G=100"

xtable::xtable(mat.res, include.rownames = FALSE,
               digits = c(0, 0, rep(2,4), 0, rep(2,4), 0, rep(2,4), 0))
