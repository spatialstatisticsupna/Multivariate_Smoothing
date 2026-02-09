################################################################################
################################################################################
##########               Within Simulation Study                      ##########
################################################################################
################################################################################
rm(list=ls())
library(spdep)
# library(Matrix)
library(MASS)
library(sf)
library(matlib)

#################################################
## Load the cartography file  ##
#################################################
load("../../Data/Carto_Spain_100areas_3Years.Rdata")
carto <- Carto.areas

sf::sf_use_s2(FALSE)
carto.nb <- poly2nb(carto)
W.nb <- nb2mat(carto.nb, style="B")
nbInfo <- nb2WB(carto.nb)
nbInfo$num


################################################################################
################################################################################
###############                      TABLES                      ###############
################################################################################
################################################################################
################################################################################
###############                 Table 1 (iCAR)                   ###############
################################################################################
S.area <- 101
c.save <- seq(25, 1000, by = 25)
nsim <- 25
var1 <- c("0.0025", "0.04", "0.25")
var2 <- c("0.0025", "0.04", "0.25")
ro <- c(0, 0.7)

mat.res <- matrix(NA, ncol = 5*length(ro)+1, nrow = length(var1)*length(var2)*4)
for (s in 1:length(ro)) {
  for (v1 in 1:length(var1)) {
    for (v2 in 1:length(var2)) {
      meanRMSS.dis1 <- 0
      meanRMSS.dis2 <- 0
      
      maxRMSS.dis1 <- 0
      maxRMSS.dis2 <- 0
      
      
      meanRS.dis1 <- c()
      meanRS.dis2 <- c()
      meanSP <- c()
      
      
      #######Separable
      TCVi <- c()
      TCVi.tr <- c()
      for (i in 1:101) {
        #For each iteration of each chain: inverse of number of neighbours of area i product Sigma_b inverse.
        eval(parse(text = paste0("mat <- ginv(nbInfo$num[i]*ginv(matrix(c(as.numeric(var1[v1]),
				sqrt(as.numeric(var1[v1])*as.numeric(var2[v2]))*as.numeric(ro[s]),
				sqrt(as.numeric(var1[v1])*as.numeric(var2[v2]))*as.numeric(ro[s]),as.numeric(var2[v2])), nrow = 2)))")))
        
        #For each iteration of each chain: determinant for each area i
        TCVi[i] <- det(mat)
        TCVi.tr[i] <- tr(mat)
        
      }
      TCV <- sum(TCVi)
      TCV.tr <- sum(TCVi.tr)
      
      for (cs in 1:length(c.save)) {
        eval(parse(text = paste0("load('./SimulationStudy_iCAR/Results_SimulationStudy_100_S1_",var1[v1],"_",var2[v2],"_",ro[s],"_",c.save[cs],".Rdata')")))
        eval(parse(text = paste0("load('../../Data/Data_SimulationStudy_S1_100.Rdata')")))
        
        for (ns in 1:nsim) {
          res <- iCAR.res[[ns]]
          
          sim.n <- (25*(cs-1)+ns)
          mean.r1 <- sum(DataSIM$observed[which(DataSIM$sim==sim.n & DataSIM$disease==1)])/sum(DataSIM$population[which(DataSIM$sim==sim.n & DataSIM$disease==1)])*10^5
          mean.r2 <- sum(DataSIM$observed[which(DataSIM$sim==sim.n & DataSIM$disease==2)])/sum(DataSIM$population[which(DataSIM$sim==sim.n & DataSIM$disease==2)])*10^5
          mean.r <- sum(DataSIM$observed[which(DataSIM$sim==sim.n)])/sum(DataSIM$population[which(DataSIM$sim==sim.n)])*10^5
          
          ##disease1
          meanRMSS.dis1 <- meanRMSS.dis1 + res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5
          
          maxRMSS.dis1 <- maxRMSS.dis1 + max(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5)
          
          meanRS.dis1[25*(cs-1)+ns] <- sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r1 - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==1)])^2)
          
          
          ##disease2
          meanRMSS.dis2 <- meanRMSS.dis2 + res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5
          
          maxRMSS.dis2 <- maxRMSS.dis2 + max(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5)
          
          meanRS.dis2[25*(cs-1)+ns] <- sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r2 - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==2)])^2)
          
          ###Total SP
          meanSP[25*(cs-1)+ns] <- (sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area, "]"), "Mean"]*10^10) + sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area, "]"), "Mean"]*10^10))/(sum((mean.r - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==1)])^2)+sum((mean.r - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==2)])^2))
          
        }
      }
      
      
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+1,4+6*(s-1)] <- round(sum(meanRMSS.dis1)/1000,2)
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+1,5+6*(s-1)] <- round(maxRMSS.dis1/1000,2)
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+1,3+6*(s-1)] <- round(sum(meanRS.dis1)/1000,2)
      
      
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+2,4+6*(s-1)] <- round(sum(meanRMSS.dis2)/1000,2)
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+2,5+6*(s-1)] <- round(maxRMSS.dis2/1000,2)
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+2,3+6*(s-1)] <- round(sum(meanRS.dis2)/1000,2)
      
      
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+3,4+6*(s-1)] <- round(sum(meanRMSS.dis1 + meanRMSS.dis2)/1000,2)
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+3,3+6*(s-1)] <- round(sum(meanRS.dis1 + meanRS.dis2)/1000,2)
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+3,2+6*(s-1)] <- round(sum(meanSP)/1000,2)
      mat.res[4*(v2-1)+4*length(var2)*(v1-1)+3,1+6*(s-1)] <- round(TCV,4)
    }
  }
  
}


mat.res <- cbind(rep(c("0.0025", rep(NA,3), "0.04", rep(NA,3), "0.25", rep(NA,3)), 
                     length(var1)), 
                 rep(c("Disease 1", "Disease 2", "Total", NA),length(var1)*length(var2)*4),
                 mat.res)
mat.res <- rbind(rep(NA, 13), mat.res)

mat.res <- cbind(c("0.0025", rep(NA,length(var2)*4), "0.04", rep(NA,length(var2)*4),
                   "0.25", rep(NA,length(var2)*4)),
                 mat.res)

xtable::xtable(mat.res, include.rownames = FALSE,
               digits = c(0, 0, 0, 0, 4, 2, 2, 2, 2, 0, 4, 2, 2, 2, 2))


################################################################################
###############                Table 2 (LCAR)                    ###############
################################################################################
S.area <- 101
c.save <- seq(25,1000, by = 25)
nsim <- 25
var1 <- c("0.0025", "0.25")
var2 <- c("0.0025", "0.25")
lambda <- c(0.2, 0.8)
ro <- c(0, 0.7)

mat.res <- matrix(NA, ncol = 5*length(ro) + 1, nrow = length(var1)*length(var2)*length(lambda)*4)
for (s in 1:length(ro)) {
  for (v1 in 1:length(var1)) {
    for (v2 in 1:length(var2)) {
      for(l in 1:length(lambda)) {
        meanRMSS.dis1 <- 0
        meanRMSS.dis2 <- 0
        
        maxRMSS.dis1 <- 0
        maxRMSS.dis2 <- 0
        
        meanRS.dis1 <- c()
        meanRS.dis2 <- c()
        meanSP <- c()
        
        
        ######Separable
        TCVi <- c()
        TCVi.tr <- c()
        for (i in 1:101) {
          #For each iteration of each chain: inverse of number of neighbours of area i product Sigma_b inverse.
          eval(parse(text = paste0("mat <- ginv((lambda[l]*nbInfo$num[i] + (1-lambda[l]))*ginv(matrix(c(as.numeric(var1[v1]),
					sqrt(as.numeric(var1[v1])*as.numeric(var2[v2]))*as.numeric(ro[s]),
					sqrt(as.numeric(var1[v1])*as.numeric(var2[v2]))*as.numeric(ro[s]),as.numeric(var2[v2])), nrow = 2)))")))
          #For each iteration of each chain: determinant for each area i
          TCVi[i] <- det(mat)
          # TCVi.tr[i] <- tr(mat)
          
        }
        TCV <- sum(TCVi)
        # TCV.tr <- sum(TCVi.tr)
        
        
        
        for (cs in 1:length(c.save)) {
          eval(parse(text = paste0("load('./SimulationStudy_LCAR/Results_SimulationStudy_100_S1_",var1[v1],"_",var2[v2],"_",ro[s],"_",lambda[l],"_",c.save[cs],".Rdata')")))
          eval(parse(text = paste0("load('../../Data/Data_SimulationStudy_S1_100.Rdata')")))
          
          for (ns in 1:nsim) {
            res <- LCAR.res[[ns]]
            
            sim.n <- (25*(cs-1)+ns)
            mean.r1 <- sum(DataSIM$observed[which(DataSIM$sim==sim.n & DataSIM$disease==1)])/sum(DataSIM$population[which(DataSIM$sim==sim.n & DataSIM$disease==1)])*10^5
            mean.r2 <- sum(DataSIM$observed[which(DataSIM$sim==sim.n & DataSIM$disease==2)])/sum(DataSIM$population[which(DataSIM$sim==sim.n & DataSIM$disease==2)])*10^5
            mean.r <- sum(DataSIM$observed[which(DataSIM$sim==sim.n)])/sum(DataSIM$population[which(DataSIM$sim==sim.n)])*10^5
            
            ##disease1
            meanRMSS.dis1 <- meanRMSS.dis1 + res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5
            
            maxRMSS.dis1 <- maxRMSS.dis1 + max(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5)
            
            
            meanRS.dis1[25*(cs-1)+ns] <- sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r1 - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==1)])^2)
            
            
            ##disease2
            meanRMSS.dis2 <- meanRMSS.dis2 + res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5
            
            maxRMSS.dis2 <- maxRMSS.dis2 + max(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5)
            
            meanRS.dis2[25*(cs-1)+ns] <- sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r2 - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==2)])^2)
            
            ###Total SP
            meanSP[25*(cs-1)+ns] <- (sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area, "]"), "Mean"]*10^10) + sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area, "]"), "Mean"]*10^10))/(sum((mean.r - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==1)])^2)+sum((mean.r - DataSIM$crude.rate[which(DataSIM$sim==sim.n & DataSIM$disease==2)])^2))
            
          }
        }
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+1,4 + 6*(s-1)] <- round(sum(meanRMSS.dis1)/1000,2)
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+1,5 + 6*(s-1)] <- round(maxRMSS.dis1/1000,2)
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+1,3 + 6*(s-1)] <- round(sum(meanRS.dis1)/1000,2)
        
        
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+2,4 + 6*(s-1)] <- round(sum(meanRMSS.dis2)/1000,2)
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+2,5 + 6*(s-1)] <- round(maxRMSS.dis2/1000,2)
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+2,3 + 6*(s-1)] <- round(sum(meanRS.dis2)/1000,2)
        
        
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+3,4 + 6*(s-1)] <- round(sum(meanRMSS.dis1 + meanRMSS.dis2)/1000,2)
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+3,3 + 6*(s-1)] <- round(sum(meanRS.dis1 + meanRS.dis2)/1000,2)
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+3,2 + 6*(s-1)] <- round(sum(meanSP)/1000,2)
        mat.res[4*(l-1)+4*length(lambda)*(v2-1)+4*length(var2)*length(lambda)*(v1-1)+3,1 + 6*(s-1)] <- round(TCV,4)
        # mat.res[4*(v2-1)+12*(v1-1)+3,2] <- round(TCV.tr,4)
      }
      
    }
  }
  
}

mat.res <- cbind(rep(c("0.0025", rep(NA,7), "0.25", rep(NA,7)), 
                     length(var1)), 
                 rep(c("0.2", rep(NA,3), "0.8", rep(NA,3)), 
                     length(var1)*length(var2)),
                 rep(c("Disease 1", "Disease 2", "Total", NA),length(var1)*length(var2)*length(lambda)*4),
                 mat.res)
mat.res <- rbind(rep(NA, 13), mat.res)

mat.res <- cbind(c("0.0025", rep(NA,length(var2)*length(lambda)*4),
                   "0.25", rep(NA,length(var2)*length(lambda)*4)),
                 mat.res)

xtable::xtable(mat.res, include.rownames = FALSE,
               digits = c(0, 0, 0, 0, 4, 2, 2, 2, 2, 0, 4, 2, 2, 2, 2))




################################################################################
################################################################################
###############                     FIGURES                      ###############
################################################################################
################################################################################
################################################################################
################################################################################
###############                 Figure 1 (iCAR)                  ###############
################################################################################
S.area <- 101
c.save <- seq(25,1000, by = 25)
nsim <- 25
scenarios <- c("S1", "S2", "S3", "S4")

var1 <- c("0.04")
var2 <- c("0.0025", "0.04", "0.25")
ro <- c(0, 0.7)

mat <- NULL
for (sc in 1:length(scenarios)) {
  for (s in 1:length(ro)) {
    for (v1 in 1:length(var1)) {
      for (v2 in 1:length(var2)) {
        meanRMSS.dis1 <- 0
        meanRMSS.dis2 <- 0
        
        
        for (cs in 1:length(c.save)) {
          eval(parse(text = paste0("load('../../Data/Data_SimulationStudy_",scenarios[sc],"_100.Rdata')")))
          
          eval(parse(text = paste0("load('./SimulationStudy_iCAR/Results_SimulationStudy_100_",scenarios[sc],"_",var1[v1],"_",var2[v2],"_",ro[s],"_",c.save[cs],".Rdata')")))
        
          for (ns in 1:nsim) {
            res <- iCAR.res[[ns]]
            
            sim.n <- (25*(cs-1)+ns)
            
            ##disease1
            meanRMSS.dis1 <- meanRMSS.dis1 + res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5
            
            ##disease2
            meanRMSS.dis2 <- meanRMSS.dis2 + res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5
            
          }
        }
        
        mat <- rbind(mat,
                     cbind(scenarios[sc], ro[s], var1[v1], var2[v2],
                           round(sum(meanRMSS.dis1)/1000,2), 
                           round(sum(meanRMSS.dis2)/1000,2),
                           round(sum(meanRMSS.dis1 + meanRMSS.dis2)/1000,2)))
        
      }
    }
  }
}


datos <- data.frame(mat)
colnames(datos) <- c("Scenario", "rho", "Sigma11", "Sigma22", "RMSS1", "RMSS2", "Total")
datos$rho <- as.factor(datos$rho)
datos$RMSS1 <- as.numeric(datos$RMSS1)
datos$RMSS2 <- as.numeric(datos$RMSS2)
datos$Total <- as.numeric(datos$Total)


nuevos_nombres <- c(
  "0.0025" = "Sigma[22] == 0.0025",
  "0.04" = "Sigma[22] == 0.04",
  "0.25" = "Sigma[22] == 0.25"
)

# OPCIÓN 1: Gráfico de pendientes (Slope Chart)
p1 <-ggplot(datos, aes(x = rho, y = Total, group = Scenario, color = Scenario)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(labels = parse(text = c("rho == 0", "rho == 0.7")))+
  facet_wrap(~Sigma22,  nrow = 1,
             labeller = labeller(Sigma22 = as_labeller(nuevos_nombres, label_parsed))) +
  labs(
    title = expression(paste(Sigma[11], " = 0.04")),
    x = "Correlation",
    y = "Total RMSS",
    color = "Scenario "
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.title.x = element_text(size = 12),        
    axis.title.y = element_text(size = 12),    
    axis.text.x = element_text(size = 11),       
    axis.text.y = element_text(size = 11),      
    legend.title = element_text(size = 12),       
    legend.text = element_text(size = 11),     
    
  )


pdf("./Figure/Figure1.pdf", width = 10, height = 5)
p1
dev.off()


################################################################################
###############                 Figure 2 (LCAR)                  ###############
################################################################################
S.area <- 101
c.save <- seq(25,1000, by = 25)
nsim <- 25
scenarios <- c("S1", "S2", "S3", "S4")

var1 <- c("0.04")
var2 <- c("0.0025", "0.04", "0.25")
ro <- c(0)#, 0.7)
lambda <- c(0.2, 0.8)

mat <- NULL
for (sc in 1:length(scenarios)) {
  for (s in 1:length(lambda)) {
    for (v1 in 1:length(var1)) {
      for (v2 in 1:length(var2)) {
        meanRMSS.dis1 <- 0
        meanRMSS.dis2 <- 0
        
        
        for (cs in 1:length(c.save)) {
          eval(parse(text = paste0("load('../../Data/Data_SimulationStudy_",scenarios[sc],"_100.Rdata')")))
          eval(parse(text = paste0("load('./SimulationStudy_LCAR/Results_SimulationStudy_100_",scenarios[sc],"_",var1[v1],"_",var2[v2],"_",ro,"_",lambda[s],"_",c.save[cs],".Rdata')")))
          
          
          for (ns in 1:nsim) {
            res <- LCAR.res[[ns]]
            
            sim.n <- (25*(cs-1)+ns)
            
            ##disease1
            meanRMSS.dis1 <- meanRMSS.dis1 + res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5
            
            ##disease2
            meanRMSS.dis2 <- meanRMSS.dis2 + res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5
            
          }
        }
        
        mat <- rbind(mat,
                     cbind(scenarios[sc], lambda[s], var1[v1], var2[v2],
                           round(sum(meanRMSS.dis1)/1000,2), 
                           round(sum(meanRMSS.dis2)/1000,2),
                           round(sum(meanRMSS.dis1 + meanRMSS.dis2)/1000,2)))
        
      }
    }
  }
}


datos <- data.frame(mat)
colnames(datos) <- c("Scenario", "lambda", "Sigma11", "Sigma22", "RMSS1", "RMSS2", "Total")
datos$lambda <- as.factor(datos$lambda)
datos$RMSS1 <- as.numeric(datos$RMSS1)
datos$RMSS2 <- as.numeric(datos$RMSS2)
datos$Total <- as.numeric(datos$Total)


nuevos_nombres <- c(
  "0.0025" = "Sigma[22] == 0.0025",
  "0.04" = "Sigma[22] == 0.04",
  "0.25" = "Sigma[22] == 0.25"
)

# OPCIÓN 1: Gráfico de pendientes (Slope Chart)
p2 <-ggplot(datos, aes(x = lambda, y = Total, group = Scenario, color = Scenario)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(labels = parse(text = c("lambda == 0.2", "lambda == 0.8")))+
  facet_wrap(~Sigma22,  nrow = 1,
             labeller = labeller(Sigma22 = as_labeller(nuevos_nombres, label_parsed))) +
  labs(
    title = expression(paste(Sigma[11], " = 0.04")),
    x = "Correlation",
    y = "Total RMSS",
    color = "Scenario "
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.title.x = element_text(size = 12),       
    axis.title.y = element_text(size = 12),  
    axis.text.x = element_text(size = 11),       
    axis.text.y = element_text(size = 11),        
    legend.title = element_text(size = 12),        
    legend.text = element_text(size = 11),       
    
  )


pdf("./Figure/Figure2.pdf", width = 10, height = 5)
p2
dev.off()


