###############################################################################
################################################################################
##########                   Empirical metrics                        ##########
################################################################################
################################################################################
rm(list=ls())
library(spdep)
library(MASS)
library(sf)
library(matlib)

################################################################################
###############                   Table SuppMAtt                 ###############
################################################################################
N.area <- c(47,100,300)
model.type <- c("iCAR", "LCAR", "LjCAR")
dis1 <- c("Colon", "Colon", "Pancreas")
dis2 <- c("Pancreas", "Stomach", "Stomach")


for (sa in 1:length(N.area)) {
  load(paste0("../../Data/Carto_Spain_",N.area[sa],"areas_3Years.Rdata"))
  if (N.area[sa]==47) {carto <- Carto_prov}
  else  {carto <- Carto.areas}
  
  
  sf::sf_use_s2(FALSE)
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  nbInfo$num
  
  mat.res <- matrix(NA, nrow = 7*length(dis1)+4, ncol = length(model.type)*4)
  for (d1 in 1:length(dis1)) {
    eval(parse(text = paste0("load('./Results/Results_",dis1[d1],"_",dis2[d1],"_",N.area[sa],".Rdata')")))
    
    load(paste0("../../Data/Data_Spain_",N.area[sa],"areas_3Years.Rdata"))
    data <- Data.areas[which(Data.areas$Type%in%c(dis1[d1], dis2[d1])),]
    data$crude.rates <- data$O/data$Pop*10^5
    S.area <- length(unique(data$ID))
    
    
    
    mean.r1 <- sum(data$O[which(data$Type==dis1[d1])])/sum(data$Pop[which(data$Type==dis1[d1])])*10^5
    mean.r2 <- sum(data$O[which(data$Type==dis2[d1])])/sum(data$Pop[which(data$Type==dis2[d1])])*10^5
    
    mean.r <- sum(data$O)/sum(data$Pop)*10^5
    
    for (mt in 1:length(model.type)) {
      eval(parse(text = paste0("res <- results$",model.type[mt],".res")))
      
      ##disease1
      meanRMSS.dis1 <- res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5
      maxRMSS.dis1 <- max(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5)
      
      meanRS.dis1 <- sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r1 - data$crude.rates[which(data$Type==dis1[d1])])^2)
      
      ##disease2
      meanRMSS.dis2 <- res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5
      maxRMSS.dis2 <- max(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5)
      
      meanRS.dis2<- sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r2 - data$crude.rates[which(data$Type==dis2[d1])])^2)
      
      ###Total SP
      meanSP <- (sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area, "]"), "Mean"]*10^10) + sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area, "]"), "Mean"]*10^10))/(sum((mean.r1 - data$crude.rates[which(data$Type==dis1[d1])])^2)+sum((mean.r2 - data$crude.rates[which(data$Type==dis2[d1])])^2))
      
      
      TCVi <- c()
      TCV <- matrix(NA, ncol = 3, nrow=333)
      cor.it <- matrix(NA, ncol = 3, nrow=333)
      var.s1.it <- matrix(NA, ncol = 3, nrow=333)
      var.s2.it <- matrix(NA, ncol = 3, nrow=333)
      for (k in 1:3) {
        for (j in 1:333) {
          
          if(model.type[mt]=="LjCAR"){
            #For each iteration of each chain: inverse of number of neighbours of area i product Sigma_b inverse.
            eval(parse(text = paste0("mat.1 <- kronecker(diag(,S.area), matrix(c(res$samples$chain",k,"[j,'M[1, 1]'],
						res$samples$chain",k,"[j,'M[2, 1]'],
						res$samples$chain",k,"[j,'M[1, 2]'],
						res$samples$chain",k,"[j,'M[2, 2]']), nrow = 2))")))
            eval(parse(text = paste0("mat.w1 <- inv(res$samples$chain",k,"[j,'rho[1]']*(diag(nbInfo$num) - W.nb) + (1-res$samples$chain",k,"[j,'rho[1]'])*diag(S.area))")))
            eval(parse(text = paste0("mat.w2 <- inv(res$samples$chain",k,"[j,'rho[2]']*(diag(nbInfo$num) - W.nb) + (1-res$samples$chain",k,"[j,'rho[2]'])*diag(S.area))")))
            
            mat.2 <- matrix(0, ncol = S.area*2, nrow = S.area*2)
            for (i in 1:S.area) {
              mat.2[2*(i-1)+1, seq(1,S.area*2,2) ] <- mat.w1[i,]
              mat.2[2*(i-1)+2, seq(2,S.area*2,2) ] <- mat.w2[i,]
            }
            
            
            eval(parse(text = paste0("mat.3 <- kronecker(diag(S.area), matrix(c(res$samples$chain",k,"[j,'M[1, 1]'],
						res$samples$chain",k,"[j,'M[1, 2]'],
						res$samples$chain",k,"[j,'M[2, 1]'],
						res$samples$chain",k,"[j,'M[2, 2]']), nrow = 2))")))
            
            Sigma <- tryCatch({ginv(mat.1%*%mat.2%*%mat.3)}, error = function(msg) {
              return(matrix(NA, ncol = S.area*2, nrow = S.area*2))})
            
          }
          for (i in 1:S.area) {
            #For each iteration of each chain: inverse of number of neighbours of area i product Sigma_b inverse.
            if (model.type[mt]=="iCAR") {
              #iCAR
              eval(parse(text = paste0("mat <- inv(nbInfo$num[i]*inv(matrix(c(res$samples$chain",k,"[j,'Sigma.b[1, 1]'],
							res$samples$chain",k,"[j,'Sigma.b[2, 1]'],
							res$samples$chain",k,"[j,'Sigma.b[1, 2]'],
							res$samples$chain",k,"[j,'Sigma.b[2, 2]']), nrow = 2)))")))
            }
            else if(model.type[mt]=="LCAR"){
              #LCAR
              eval(parse(text = paste0("mat <- inv((res$samples$chain",k,"[j,'rho']*nbInfo$num[i] + (1-res$samples$chain",k,"[j,'rho']))*inv(matrix(c(res$samples$chain",k,"[j,'Sigma.b[1, 1]'],
							res$samples$chain",k,"[j,'Sigma.b[2, 1]'],
							res$samples$chain",k,"[j,'Sigma.b[1, 2]'],
							res$samples$chain",k,"[j,'Sigma.b[2, 2]']), nrow = 2)))")))
            }
            else{
              mat <- tryCatch({inv(Sigma[(2*(i-1)+1:2),(2*(i-1)+1:2)])}, error = function(msg) {
                return(matrix(NA, ncol = 2, nrow = 2))})
            }
            
            #For each iteration of each chain: determinant for each area i
            TCVi[i] <- det(mat)
            
          }
          ##For each iteration of each chain: Sum of the determinants for all areas
          TCV[j,k] <- sum(TCVi)
          
          
          ##cor
          eval(parse(text =paste0("cor.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[1, 2]']/sqrt(res$samples$chain",k,"[j,'Sigma.b[1, 1]']*res$samples$chain",k,"[j,'Sigma.b[2, 2]'])")))
          
          eval(parse(text =paste0("var.s1.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[1, 1]']")))
          eval(parse(text =paste0("var.s2.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[2, 2]']")))
        }
      }
      
      mat.res[1+9*(d1-1), 4*(mt-1)+1] <- round(sum(meanRMSS.dis1),2)
      mat.res[2+9*(d1-1), 4*(mt-1)+1] <- round(maxRMSS.dis1,2)
      mat.res[3+9*(d1-1), 4*(mt-1)+1] <- round(sum(meanRS.dis1),2)
      
      mat.res[1+9*(d1-1), 4*(mt-1)+2] <- round(sum(meanRMSS.dis2),2)
      mat.res[2+9*(d1-1), 4*(mt-1)+2] <- round(maxRMSS.dis2,2)
      mat.res[3+9*(d1-1), 4*(mt-1)+2] <- round(sum(meanRS.dis2),2)
      
      
      mat.res[1+9*(d1-1), 4*(mt-1)+3] <- round(sum(meanRMSS.dis1 + meanRMSS.dis2),2)
      # mat.res[2, 4*(mt-1)+3] <- round(sum(meanRS.dis1 + meanRS.dis2),2)
      mat.res[4+9*(d1-1), 4*(mt-1)+3] <- round(sum(meanSP),2)
      
      mat.res[5+9*(d1-1), 4*(mt-1)+3] <- round(mean(as.vector(TCV)),3)
      mat.res[6+9*(d1-1), 4*(mt-1)+1] <- round(mean(as.vector(var.s1.it)),3)
      mat.res[6+9*(d1-1), 4*(mt-1)+2] <- round(mean(as.vector(var.s2.it)),3)
      mat.res[7+9*(d1-1), 4*(mt-1)+3]<- round(mean(as.vector(cor.it)),3)
    }
  }
  mat.res <- cbind(rep(c("RMSS", "maxRMSS", "RSP", "SP", "MultiTCV", "Sigma[jj]", 
                         "rho", NA, NA),length(d1)),
                   mat.res)
  mat.res <- rbind(rep(NA, length(model.type)*4), rep(NA, length(model.type)*4), mat.res)
  
  mat.res <- cbind(c(paste("G=",N.area[sa]), paste(dis1[1]," - ",dis2[1]),
                     rep(NA,8), paste(dis1[2]," - ",dis2[2]), rep(NA,8),
                     paste(dis1[3]," - ",dis2[3]),rep(NA,7)),
                   mat.res)
  
  xtable::xtable(mat.res, include.rownames = FALSE,
                 digits = c(0, 0, 3,3,3,0,3,3,3,0,3,3,3,0,0))
}


################################################################################
###############                    Fig SuppMAtt                  ###############
################################################################################

# Create dataframe for plotting
N.area <- c(47,100,300)
model.type <- c("iCAR", "LCAR", "LjCAR")

mat <- NULL
for (na in 1:length(N.area)) {
  load(paste0("../../Data/Data_Spain_",N.area[na],"areas_3Years.Rdata"))
  data <- Data.areas
  data$crude.rates <- data$O/data$Pop*10^5
  ID.area <- unique(data$ID)
  S.area <- length(ID.area)

  
  for(mt in 1:length(model.type)){
    #####################################
    dis1 <- "Colon"
    dis2 <- "Stomach"
    eval(parse(text = paste0("load('./Results/Results_",dis1,"_",dis2,"_",N.area[na],".Rdata')")))

    eval(parse(text = paste0("res <- results$",model.type[mt],".res")))

    dis1 <- as.numeric(res$summary$all.chains[paste0("r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    dis2 <- rep(NA,S.area)
    dis3 <- as.numeric(res$summary$all.chains[paste0("r[2, ", 1:S.area, "]"), "Mean"]*10^5)


    RMSS1 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    RMSS2 <- rep(NA,S.area)
    RMSS3 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5)

    mat <- rbind(mat,
                 cbind(cbind(rep(ID.area, 3), rep(c("Colon", "Pancreas", "Stomach"), each = S.area),
                             data$crude.rates),
                       c(dis1, dis2, dis3),
                       c(RMSS1, RMSS2, RMSS3),
                       rep("Colon-Stomach", 3*S.area),
                       rep(model.type[mt], 3*S.area),
                       rep(N.area[na], 3*S.area)))

    
    #####################################
    dis1 <- "Colon" 
    dis2 <- "Pancreas"
    eval(parse(text = paste0("load('./Results/Results_",dis1,"_",dis2,"_",N.area[na],".Rdata')")))
    
    eval(parse(text = paste0("res <- results$",model.type[mt],".res")))
    
    dis1 <- as.numeric(res$summary$all.chains[paste0("r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    dis2 <- as.numeric(res$summary$all.chains[paste0("r[2, ", 1:S.area, "]"), "Mean"]*10^5)
    dis3 <- rep(NA,S.area)
    
    RMSS1 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    RMSS3 <- rep(NA,S.area)
    RMSS2 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5)
    
    mat <- rbind(mat, 
                 cbind(cbind(rep(ID.area, 3), rep(c("Colon", "Pancreas", "Stomach"), each = S.area),
                             data$crude.rates),
                       c(dis1, dis2, dis3),
                       c(RMSS1, RMSS2, RMSS3),
                       rep("Colon-Pancreas", 3*S.area),
                       rep(model.type[mt], 3*S.area),
                       rep(N.area[na], 3*S.area)))
    
    #####################################
    dis1 <- "Pancreas"
    dis2 <- "Stomach"
    eval(parse(text = paste0("load('./Results/Results_",dis1,"_",dis2,"_",N.area[na],".Rdata')")))

    eval(parse(text = paste0("res <- results$",model.type[mt],".res")))

    dis2 <- as.numeric(res$summary$all.chains[paste0("r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    dis3 <- as.numeric(res$summary$all.chains[paste0("r[2, ", 1:S.area, "]"), "Mean"]*10^5)
    dis1 <- rep(NA,S.area)

    RMSS2 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    RMSS1 <- rep(NA,S.area)
    RMSS3 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5)

    mat <- rbind(mat,
                 cbind(cbind(rep(ID.area, 3), rep(c("Colon", "Pancreas", "Stomach"), each = S.area),
                             data$crude.rates),
                       c(dis1, dis2, dis3),
                       c(RMSS1, RMSS2, RMSS3),
                       rep("Pancreas-Stomach", 3*S.area),
                       rep(model.type[mt], 3*S.area),
                       rep(N.area[na], 3*S.area)))
  }
  
}

datos_long <- as.data.frame(mat)
colnames(datos_long) <- c("ID.area", "Disease", "Observed", "Predicted", 
                          "RMSS", "Model", "Spatial", "N.area")

datos_long$Observed <- as.numeric(datos_long$Observed)
datos_long$Predicted <- as.numeric(datos_long$Predicted)
datos_long$RMSS <- as.numeric(datos_long$RMSS)



# packages for plotting
library(ggplot2)
library(sf)
library(dplyr)
library(patchwork)


# Function to create a choropleth map
create_choropleth <- function(carto, data, dis_name, model, spatial) {
  
  s.data <- subset(data, Disease == dis_name &
                     Model == model &
                     Spatial == spatial)
  break_points <- as.vector(round(quantile(s.data$RMSS, 
                                           probs = c(0.5,0.75,0.85,0.9,0.95,0.99)),2))

  cols <- c("#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594")
  colnames(s.data)[1] <- c("ID")
  map_data <- left_join(carto, s.data)
  
  
  ggplot(map_data)+
    geom_sf(aes(fill = RMSS), 
            color = "white", size = 0.1) +
    theme_void()+
    labs(title = dis_name,
         x = NULL,
         y = NULL) +
    theme(
      plot.title = element_text(size = 20, face = 'bold',
                                hjust=0.5,color='gray35'),
      plot.subtitle = element_text(size = 20,color='gray35'),
      # plot.caption = element_text(size=10,color='gray35'),
      strip.text.x = element_text(size=25,hjust=0.1,vjust=0, face = 'bold',
                                  color='gray35', margin = margin(0, 0, 20, 0, 'pt')),
      plot.background = element_rect(fill = '#f5f5f2', color = NA), 
      # plot.margin = margin(0.8, 0.5, 0.5, 0.5, 'cm'),
      panel.background = element_rect(fill = '#f5f5f2', color = NA),
      panel.border = element_blank(),
      legend.position = 'bottom',
      legend.background = element_rect(fill = '#f5f5f2', color = NA),
      legend.title = element_text(size = 20, color='gray35'),
      legend.text = element_text(angle=80, vjust = 0.5,size = 15, color='gray35'),
      legend.key = element_rect()) +
    # scale_fill_binned( type = 'viridis', breaks = break_points) +
    binned_scale(aesthetics = 'fill', scale_name = 'custom', 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols)),
                 guide = 'bins',
                 breaks = break_points)+
    guides(fill = guide_colourbar(direction = 'horizontal',  ## transform legend
                                  title=' ',  ##rename default legend
                                  title.position='top',
                                  title.hjust=0.5,
                                  ticks.colour='#f5f5f2',
                                  ticks.linewidth=3,
                                  barwidth = 19,
                                  barheight = 1))
}

# Create the map
load('../../Data/Carto_Spain_47areas_3Years.Rdata')
map1_higher <- create_choropleth(Carto_prov, subset(datos_long, N.area == "47"),
                                 "Colon", "Colon-Pancreas", "LCAR")
map2_higher <- create_choropleth(Carto_prov, subset(datos_long, N.area == "47"),
                                 "Pancreas", "Colon-Pancreas", "LCAR")
# Or combine all maps in a grid
all_maps <- (map1_higher | map2_higher) 
# Save the combined plot
ggsave("RMSS_maps_47_LCAR_Colon-Pancreas.png", all_maps, width = 11.5, height = 6.5, dpi = 300)


load('../../Data/Carto_Spain_100areas_3Years.Rdata')
map1_higher <- create_choropleth(Carto.areas, subset(datos_long, N.area == "100"),
                                 "Colon", "Colon-Pancreas", "LCAR")
map2_higher <- create_choropleth(Carto.areas, subset(datos_long, N.area == "100"),
                                 "Pancreas", "Colon-Pancreas", "LCAR")
# Or combine all maps in a grid
all_maps <- (map1_higher | map2_higher) 
# Save the combined plot
ggsave("RMSS_maps_100_LCAR_Colon-Pancreas.png", all_maps, width = 11.5, height = 6.5, dpi = 300)


load('../../Data/Carto_Spain_300areas_3Years.Rdata')
map1_higher <- create_choropleth(Carto.areas, subset(datos_long, N.area == "300"), 
                                "Colon", "Colon-Pancreas", "LCAR")
map2_higher <- create_choropleth(Carto.areas, subset(datos_long, N.area == "300"),
                                "Pancreas", "Colon-Pancreas", "LCAR")
# Or combine all maps in a grid
all_maps <- (map1_higher | map2_higher) 
# Save the combined plot
ggsave("RMSS_maps_300_LCAR_Colon-Pancreas.png", all_maps, width = 11.5, height = 6.5, dpi = 300)


