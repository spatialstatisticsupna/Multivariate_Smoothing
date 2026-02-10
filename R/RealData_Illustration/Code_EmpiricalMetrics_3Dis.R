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
dis1 <- "Colon"
dis2 <- "Pancreas"
dis3 <- "Stomach"

for (sa in 1:length(N.area)) {
  load(paste0("../../Data/Carto_Spain_",N.area[sa],"areas_3Years.Rdata"))
  if (N.area[sa]==47) {carto <- Carto_prov}
  else  {carto <- Carto.areas}
  
  
  sf::sf_use_s2(FALSE)
  carto.nb <- poly2nb(carto)
  W.nb <- nb2mat(carto.nb, style="B")
  nbInfo <- nb2WB(carto.nb)
  nbInfo$num
  
  mat.res <- matrix(NA, nrow = 9*3+4, ncol = length(model.type)*4)
  eval(parse(text = paste0("load('./Results/Results_",N.area[sa],".Rdata')")))
  
  load(paste0("../../Data/Data_Spain_",N.area[sa],"areas_3Years.Rdata"))
  data <- Data.areas
  data$crude.rates <- data$O/data$Pop*10^5
  S.area <- length(unique(data$ID))
  
  
  mean.r1 <- sum(data$O[which(data$Type==dis1)])/sum(data$Pop[which(data$Type==dis1)])*10^5
  mean.r2 <- sum(data$O[which(data$Type==dis2)])/sum(data$Pop[which(data$Type==dis2)])*10^5
  mean.r3 <- sum(data$O[which(data$Type==dis3)])/sum(data$Pop[which(data$Type==dis3)])*10^5
  
  mean.r <- sum(data$O)/sum(data$Pop)*10^5
  
  for (mt in 1:length(model.type)) {
    eval(parse(text = paste0("res <- results$",model.type[mt],".res")))
    
    ##disease1
    meanRMSS.dis1 <- res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5
    maxRMSS.dis1 <- max(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    
    meanRS.dis1 <- sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r1 - data$crude.rates[which(data$Type==dis1)])^2)
    
    ##disease2
    meanRMSS.dis2 <- res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5
    maxRMSS.dis2 <- max(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5)
    
    meanRS.dis2 <- sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r2 - data$crude.rates[which(data$Type==dis2)])^2)
    
    ##disease3
    meanRMSS.dis3 <- res$summary$all.chains[paste0("RMSS.r[3, ", 1:S.area, "]"), "Mean"]*10^5
    maxRMSS.dis3 <- max(res$summary$all.chains[paste0("RMSS.r[3, ", 1:S.area, "]"), "Mean"]*10^5)
    
    meanRS.dis3 <- sum(res$summary$all.chains[paste0("MSS.r[3, ", 1:S.area, "]"), "Mean"]*10^10)/sum((mean.r3 - data$crude.rates[which(data$Type==dis3)])^2)
    
    ###Total SP
    meanSP <- (sum(res$summary$all.chains[paste0("MSS.r[1, ", 1:S.area, "]"), "Mean"]*10^10) + 
                 sum(res$summary$all.chains[paste0("MSS.r[2, ", 1:S.area, "]"), "Mean"]*10^10) + 
                 sum(res$summary$all.chains[paste0("MSS.r[3, ", 1:S.area, "]"), "Mean"]*10^10))/(sum((mean.r1 - data$crude.rates[which(data$Type==dis1)])^2)+
                                                                                                   sum((mean.r2 - data$crude.rates[which(data$Type==dis2)])^2)+
                                                                                                   sum((mean.r3 - data$crude.rates[which(data$Type==dis3)])^2))
    
    
    TCVi <- c()
    TCV <- matrix(NA, ncol = 3, nrow=333)
    cor1.it <- matrix(NA, ncol = 3, nrow=333)
    cor2.it <- matrix(NA, ncol = 3, nrow=333)
    cor3.it <- matrix(NA, ncol = 3, nrow=333)
    var.s1.it <- matrix(NA, ncol = 3, nrow=333)
    var.s2.it <- matrix(NA, ncol = 3, nrow=333)
    var.s3.it <- matrix(NA, ncol = 3, nrow=333)
    for (k in 1:3) {
      for (j in 1:333) {
        if(model.type[mt]=="LjCAR"){
          #For each iteration of each chain: inverse of number of neighbours of area i product Sigma_b inverse.
          eval(parse(text = paste0("mat.1 <- kronecker(diag(S.area), matrix(c(res$samples$chain",k,"[j,'M[1, 1]'],
					res$samples$chain",k,"[j,'M[2, 1]'], res$samples$chain",k,"[j,'M[3, 1]'],
					res$samples$chain",k,"[j,'M[1, 2]'], res$samples$chain",k,"[j,'M[2, 2]'],
					res$samples$chain",k,"[j,'M[3, 2]'],
					res$samples$chain",k,"[j,'M[1, 3]'], res$samples$chain",k,"[j,'M[2, 3]'],
					res$samples$chain",k,"[j,'M[3, 3]']), nrow = 3))")))
          eval(parse(text = paste0("mat.w1 <- inv(res$samples$chain",k,"[j,'rho[1]']*(diag(nbInfo$num) - W.nb) + (1-res$samples$chain",k,"[j,'rho[1]'])*diag(S.area))")))
          eval(parse(text = paste0("mat.w2 <- inv(res$samples$chain",k,"[j,'rho[2]']*(diag(nbInfo$num) - W.nb) + (1-res$samples$chain",k,"[j,'rho[2]'])*diag(S.area))")))
          eval(parse(text = paste0("mat.w3 <- inv(res$samples$chain",k,"[j,'rho[3]']*(diag(nbInfo$num) - W.nb) + (1-res$samples$chain",k,"[j,'rho[3]'])*diag(S.area))")))
          
          mat.2 <- matrix(0, ncol = S.area*3, nrow = S.area*3)
          for (i in 1:S.area) {
            mat.2[3*(i-1)+1, seq(1,S.area*3,3) ] <- mat.w1[i,]
            mat.2[3*(i-1)+2, seq(2,S.area*3,3) ] <- mat.w2[i,]
            mat.2[3*(i-1)+3, seq(3,S.area*3,3) ] <- mat.w3[i,]
          }
          
          
          eval(parse(text = paste0("mat.3 <- kronecker(diag(S.area), matrix(c(res$samples$chain",k,"[j,'M[1, 1]'],
					res$samples$chain",k,"[j,'M[2, 1]'], res$samples$chain",k,"[j,'M[3, 1]'],
					res$samples$chain",k,"[j,'M[1, 2]'], res$samples$chain",k,"[j,'M[2, 2]'],
					res$samples$chain",k,"[j,'M[3, 2]'],
					res$samples$chain",k,"[j,'M[1, 3]'], res$samples$chain",k,"[j,'M[2, 3]'],
					res$samples$chain",k,"[j,'M[3, 3]']), nrow = 3))")))
          
          Sigma <- tryCatch({ginv(mat.1%*%mat.2%*%mat.3)}, error = function(msg) {
            return(matrix(NA, ncol = S.area*3, nrow = S.area*3))})
          
        }
        for (i in 1:S.area) {
          #For each iteration of each chain: inverse of number of neighbours of area i product Sigma_b inverse.
          if (model.type[mt]=="iCAR") {
            #iCAR
            eval(parse(text = paste0("mat <- inv(nbInfo$num[i]*inv(matrix(c(res$samples$chain",k,"[j,'Sigma.b[1, 1]'],
						res$samples$chain",k,"[j,'Sigma.b[2, 1]'], res$samples$chain",k,"[j,'Sigma.b[3, 1]'],
						res$samples$chain",k,"[j,'Sigma.b[1, 2]'], res$samples$chain",k,"[j,'Sigma.b[2, 2]'], 
						res$samples$chain",k,"[j,'Sigma.b[3, 2]'],
						res$samples$chain",k,"[j,'Sigma.b[1, 3]'], res$samples$chain",k,"[j,'Sigma.b[2, 3]'],   
						res$samples$chain",k,"[j,'Sigma.b[3, 3]']), nrow = 3)))")))
          }
          else if(model.type[mt]=="LCAR"){
            #LCAR
            eval(parse(text = paste0("mat <- inv((res$samples$chain",k,"[j,'rho']*nbInfo$num[i] + (1-res$samples$chain",k,"[j,'rho']))*inv(matrix(c(res$samples$chain",k,"[j,'Sigma.b[1, 1]'],
						res$samples$chain",k,"[j,'Sigma.b[2, 1]'], res$samples$chain",k,"[j,'Sigma.b[3, 1]'],
						res$samples$chain",k,"[j,'Sigma.b[1, 2]'], res$samples$chain",k,"[j,'Sigma.b[2, 2]'], 
						res$samples$chain",k,"[j,'Sigma.b[3, 2]'],
						res$samples$chain",k,"[j,'Sigma.b[1, 3]'], res$samples$chain",k,"[j,'Sigma.b[2, 3]'],   
						res$samples$chain",k,"[j,'Sigma.b[3, 3]']), nrow = 3)))")))
          }
          else{
            mat <- tryCatch({inv(Sigma[(3*(i-1)+1:3),(3*(i-1)+1:3)])}, error = function(msg) {
              return(matrix(NA, ncol = 3, nrow = 3))})
          }
          #For each iteration of each chain: determinant for each area i
          TCVi[i] <- det(mat)
          
        }
        ##For each iteration of each chain: Sum of the determinants for all areas
        TCV[j,k] <- sum(TCVi)
        
        ##cor
        eval(parse(text =paste0("cor1.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[1, 2]']/sqrt(res$samples$chain",k,"[j,'Sigma.b[1, 1]']*res$samples$chain",k,"[j,'Sigma.b[2, 2]'])")))
        eval(parse(text =paste0("cor2.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[1, 3]']/sqrt(res$samples$chain",k,"[j,'Sigma.b[1, 1]']*res$samples$chain",k,"[j,'Sigma.b[3, 3]'])")))
        eval(parse(text =paste0("cor3.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[2, 3]']/sqrt(res$samples$chain",k,"[j,'Sigma.b[2, 2]']*res$samples$chain",k,"[j,'Sigma.b[3, 3]'])")))
        
        eval(parse(text =paste0("var.s1.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[1, 1]']")))
        eval(parse(text =paste0("var.s2.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[2, 2]']")))
        eval(parse(text =paste0("var.s3.it[j,k] <- res$samples$chain",k,"[j,'Sigma.b[3, 3]']")))
      }
    }
    
    mat.res[1+11*(na-1), 5*(mt-1)+1] <- round(sum(meanRMSS.dis1),2)
    mat.res[2+11*(na-1), 5*(mt-1)+1] <- round(maxRMSS.dis1,2)
    mat.res[3+11*(na-1), 5*(mt-1)+1] <- round(sum(meanRS.dis1),2)
    
    mat.res[1+11*(na-1), 5*(mt-1)+2] <- round(sum(meanRMSS.dis2),2)
    mat.res[2+11*(na-1), 5*(mt-1)+2] <- round(maxRMSS.dis2,2)
    mat.res[3+11*(na-1), 5*(mt-1)+2] <- round(sum(meanRS.dis2),2)
    
    mat.res[1+11*(na-1), 5*(mt-1)+3] <- round(sum(meanRMSS.dis3),2)
    mat.res[2+11*(na-1), 5*(mt-1)+3] <- round(maxRMSS.dis3,2)
    mat.res[3+11*(na-1), 5*(mt-1)+3] <- round(sum(meanRS.dis3),2)
    
    mat.res[1+11*(na-1), 5*(mt-1)+4] <- round(sum(meanRMSS.dis1 + meanRMSS.dis2 + meanRMSS.dis3),2)
    mat.res[4+11*(na-1), 5*(mt-1)+4] <- round(sum(meanSP),2)
    
    mat.res[5+11*(na-1), 4*(mt-1)+4] <- round(mean(as.vector(TCV)),3)
    mat.res[6+11*(na-1), 4*(mt-1)+1] <- round(mean(as.vector(var.s1.it)),3)
    mat.res[6+11*(na-1), 4*(mt-1)+2] <- round(mean(as.vector(var.s2.it)),3)
    mat.res[6+11*(na-1), 4*(mt-1)+3] <- round(mean(as.vector(var.s3.it)),3)
    mat.res[7+11*(na-1), 4*(mt-1)+4]<- round(mean(as.vector(cor1.it)),3)
    mat.res[8+11*(na-1), 4*(mt-1)+4]<- round(mean(as.vector(cor3.it)),3)
    mat.res[9+11*(na-1), 4*(mt-1)+4]<- round(mean(as.vector(cor2.it)),3)
  }
}
mat.res <- cbind(rep(c("RMSS", "maxRMSS", "RSP", "SP", "MultiTCV", "Sigma[jj]", 
                       "rho[12]", "rho[23]", "rho[13]", NA, NA),length(N.area)),
                 mat.res)
mat.res <- rbind(rep(NA, length(model.type)*4), rep(NA, length(model.type)*4), mat.res)

mat.res <- cbind(c("G=47", rep(NA,10), "G=100", rep(NA,10), "G=300", rep(NA,9)),
                 mat.res)

xtable::xtable(mat.res, include.rownames = FALSE,
               digits = c(0, 0, 3,3,3,3,0,3,3,3,3,0,3,3,3,3,0,0))




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
    dis1 <- "Colon"
    dis2 <- "Pancreas" 
    dis3 <- "Stomach"
    
    eval(parse(text = paste0("load('./Results/Results_",N.area[na],".Rdata')")))
    
    eval(parse(text = paste0("res <- results$",model.type[mt],".res")))
    dis1 <- as.numeric(res$summary$all.chains[paste0("r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    dis2 <- as.numeric(res$summary$all.chains[paste0("r[2, ", 1:S.area, "]"), "Mean"]*10^5)
    dis3 <- as.numeric(res$summary$all.chains[paste0("r[3, ", 1:S.area, "]"), "Mean"]*10^5)
    
    
    RMSS1 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[1, ", 1:S.area, "]"), "Mean"]*10^5)
    RMSS2 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[2, ", 1:S.area, "]"), "Mean"]*10^5)
    RMSS3 <- as.numeric(res$summary$all.chains[paste0("RMSS.r[3, ", 1:S.area, "]"), "Mean"]*10^5)
    
    mat <- rbind(mat, 
                 cbind(cbind(rep(ID.area, 3), rep(c("Colon", "Pancreas", "Stomach"), each = S.area),
                             data$crude.rates),
                       c(dis1, dis2, dis3),
                       
                       c(RMSS1, RMSS2, RMSS3),
                       rep("Colon-Pancreas-Stomach", 3*S.area),
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



###################################### Figure 2
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
                                 "Colon", "Colon-Pancreas-Stomach", "LCAR")
map2_higher <- create_choropleth(Carto_prov, subset(datos_long, N.area == "47"), 
                                 "Pancreas", "Colon-Pancreas-Stomach", "LCAR")
map3_higher <- create_choropleth(Carto_prov, subset(datos_long, N.area == "47"), 
                                 "Stomach", "Colon-Pancreas-Stomach", "LCAR")
# Or combine all maps in a grid
all_maps <- (map1_higher | map2_higher | map3_higher) 
# Save the combined plot
ggsave("RMSS_maps_47_LCAR.png", all_maps, width = 15, height = 6.5, dpi = 300)




###################################### Figure SuppMAtt

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
      legend.text = element_text(angle=45, vjust = 0.5,size = 15, color='gray35'),
      legend.key = element_rect()) +
    # scale_fill_binned( type = 'viridis', breaks = break_points) +
    binned_scale(aesthetics = 'fill', scale_name = 'custom', 
                 palette = ggplot2:::pal_binned(scales::manual_pal(values = cols)),
                 # guide = 'bins',
                 breaks = break_points,
                 guide = guide_bins(
                   direction = "horizontal",
                   title = " "
                 ))+ 
    theme(
      # legend.title = element_text(0.5),
      legend.position = "bottom",
      legend.key.width  = unit(1.3, "cm"),
      # legend.key.height = unit(1, "cm"),
      legend.key = element_rect(colour = "#f5f5f2", linewidth = 0.3)
  )
}





# Create the map
load('../../Data/Carto_Spain_100areas_3Years.Rdata')
map1_higher <- create_choropleth(Carto.areas, subset(datos_long, N.area == "100"), 
                                 "Colon", "Colon-Pancreas-Stomach", "LCAR")
map2_higher <- create_choropleth(Carto.areas, subset(datos_long, N.area == "100"), 
                                 "Pancreas", "Colon-Pancreas-Stomach", "LCAR")
map3_higher <- create_choropleth(Carto.areas, subset(datos_long, N.area == "100"), 
                                 "Stomach", "Colon-Pancreas-Stomach", "LCAR")

load('../../Data/Carto_Spain_300areas_3Years.Rdata')
map1_lower <- create_choropleth(Carto.areas, subset(datos_long, N.area == "300"), 
                                "Colon", "Colon-Pancreas-Stomach", "LCAR")
map2_lower <- create_choropleth(Carto.areas, subset(datos_long, N.area == "300"), 
                                "Pancreas", "Colon-Pancreas-Stomach", "LCAR")
map3_lower <- create_choropleth(Carto.areas, subset(datos_long, N.area == "300"),
                                "Stomach", "Colon-Pancreas-Stomach", "LCAR")

# # Or combine all maps in a grid
all_maps <- (map1_higher | map2_higher | map3_higher) / 
  (map1_lower | map2_lower | map3_lower) 
# Save the combined plot
ggsave("RMSS_maps_100_300.png", all_maps, width = 15, height = 11, dpi = 300)

