# Pipeline for the Quantification and plotting of lipid droplets along gut length and intensity along length
# written by Laura Blackie, Miguel-Aliaga lab, Francis Crick Institute 2025
# Part 3 intensity along length functions

library(tidyverse)
library(fields)
library(Cairo)

extractIntensityData <- function(centrelinefiles,pxSize,MaleGeno1,MaleGeno2,MaleGeno3,FemaleGeno1,FemaleGeno2,FemaleGeno3) {
  
  AllProfiles <- data.frame(matrix(ncol = 5, nrow = 0))
  LMAdjProfiles <- data.frame(matrix(ncol = 5, nrow = 0))

  for (i in 1:length(centrelinefiles)){
    gutcentreline <- suppressMessages(read_csv(file=centrelinefiles[i]))
    gutcentreline <- select(gutcentreline, "Position", "Intensity")
    #adds columns for name, sex and genotype
    name <- strsplit(centrelinefiles[i], "_Intensity") %>% sapply('[',1)
    if (grepl(MaleGeno1,name)){
      sexGenotype = rep(c("MaleGeno1"),dim(gutcentreline)[1])
    }
    if (grepl(FemaleGeno1, name)){
      sexGenotype = rep(c("FemaleGeno1"),dim(gutcentreline)[1])
    }
    if (grepl(MaleGeno2,name)){
      sexGenotype = rep(c("MaleGeno2"),dim(gutcentreline)[1])
    }
    if (grepl(FemaleGeno2, name)){
      sexGenotype = rep(c("FemaleGeno2"),dim(gutcentreline)[1])
    }
    if (grepl(MaleGeno3,name)){
      sexGenotype = rep(c("MaleGeno3"),dim(gutcentreline)[1])
    }
    if (grepl(FemaleGeno3, name)){
      sexGenotype = rep(c("FemaleGeno3"),dim(gutcentreline)[1])
    }
    namelist <- rep(c(name),dim(gutcentreline)[1])
    profile <- cbind(gutcentreline,sexGenotype,namelist)
    lmadjprofile <- cbind(gutcentreline,sexGenotype,namelist)
    cllength <- gutcentreline[dim(gutcentreline)[1],1]
    
    #finds position of landmarks and measures ditances between them along centreline
    xy <- suppressMessages(read_csv(file=paste0(name,"_XYCoordinatesOfLengthandDistances.csv")))
    Landmark <- suppressMessages(read_csv(file=paste0(name,"_R3coords.csv")))
    distance <- rdist(as.matrix(Landmark[1,2:3]/pxSize),xy[,2:3])
    minpoint <- which.min(distance)
    toLMlenxy <- sum(xy[1:minpoint,4])
    Landmarkpos <- which.min(abs(gutcentreline$Position-toLMlenxy))
    Landmarkdis <- gutcentreline$Position[Landmarkpos]
    lmadjprofile$Position[1:Landmarkpos] <- as.numeric(lmadjprofile$Position[1:Landmarkpos]) / as.numeric(Landmarkdis)*100
    lmadjprofile$Position[(Landmarkpos+1):length(lmadjprofile$Position)] <- (((as.numeric(lmadjprofile$Position[(Landmarkpos+1):length(lmadjprofile$Position)])-Landmarkdis) / as.numeric((cllength - Landmarkdis)))*100+100)
    
    #normalises for gut length for intensity along length plots
    profile$Position <- as.numeric(profile$Position) / as.numeric(cllength)*100
    AllProfiles <- rbind(AllProfiles,profile)
    LMAdjProfiles <- rbind(LMAdjProfiles,lmadjprofile)
  }
  
  write.csv(AllProfiles,file = "AllProfiles.csv")
  write.csv(LMAdjProfiles,file = "RegisteredProfiles_R3corrected.csv")
  
  return(list(AllProfiles=AllProfiles, LMAdjProfiles=LMAdjProfiles))
}

IntensityNnums <- function(centrelinefiles,genotype,MaleGeno1,MaleGeno2,MaleGeno3,FemaleGeno1,FemaleGeno2,FemaleGeno3) {
  
  nnums <- data.frame(matrix(ncol = 1, nrow = 0))
  for (i in 1:length(centrelinefiles)){
    name2 <- centrelinefiles[i]
    if (grepl(MaleGeno1,name2)){
      sexGenotype2 = "MaleGeno1"
    }
    if (grepl(FemaleGeno1, name2)){
      sexGenotype2 = "FemaleGeno1"
    }
    if (grepl(MaleGeno2,name2)){
      sexGenotype2 = "MaleGeno2"
    }
    if (grepl(FemaleGeno2, name2)){
      sexGenotype2 = "FemaleGeno2"
    }
    if (grepl(MaleGeno3,name2)){
      sexGenotype2 = "MaleGeno3"
    }
    if (grepl(FemaleGeno3, name2)){
      sexGenotype2 = "FemaleGeno3"
    }
    nnums <- rbind(nnums,sexGenotype2)}
  write.csv(as.data.frame(table(nnums)),file = paste0(genotype,"_Nnumbers.csv"))
}

intensityPlots <- function(AllProfiles, LMAdjProfiles, genotype){
  #plot - mean and sd of intensity profile unregistered
  binsize = 2
  VirginControlPlot <- AllProfiles %>% dplyr::filter(sexGenotype == "FemaleGeno1") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleControlPlot <- AllProfiles %>% dplyr::filter(sexGenotype == "MaleGeno1") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  VirginTestPlot <- AllProfiles %>% dplyr::filter(sexGenotype == "FemaleGeno2") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleTestPlot <- AllProfiles %>% dplyr::filter(sexGenotype == "MaleGeno2") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  VirginTestPlot2 <- AllProfiles %>% dplyr::filter(sexGenotype == "FemaleGeno3") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleTestPlot2 <- AllProfiles %>% dplyr::filter(sexGenotype == "MaleGeno3") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  all <- bind_rows(list("FemaleGeno1" = VirginControlPlot,"FemaleGeno2"=VirginTestPlot, "FemaleGeno3" =VirginTestPlot2,"MaleGeno1" = MaleControlPlot,"MaleGeno2"=MaleTestPlot, "MaleGeno3" =MaleTestPlot2),.id = "id")
  eb <- aes(ymax = Intensity_mean + Intensity_sd, ymin = Intensity_mean - Intensity_sd)
  CairoPDF(paste0(genotype,"_intensityAlongGutLength.pdf"), width = 8, height = 6)
  p1 = ggplot(data = all, aes(x = Position_mean, y = Intensity_mean, color=id, fill = id)) + 
    geom_line(size = 2) + 
    geom_ribbon(eb, alpha = 0.2, fill = "grey70", linetype = 0) +
    theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
          panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
          axis.line=element_line(colour='grey20'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())#,
  #legend.position="none")
  print(p1)
  dev.off()
  
  #plot - mean and sd of blue intensity profile REGISTERED by R3 position
  binsize = 2
  VirginControlPlot <- LMAdjProfiles %>% dplyr::filter(sexGenotype == "FemaleGeno1") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleControlPlot <- LMAdjProfiles %>% dplyr::filter(sexGenotype == "MaleGeno1") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  VirginTestPlot <- LMAdjProfiles %>% dplyr::filter(sexGenotype == "FemaleGeno2") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleTestPlot <- LMAdjProfiles %>% dplyr::filter(sexGenotype == "MaleGeno2") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  VirginTestPlot2 <- LMAdjProfiles %>% dplyr::filter(sexGenotype == "FemaleGeno3") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleTestPlot2 <- LMAdjProfiles %>% dplyr::filter(sexGenotype == "MaleGeno3") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  all <- bind_rows(list("FemaleGeno1" = VirginControlPlot,"FemaleGeno2"=VirginTestPlot, "FemaleGeno3" =VirginTestPlot2,"MaleGeno1" = MaleControlPlot,"MaleGeno2"=MaleTestPlot, "MaleGeno3" =MaleTestPlot2),.id = "id")
  eb <- aes(ymax = Intensity_mean + Intensity_sd, ymin = Intensity_mean - Intensity_sd)
  CairoPDF(paste0(genotype,"_intensityAlongGutLength_R3adjust.pdf"), width = 8, height = 6)
  p3 = ggplot(data = all, aes(x = Position_mean, y = Intensity_mean, color=id, fill = id)) + 
    geom_line(size = 2) + 
    geom_ribbon(eb, alpha = 0.2, fill = "grey70", linetype = 0) +
    theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
          panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
          axis.line=element_line(colour='grey20'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())#,
  #legend.position="none")
  print(p3)
  dev.off()
}
