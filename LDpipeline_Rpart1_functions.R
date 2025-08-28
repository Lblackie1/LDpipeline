# Pipeline for the Quantification and plotting of lipid droplets along gut length
# written by Laura Blackie, Miguel-Aliaga lab, Francis Crick Institute 2025
# Part 1 functions

library(tidyverse)
library(Cairo)
library(pracma)

rangesing <- function(x, na.rm=TRUE){
  minmax <- range(x, na.rm = TRUE)
  return(minmax[2]-minmax[1])
}

extractSpotData <- function(fileList) {

  allFliesR3 <- array()
  totalDrops <- array()
  fliestoremove <- array()
  flieswithmorebins <- array()

  for (i in 1:length(fileList)){
  
    savename = unlist(strsplit(fileList[i], "_DTspots"))
    savename = savename[1]
    savename = unlist(strsplit(savename, "Q_"))
    savename = savename[2]
  
    #reads in spot count data and tallys per bin
    boneJresults = suppressMessages(read_csv(fileList[i],show_col_types = FALSE))
    boneJmode = boneJresults$Mode
    spotCount <- aggregate(data.frame(count = boneJmode), list(value = boneJmode), length)
  
    #read in volumes and intensities
    spotVolumes <- suppressMessages(read_csv(paste0("M_",savename,"_DTspots_measurement_Results.csv"),show_col_types = FALSE))
    spotIntensity <- suppressMessages(read_csv(paste0("Q_",savename,"_DTspots_originalImage_quantification_Results.csv"),show_col_types = FALSE))
    spotQuants <- data.frame(Volume = spotVolumes$`Vol (unit)`)
    spotQuants$Intensity <- spotIntensity$Mean
    spotQuants$mode <- boneJresults$Mode
    spotMeans <- aggregate(data.frame(spotQuants), list(value = spotQuants$mode), mean)
    spotMeans <- dplyr::select(spotMeans, "value"="value", "Volume_mean"="Volume","Intensity_mean"="Intensity")
    spotVar <- aggregate(data.frame(spotQuants), list(value = spotQuants$mode), var)
    spotVar <- dplyr::select(spotVar, "value"="value", "Volume_var"="Volume","Intensity_var"="Intensity")
    spotRange <- aggregate(data.frame(spotQuants), list(value = spotQuants$mode), rangesing)
    spotRange <- dplyr::select(spotRange, "value"="value", "Volume_range"="Volume")
  
    #read in watershed profile and tallys per bin to calculate length and order of each bin
    WatershedProfile <- suppressMessages(read.delim(paste0(savename,"_gutMaskWatershedProfile.txt")))
    lineInt <- WatershedProfile$Value
    linInt16uni <- floor(0.5 + (unique(lineInt[!(lineInt %% 1)]) /max(lineInt)*65535))
    linInt16 <- floor(0.5 + (lineInt[!(lineInt %% 1)] /max(lineInt)*65535))
    linIntTally = data.frame(gutlentally = summary(factor(linInt16, level = linInt16uni)))
    linIntTally$value <- as.numeric(rownames(linIntTally))
    merged <- left_join(linIntTally, spotCount, by="value")
    merged <- left_join(merged, spotMeans, by="value")
    merged <- left_join(merged, spotVar, by="value")
    merged <- left_join(merged, spotRange, by="value")
    merged[is.na(merged)] <- 0
  
    final <- mutate(merged, gutlen=cumsum(merged$gutlentally), gutlenperc = cumsum(merged$gutlentally)/sum(gutlentally)*100, flyname = savename)
  
    if(any(linInt16uni==0)){
      print(paste0("Fly ",fileList[i]," has bin of 0"))
    }
  
    #adjusts for R3 position - find R3 bin and adjusts all the bins before it and all bins after it
    R3 <- suppressMessages(read_csv(paste0(savename,"_R3MaskVal.csv"),show_col_types = FALSE))
    R3val <- R3$Mode
    ind = which(unique(lineInt[!(lineInt %% 1)])==unique(R3val[!(R3val %% 1)]))
      
    gutlenperc1 = cumsum(merged$gutlentally[1:ind-1])/sum(merged$gutlentally[1:ind-1])*100
    gutlenperc2 = cumsum(merged$gutlentally[ind:length(merged$gutlentally)])/sum(merged$gutlentally[ind:length(merged$gutlentally)])*100 +100
    gutlenpercR3 = c(gutlenperc1,gutlenperc2)/2
    
    final2 <- cbind(final,gutlenpercR3)
    
    #saves table of bin lengths and no of lipid droplets per bin along gut length for each fly
    savetable <- dplyr::select(final, gutBinLength = gutlentally, NoOfLipidDroplets = count, MeanIntensity = Intensity_mean, MeanVolume.um = Volume_mean, VarIntensity=Intensity_var, VarVolume=Volume_var)
    #write.csv(savetable, file=paste0(savename, "_LipidDropletCountsIntensityVolumeBinnedAlongGutLength_largeDroplets.csv"))
    savetable2 <- dplyr::select(final2, gutBinLength = gutlentally, gutBinPercR3 = gutlenpercR3, NoOfLipidDroplets = count, MeanIntensity = Intensity_mean, MeanVolume.um = Volume_mean, VarIntensity=Intensity_var, VarVolume=Volume_var)
    #write.csv(savetable2, file=paste0(savename, "_LipidDropletCountsIntensityVolumeBinnedAlongGutLength_adjustR3_largeDroplets.csv"))
    
    #populates arrays for later plotting
    totalDrops <- rbind(totalDrops,c(savename, sum(final$count)))
    if(length(unique(lineInt[!(lineInt %% 1)])[!(unique(lineInt[!(lineInt %% 1)]) %in% 0)])==20){
      allFliesR3 <- rbind(allFliesR3, final2)
    } else if(length(unique(lineInt[!(lineInt %% 1)])[!(unique(lineInt[!(lineInt %% 1)]) %in% 0)])>=21){
      print(paste0("Fly ",fileList[i]," has ",length(unique(lineInt[!(lineInt %% 1)])), " bins"))
      flieswithmorebins <- rbind(flieswithmorebins,fileList[i])
    } else if(length(unique(lineInt[!(lineInt %% 1)])[!(unique(lineInt[!(lineInt %% 1)]) %in% 0)])<20){
      print(paste0("Fly ",fileList[i]," does not have enough bins"))
      fliestoremove <- rbind(fliestoremove,fileList[i])
    
    }
  }

  allFliesR3 <- na.omit(allFliesR3)
  if(length(fliestoremove)>1){
    fliestoremove <- na.omit(fliestoremove)
    write.csv(fliestoremove, file=paste0(genotype, "_fliesRemoved.csv"))
  }

  #saves csv file of total number of lipid droplets per fly
  totalDrops <- na.omit(totalDrops)
  colnames(totalDrops) <- c("FlyName", "TotalNoOfLipidDroplets")
  write.csv(totalDrops, file="TotalNoLipidDroplets_DTspots.csv")
  
  return(list(allFliesR3=allFliesR3, fliestoremove=fliestoremove))
}

assignGenoNnums <- function(allFliesR3, fliestoremove, MaleGeno1, FemaleGeno1, MaleGeno2, FemaleGeno2, MaleGeno3, FemaleGeno3) {
for(i in 1:nrow(allFliesR3)){
  if (grepl(MaleGeno2,allFliesR3$flyname[i])){ 
    allFliesR3$sex[i] = "MaleGeno2"
  } else if (grepl(FemaleGeno2,allFliesR3$flyname[i])){ 
    allFliesR3$sex[i] = "FemaleGeno2"
  } else if (grepl(MaleGeno3,allFliesR3$flyname[i])){ 
    allFliesR3$sex[i] = "MaleGeno3"
  } else if (grepl(FemaleGeno3,allFliesR3$flyname[i])){ 
    allFliesR3$sex[i] = "FemaleGeno3"
  } else if (grepl(MaleGeno1,allFliesR3$flyname[i])){ 
    allFliesR3$sex[i] = "MaleGeno1"
  } else if (grepl(FemaleGeno1,allFliesR3$flyname[i])){ 
    allFliesR3$sex[i] = "FemaleGeno1"
  }
}
allFliesR3 = allFliesR3[!(allFliesR3$value %in% 0),]
if(!is.na(fliestoremove[1])){
  allFliesR3$bins <- c(rep(c(1:20),(length(fileList)-length(fliestoremove))))
} else{
  allFliesR3$bins <- c(rep(c(1:20),(length(fileList))))
}
write.csv(as.data.frame(table(allFliesR3$sex[allFliesR3$bins == 20])),file = paste0(genotype,"_LDcounts_Nnumbers.csv"))

return(allFliesR3)
}

plotwithoutzero <- function(allFliesR3, genotype){
#plot with zeros removed
for (plotval in c("count", "Volume_mean", "Intensity_mean", "Volume_var", "Intensity_var", "Volume_range")){
  #plotval = "count"
  binsize = 5
  FemaleGeno1Plot <- allFliesR3 %>% filter(sex == "FemaleGeno1") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  MaleGeno1Plot <- allFliesR3 %>% filter(sex == "MaleGeno1") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  FemaleGeno2Plot <- allFliesR3 %>% filter(sex == "FemaleGeno2") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  MaleGeno2Plot <- allFliesR3 %>% filter(sex == "MaleGeno2") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  FemaleGeno3Plot <- allFliesR3 %>% filter(sex == "FemaleGeno3") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  MaleGeno3Plot <- allFliesR3 %>% filter(sex == "MaleGeno3") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  all <- bind_rows(list("FemaleGeno1" = FemaleGeno1Plot,"MaleGeno1" = MaleGeno1Plot,"FemaleGeno2" = FemaleGeno2Plot,"MaleGeno2" = MaleGeno2Plot,"FemaleGeno3" = FemaleGeno3Plot,"MaleGeno3" = MaleGeno3Plot),.id = "id")
  all[is.na(all)] <- 0
  eb <- aes(ymax = get(paste0(plotval,"_mean")) + get(paste0(plotval,"_sd")), ymin = get(paste0(plotval,"_mean")) - get(paste0(plotval,"_sd")))
  
  CairoPDF(paste0(genotype,"_meanLipidDroplet",plotval,"_DTspot_R3adjust_noZero.pdf"), width = 8, height = 6)
  p1 <- ggplot(data = all, aes(x = gutlenpercR3_mean, y = get(paste0(plotval,"_mean")), color=id, fill = id)) + 
    geom_line() + 
    geom_ribbon(eb, alpha = 0.2, fill = "grey70", linetype = 0) +
    scale_x_continuous(name="Distance along midgut length binned (% of total length adjusted by R3 alignment)") +
    scale_y_continuous(name=paste0("Mean ",plotval," of lipid droplets per bin")) +
    theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
          panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
          axis.line=element_line(colour='grey20'))#,
  #legend.position="none")
  print(p1)
  dev.off()
  
}
}

plotwithzero <- function(allFliesR3, genotype){
#plot with zero included
for (plotval in c("count", "Volume_mean", "Intensity_mean", "Volume_var", "Intensity_var", "Volume_range")){
  binsize = 5
  FemaleGeno1Plot <- allFliesR3 %>% filter(sex == "FemaleGeno1") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleGeno1Plot <- allFliesR3 %>% filter(sex == "MaleGeno1") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  FemaleGeno2Plot <- allFliesR3 %>% filter(sex == "FemaleGeno2") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleGeno2Plot <- allFliesR3 %>% filter(sex == "MaleGeno2") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  FemaleGeno3Plot <- allFliesR3 %>% filter(sex == "FemaleGeno3") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleGeno3Plot <- allFliesR3 %>% filter(sex == "MaleGeno3") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  all <- bind_rows(list("FemaleGeno1" = FemaleGeno1Plot,"MaleGeno1" = MaleGeno1Plot,"FemaleGeno2" = FemaleGeno2Plot,"MaleGeno2" = MaleGeno2Plot,"FemaleGeno3" = FemaleGeno3Plot,"MaleGeno3" = MaleGeno3Plot),.id = "id")
  all[is.na(all)] <- 0
  eb <- aes(ymax = get(paste0(plotval,"_mean")) + get(paste0(plotval,"_sd")), ymin = get(paste0(plotval,"_mean")) - get(paste0(plotval,"_sd")))
  
  CairoPDF(paste0(genotype,"_meanLipidDroplet",plotval,"_DTspot_R3adjust_withZero.pdf"), width = 8, height = 6)
  p1 <- ggplot(data = all, aes(x = gutlenpercR3_mean, y = get(paste0(plotval,"_mean")), color=id, fill = id)) + 
    geom_line() + 
    geom_ribbon(eb, alpha = 0.2, fill = "grey70", linetype = 0) +
    scale_x_continuous(name="Distance along midgut length binned (% of total length adjusted by R3 alignment)") +
    scale_y_continuous(name=paste0("Mean ",plotval," of lipid droplets per bin")) +
    theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
          panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
          axis.line=element_line(colour='grey20'))#,
  #legend.position="none")
  print(p1)
  dev.off()
}
}

exampleHistogram <- function(fileList, genotype) {
  i=2
  savenamea = unlist(strsplit(fileList[i], "_DTspots"))
  savenamea = savenamea[1]
  savenamea = unlist(strsplit(savenamea, "Q_"))
  savenamea = savenamea[2]
  #read in volumes and intensities for one gut
  spotVolumes2 <- suppressMessages(read_csv(paste0("M_",savenamea,"_DTspots_measurement_Results.csv"),show_col_types = FALSE))
  spotIntensity2 <- suppressMessages(read_csv(paste0("Q_",savenamea,"_DTspots_originalImage_quantification_Results.csv"),show_col_types = FALSE))
  spotQuants2 <- data.frame(Volume = spotVolumes2$`Vol (unit)`)
  spotQuants2$Intensity <- spotIntensity2$Mean

   #plot histogram of spot volumes
   ggplot(spotQuants2, aes(x=Volume))+
     geom_histogram()
  #plot scatter plot of volume and intensity
  CairoPNG(paste0(savenamea,"_scatterVolumeVSIntensity.png"), width = 800, height = 600)
  p5 <- 
    ggplot(spotQuants2, aes(x=Volume,y=Intensity))+
    geom_point() +
    scale_x_continuous(limits=c(0,1000), expand = c(0,0))
  print(p5)
  dev.off()
  return(p5)
}

extractSpotDataLarge <- function(fileList, spotCutOff) {

allFliesR3 <- array()
totalDrops <- array()
fliestoremove <- array()
flieswithmorebins <- array()

for (i in 1:length(fileList)){
  
  savename = unlist(strsplit(fileList[i], "_DTspots"))
  savename = savename[1]
  savename = unlist(strsplit(savename, "Q_"))
  savename = savename[2]
  
  #reads in spot count data and tallys per bin
  boneJresults = suppressMessages(read_csv(fileList[i],show_col_types = FALSE))
  boneJmode = boneJresults$Mode
  spotCount <- aggregate(data.frame(count = boneJmode), list(value = boneJmode), length)
  
  #read in volumes and intensities
  spotVolumes <- suppressMessages(read_csv(paste0("M_",savename,"_DTspots_measurement_Results.csv"),show_col_types = FALSE))
  spotIntensity <- suppressMessages(read_csv(paste0("Q_",savename,"_DTspots_originalImage_quantification_Results.csv"),show_col_types = FALSE))
  spotQuants <- data.frame(Volume = spotVolumes$`Vol (unit)`)
  spotQuants$Intensity <- spotIntensity$Mean
  spotQuants$mode <- boneJresults$Mode
  
  #_____________________********************************
  #filter to keep only large spots
  spotQuantsAll <- spotQuants
  spotQuantsLarge <- spotQuants[spotQuants$Volume >= spotCutOff,]
  spotQuants <- spotQuantsLarge
  #_____________________********************************
  
  #spotCountAll is a count of all the spots, the rest (intensity, volume, etc) are for the large spots
  spotCountAll <- aggregate(data.frame(countAll = spotQuantsAll$mode), list(value = spotQuantsAll$mode), length)
  if (dim(spotQuantsLarge)[1] == 0){
    spotCountLarge <- data.frame(value=0,count=0)
    spotMeans <- data.frame(value=0,Volume_mean=0,Intensity_mean=0)
    spotVar <- data.frame(value=0,Volume_var=0,Intensity_var=0)
  } else {
    spotCountLarge <- aggregate(data.frame(count = spotQuantsLarge$mode), list(value = spotQuantsLarge$mode), length)
    spotMeans <- aggregate(data.frame(spotQuants), list(value = spotQuants$mode), mean)
    spotMeans <- dplyr::select(spotMeans, "value"="value", "Volume_mean"="Volume","Intensity_mean"="Intensity")
    spotVar <- aggregate(data.frame(spotQuants), list(value = spotQuants$mode), var)
    spotVar <- dplyr::select(spotVar, "value"="value", "Volume_var"="Volume","Intensity_var"="Intensity")
  }
  
  #read in watershed profile and tallys per bin to calculate length and order of each bin
  WatershedProfile <- suppressMessages(read.delim(paste0(savename,"_gutMaskWatershedProfile.txt")))
  lineInt <- WatershedProfile$Value
  linInt16uni <- floor(0.5 + (unique(lineInt[!(lineInt %% 1)]) /max(lineInt)*65535))
  linInt16 <- floor(0.5 + (lineInt[!(lineInt %% 1)] /max(lineInt)*65535))
  linIntTally = data.frame(gutlentally = summary(factor(linInt16, level = linInt16uni)))
  linIntTally$value <- as.numeric(rownames(linIntTally))
  merged <- left_join(linIntTally, spotCountAll, by="value")
  merged <- left_join(merged, spotCountLarge, by="value")
  merged <- left_join(merged, spotMeans, by="value")
  merged <- left_join(merged, spotVar, by="value")
  merged[is.na(merged)] <- 0
  
  final <- mutate(merged, gutlen=cumsum(merged$gutlentally), gutlenperc = cumsum(merged$gutlentally)/sum(gutlentally)*100, flyname = savename, countPct = count/countAll*100)
  final$countPct[is.na(final$countPct)] <- 0
  
  if(any(linInt16uni==0)){
    print(paste0("Fly ",fileList[i]," has bin of 0"))
  }
  
  #adjusts for R3 position - find R3 bin and adjusts all the bins before it and all bins after it
  R3 <- suppressMessages(read_csv(paste0(savename,"_R3MaskVal.csv"),show_col_types = FALSE))
  R3val <- R3$Mode
  ind = which(unique(lineInt[!(lineInt %% 1)])==unique(R3val[!(R3val %% 1)]))
  
  gutlenperc1 = cumsum(merged$gutlentally[1:ind-1])/sum(merged$gutlentally[1:ind-1])*100
  gutlenperc2 = cumsum(merged$gutlentally[ind:length(merged$gutlentally)])/sum(merged$gutlentally[ind:length(merged$gutlentally)])*100 +100
  gutlenpercR3 = c(gutlenperc1,gutlenperc2)/2
  
  final2 <- cbind(final,gutlenpercR3)
  
  #saves table of bin lengths and no of lipid droplets per bin along gut length for each fly
  savetable <- dplyr::select(final, gutBinLength = gutlentally, NoOfLipidDroplets = count, MeanIntensity = Intensity_mean, MeanVolume.um = Volume_mean, VarIntensity=Intensity_var, VarVolume=Volume_var)
  #write.csv(savetable, file=paste0(savename, "_LipidDropletCountsIntensityVolumeBinnedAlongGutLength_largeDroplets.csv"))
  savetable2 <- dplyr::select(final2, gutBinLength = gutlentally, gutBinPercR3 = gutlenpercR3, NoOfLipidDroplets = count, MeanIntensity = Intensity_mean, MeanVolume.um = Volume_mean, VarIntensity=Intensity_var, VarVolume=Volume_var)
  #write.csv(savetable2, file=paste0(savename, "_LipidDropletCountsIntensityVolumeBinnedAlongGutLength_adjustR3_largeDroplets.csv"))
  
  #populates arrays and plots later plotting
  totalDrops <- rbind(totalDrops,c(savename, sum(final$count)))
  if(length(unique(lineInt[!(lineInt %% 1)])[!(unique(lineInt[!(lineInt %% 1)]) %in% 0)])==20){
    allFliesR3 <- rbind(allFliesR3, final2)
  } else if(length(unique(lineInt[!(lineInt %% 1)])[!(unique(lineInt[!(lineInt %% 1)]) %in% 0)])>=21){
    print(paste0("Fly ",fileList[i]," has ",length(unique(lineInt[!(lineInt %% 1)])), " bins"))
    flieswithmorebins <- rbind(flieswithmorebins,fileList[i])
  } else if(length(unique(lineInt[!(lineInt %% 1)])[!(unique(lineInt[!(lineInt %% 1)]) %in% 0)])<20){
    print(paste0("Fly ",fileList[i]," does not have enough bins"))
    fliestoremove <- rbind(fliestoremove,fileList[i])
    
  }
}

allFliesR3Large <- na.omit(allFliesR3)
if(length(fliestoremove)>1){
  fliestoremove <- na.omit(fliestoremove)
  fliestoremove
  write.csv(fliestoremove, file=paste0(genotype, "_fliesRemoved.csv"))
}

#saves csv file of total number of lipid droplets per fly
totalDrops <- na.omit(totalDrops)
colnames(totalDrops) <- c("FlyName", "TotalNoOfLipidDroplets")
write.csv(totalDrops, file="TotalNoLipidDroplets_DTspots_LargeDroplets.csv")

return(list(allFliesR3Large=allFliesR3Large, fliestoremoveLarge=fliestoremove))
}

plotwithoutzeroLarge <- function(allFliesR3, genotype){
#plot with zero removed
for (plotval in c("count", "Volume_mean", "Intensity_mean", "Volume_var", "Intensity_var", "countPct")){
  binsize = 5
  FemaleGeno1Plot <- allFliesR3 %>% filter(sex == "FemaleGeno1") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  MaleGeno1Plot <- allFliesR3 %>% filter(sex == "MaleGeno1") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  FemaleGeno2Plot <- allFliesR3 %>% filter(sex == "FemaleGeno2") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  MaleGeno2Plot <- allFliesR3 %>% filter(sex == "MaleGeno2") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  FemaleGeno3Plot <- allFliesR3 %>% filter(sex == "FemaleGeno3") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  MaleGeno3Plot <- allFliesR3 %>% filter(sex == "MaleGeno3") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% mutate(!!sym(plotval) := ifelse(get(plotval)==0, NA, get(plotval))) %>% summarise_all(list(mean=mean,sd=sd),na.rm=TRUE)
  all <- bind_rows(list("FemaleGeno1" = FemaleGeno1Plot,"MaleGeno1" = MaleGeno1Plot,"FemaleGeno2" = FemaleGeno2Plot,"MaleGeno2" = MaleGeno2Plot,"FemaleGeno3" = FemaleGeno3Plot,"MaleGeno3" = MaleGeno3Plot),.id = "id")
  all[is.na(all)] <- 0
  eb <- aes(ymax = get(paste0(plotval,"_mean")) + get(paste0(plotval,"_sd")), ymin = get(paste0(plotval,"_mean")) - get(paste0(plotval,"_sd")))
  
  CairoPDF(paste0(genotype,"_meanLipidDroplet",plotval,"_DTspot_R3adjust_LargeDroplets_noZero_",spotCutOff,".pdf"), width = 8, height = 6)
  p1 <- ggplot(data = all, aes(x = gutlenpercR3_mean, y = get(paste0(plotval,"_mean")), color=id, fill = id)) + 
    geom_line() + 
    geom_ribbon(eb, alpha = 0.2, fill = "grey70", linetype = 0) +
    scale_x_continuous(name="Distance along midgut length binned (% of total length adjusted by R3 alignment)") +
    scale_y_continuous(name=paste0("Mean ",plotval," of large lipid droplets (>",spotCutOff,"um^3) per bin")) +
    theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
          panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
          axis.line=element_line(colour='grey20'))#,
  #legend.position="none")
  print(p1)
  dev.off()
}
}

plotwithzeroLarge <- function(allFliesR3, genotype){
#plot with zeros included
for (plotval in c("count", "Volume_mean", "Intensity_mean", "Volume_var", "Intensity_var", "countPct")){
  binsize = 5
  FemaleGeno1Plot <- allFliesR3 %>% filter(sex == "FemaleGeno1") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleGeno1Plot <- allFliesR3 %>% filter(sex == "MaleGeno1") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  FemaleGeno2Plot <- allFliesR3 %>% filter(sex == "FemaleGeno2") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleGeno2Plot <- allFliesR3 %>% filter(sex == "MaleGeno2") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  FemaleGeno3Plot <- allFliesR3 %>% filter(sex == "FemaleGeno3") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  MaleGeno3Plot <- allFliesR3 %>% filter(sex == "MaleGeno3") %>% subset(select=c(gutlenpercR3,get(plotval))) %>% group_by(Position_bin = cut(gutlenpercR3, seq(0, 100, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
  all <- bind_rows(list("FemaleGeno1" = FemaleGeno1Plot,"MaleGeno1" = MaleGeno1Plot,"FemaleGeno2" = FemaleGeno2Plot,"MaleGeno2" = MaleGeno2Plot,"FemaleGeno3" = FemaleGeno3Plot,"MaleGeno3" = MaleGeno3Plot),.id = "id")
  all[is.na(all)] <- 0
  eb <- aes(ymax = get(paste0(plotval,"_mean")) + get(paste0(plotval,"_sd")), ymin = get(paste0(plotval,"_mean")) - get(paste0(plotval,"_sd")))
  
  CairoPDF(paste0(genotype,"_meanLipidDroplet",plotval,"_DTspot_R3adjust_LargeDroplets_withZero_",spotCutOff,".pdf"), width = 8, height = 6)
  p1 <- ggplot(data = all, aes(x = gutlenpercR3_mean, y = get(paste0(plotval,"_mean")), color=id, fill = id)) + 
    geom_line() + 
    geom_ribbon(eb, alpha = 0.2, fill = "grey70", linetype = 0) +
    scale_x_continuous(name="Distance along midgut length binned (% of total length adjusted by R3 alignment)") +
    scale_y_continuous(name=paste0("Mean ",plotval," of large lipid droplets (>",spotCutOff,"um^3) per bin")) +
    theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
          panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
          axis.line=element_line(colour='grey20'))#,
  #legend.position="none")
  print(p1)
  dev.off()
  
}
}

# Function to compute AUC for a region
compute_auc_in_region <- function(df, region, plotval) {
  region_data <- df %>%
    filter(gutlenpercR3 >= region[1], gutlenpercR3 <= region[2]) %>%
    arrange(gutlenpercR3)
  if (nrow(region_data) < 2) return(NA)  # Need at least two points
  region_data_sub <- region_data %>% subset(select=c(get(plotval)))
  auc <- trapz(region_data$gutlenpercR3, t(region_data_sub))
  return(auc)
}

compute_auc_all <- function(df, region1, region2,genotype,plotval,droplet){
  gut_aucs <- df %>%
    group_by(flyname, sex) %>%
    group_modify(~{
      data.frame(
        AUC_region1 = compute_auc_in_region(.x, R2region,plotval),
        AUC_region2 = compute_auc_in_region(.x, R4region,plotval)
      )
    }) %>%
    ungroup()
  
  write.csv(gut_aucs, paste0(genotype,"_",plotval,droplet,'_AreaUnderCurveR2R4.csv'))
  return(gut_aucs)
}

plot_stats_auc <- function(df,genotype,plotval,droplet){
  CairoPDF(paste0(genotype,"_",plotval,"_",droplet,"_intensityAlongGutLength_AreaUnderCurve_R2.pdf"), width = 8, height = 6)
  p1 = ggplot(data=gut_aucs, mapping = aes(x=sex, y=AUC_region1)) +
    geom_boxplot(outlier.colour = "red") +
    geom_jitter(width = 0.2, height = 0, size=1) +
    theme(panel.background=element_rect(fill='white'), plot.background = element_blank()) +
    theme(panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank()) +
    theme(axis.line=element_line(colour='grey20')) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    stat_summary(geom = "point", fun.y = "mean", col = "red", size = 10, shape = "-") +
    scale_x_discrete(name="") +
    scale_y_continuous(name=paste0(plotval," area under the curve R2"), expand = c(0, 0)) #, limits=c(0,8)
  print(p1)
  dev.off()
  CairoPDF(paste0(genotype,"_",plotval,"_",droplet,"_intensityAlongGutLength_AreaUnderCurve_R4.pdf"), width = 8, height = 6)
  p2 = ggplot(data=gut_aucs, mapping = aes(x=sex, y=AUC_region2)) +
    geom_boxplot(outlier.colour = "red") +
    geom_jitter(width = 0.2, height = 0, size=1) +
    theme(panel.background=element_rect(fill='white'), plot.background = element_blank()) +
    theme(panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank()) +
    theme(axis.line=element_line(colour='grey20')) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    stat_summary(geom = "point", fun.y = "mean", col = "red", size = 10, shape = "-") +
    scale_x_discrete(name="") +
    scale_y_continuous(name=paste0(plotval," area under the curve R4"), expand = c(0, 0)) #, limits=c(0,8)
  print(p2)
  dev.off()
  
  # Save ANOVA results to txt file
  sink(paste0(genotype,"_",plotval,"_",droplet,"_anova_results_areaUnderCurve.txt"))
  cat("Area under curve in R2 region:\n")
  totalanova <- aov(AUC_region1 ~ sex, data = gut_aucs)
  cat("ANOVA Summary:\n")
  print(summary(totalanova))
  cat("\nTukey HSD:\n")
  print(TukeyHSD(totalanova))
  cat("Area under curve in R4 region:\n")
  totalanova <- aov(AUC_region2 ~ sex, data = gut_aucs)
  cat("ANOVA Summary:\n")
  print(summary(totalanova))
  cat("\nTukey HSD:\n")
  print(TukeyHSD(totalanova))
  sink()
}

