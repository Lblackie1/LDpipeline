# Pipeline for the Quantification and plotting of lipid droplets along gut length and intensity along length
# written by Laura Blackie, Miguel-Aliaga lab, Francis Crick Institute 2025
# Part 4 run for intensity along length

#change working directory to folder where files from FIJI macro are saved
setwd("FIJI/output/folder/IntensityAlongLength")

#Genotype (used for naming saved plots)
genotype = "Genotype"

#Allows up to three genotypes to be compared at the same time all with males and females so 6 groups
#Add text that is found in and only in the file names for each group. Control is typically geno1. If you only have two genotypes, put 'null' in geno3.  
MaleGeno1 <- "Ctrl_VirginMale"
MaleGeno2 <- "Test_VirginMale"
MaleGeno3 <- "null"
FemaleGeno1 <- "Ctrl_VirginFemale"
FemaleGeno2 <- "Test_VirginFemale"
FemaleGeno3 <- "null"

#Add the xy pixel size here
pxSize <- 0.3786027

centrelinefiles <- list.files(path=getwd(), pattern="\\Centreline.csv", full.names=FALSE)
intensities <- extractIntensityData(centrelinefiles,pxSize,MaleGeno1,MaleGeno2,MaleGeno3,FemaleGeno1,FemaleGeno2,FemaleGeno3)
AllProfiles <- intensities$AllProfiles
LMAdjProfiles <- intensities$LMAdjProfiles
IntensityNnums(centrelinefiles,genotype,MaleGeno1,MaleGeno2,MaleGeno3,FemaleGeno1,FemaleGeno2,FemaleGeno3)
intensityPlots(AllProfiles, LMAdjProfiles, genotype)
