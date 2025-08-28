# Pipeline for the Quantification and plotting of lipid droplets along gut length
# written by Laura Blackie, Miguel-Aliaga lab, Francis Crick Institute 2025
# Part 2 run script

#change working directory to folder where output from FIJI macro are saved
setwd("FIJI/output/folder/LipidDropletCount")

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

#select regions
R2region <- c(25,40)
R4region <- c(70,85)

#Gets list of files in folder
fileList <- list.files(path=getwd(), pattern="\\DTspots_GutMaskMode_Results.csv", full.names=FALSE)

spots <- extractSpotData(fileList)
allFliesR3 = spots$allFliesR3
fliestoremove = spots$fliestoremove
genos <- assignGenoNnums(allFliesR3, fliestoremove, MaleGeno1, FemaleGeno1, MaleGeno2, FemaleGeno2, MaleGeno3, FemaleGeno3)
allFliesR3 = genos
write_csv(allFliesR3,paste0(genotype,"_LipidDropletQuants.csv"))
plotwithzero(allFliesR3, genotype)
plotwithoutzero(allFliesR3, genotype)

for (plotval in c("count", "Volume_mean", "Intensity_mean", "Volume_var", "Intensity_var", "Volume_range")){
  gut_aucs <- compute_auc_all(allFliesR3, R2region, R4region,genotype,plotval,"allDroplets")
  plot_stats_auc(gut_aucs,genotype,plotval,"allDroplets")
}

#Large spots
#Add value for cut off as large spots
spotCutOff = 100
p5 = exampleHistogram(fileList, genotype)
spotsLarge <- extractSpotDataLarge(fileList, spotCutOff)
allFliesR3Large = spotsLarge$allFliesR3Large
fliestoremoveLarge = spotsLarge$fliestoremoveLarge
genosLarge <- assignGenoNnums(allFliesR3Large, fliestoremoveLarge, MaleGeno1, FemaleGeno1, MaleGeno2, FemaleGeno2, MaleGeno3, FemaleGeno3)
allFliesR3Large = genosLarge
write_csv(allFliesR3Large,paste0(genotype,"_largeDroplets_LipidDropletQuants.csv"))
plotwithzeroLarge(allFliesR3Large, genotype)
plotwithoutzeroLarge(allFliesR3Large, genotype)

for (plotval in c("count", "Volume_mean", "Intensity_mean", "Volume_var", "Intensity_var", "countPct")){
  gut_aucs <- compute_auc_all(allFliesR3Large, R2region, R4region,genotype,plotval,"largeDroplets")
  plot_stats_auc(gut_aucs,genotype,plotval,"largeDroplets")
}
