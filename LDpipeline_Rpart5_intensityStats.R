# Pipeline for the Quantification and plotting of lipid droplets along gut length and intensity along length
# written by Laura Blackie, Miguel-Aliaga lab, Francis Crick Institute 2025
# Part 5 statistics for intensity along length

library(tidyverse)
library(pracma)
library(Cairo)

# Function to compute AUC for a region
compute_auc_in_region <- function(df, region) {
  region_data <- df %>%
    filter(percentPosition >= region[1], percentPosition <= region[2]) %>%
    arrange(percentPosition)
  
  if (nrow(region_data) < 2) return(NA)  # Need at least two points
  auc <- trapz(region_data$percentPosition, region_data$Intensity)
  return(auc)
}

compute_auc_all <- function(df, region1, region2){
    gut_aucs <- LMAdjProfiles %>%
      group_by(namelist, sexGenotype) %>%
      group_modify(~{
        data.frame(
          AUC_region1 = compute_auc_in_region(.x, R2region),
          AUC_region2 = compute_auc_in_region(.x, R4region)
        )
      }) %>%
      ungroup()

  write.csv(gut_aucs, paste0(genotype,'_AreaUnderCurveR2R4.csv'))
  return(gut_aucs)
}

plot_stats_auc <- function(df,genotype){
CairoPDF(paste0(genotype,"_intensityAlongGutLength_AreaUnderCurve_R2.pdf"), width = 8, height = 6)
p1 = ggplot(data=gut_aucs, mapping = aes(x=sexGenotype, y=AUC_region1)) +
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(width = 0.2, height = 0, size=1) +
  theme(panel.background=element_rect(fill='white'), plot.background = element_blank()) +
  theme(panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_line(colour='grey20')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  stat_summary(geom = "point", fun.y = "mean", col = "red", size = 10, shape = "-") +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Area under the curve R2", expand = c(0, 0)) #, limits=c(0,8)
print(p1)
dev.off()
CairoPDF(paste0(genotype,"_intensityAlongGutLength_AreaUnderCurve_R4.pdf"), width = 8, height = 6)
p2 = ggplot(data=gut_aucs, mapping = aes(x=sexGenotype, y=AUC_region2)) +
  geom_boxplot(outlier.colour = "red") +
  geom_jitter(width = 0.2, height = 0, size=1) +
  theme(panel.background=element_rect(fill='white'), plot.background = element_blank()) +
  theme(panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank()) +
  theme(axis.line=element_line(colour='grey20')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  stat_summary(geom = "point", fun.y = "mean", col = "red", size = 10, shape = "-") +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Area under the curve R4", expand = c(0, 0)) #, limits=c(0,8)
print(p2)
dev.off()

# Save ANOVA results to txt file
sink(paste0(genotype,"_anova_results_areaUnderCurve.txt"))
cat("Area under curve in R2 region:\n")
totalanova <- aov(AUC_region1 ~ sexGenotype, data = gut_aucs)
cat("ANOVA Summary:\n")
print(summary(totalanova))
cat("\nTukey HSD:\n")
print(TukeyHSD(totalanova))
cat("Area under curve in R4 region:\n")
totalanova <- aov(AUC_region2 ~ sexGenotype, data = gut_aucs)
cat("ANOVA Summary:\n")
print(summary(totalanova))
cat("\nTukey HSD:\n")
print(TukeyHSD(totalanova))
sink()
}

#Move all "RegisteredProfiles_R3corrected.csv" files to one folder and add genotype name to front separated by underscore
setwd("Intensity/profiles/folder/")

#select regions
R2region <- c(25,40)
R4region <- c(70,85)

fileList <- list.files(path=getwd(), pattern="\\.csv", full.names=FALSE)

for(file in 1:length(fileList)){
  genotype = strsplit(fileList[file],"_")
  genotype = unlist(genotype)[1]

  LMAdjProfiles = read.csv(fileList[file])
  LMAdjProfiles$percentPosition <- LMAdjProfiles$Position/200*100

  gut_aucs <- compute_auc_all(LMAdjProfiles, R2region, R4region)
  plot_stats_auc(gut_aucs,genotype)
}

