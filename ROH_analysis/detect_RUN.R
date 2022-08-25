#plink --bfile RP1_rs_filtered --keep keep.txt --make-bed --out frequent --snps-only --maf 0.25 --hwe 0.00000001 --geno 0.05 --me 0.5 0.1
#plink --noweb --allow-no-sex --bfile frequent --indep-pairwise 5kb 5 0.1 --out prunedsnplist
#plink --noweb --allow-no-sex --bfile frequent --extract prunedsnplist.prune.in --make-bed --out ROH_analysis --snps-only
#plink --bfile ROH_analysis --homozyg

setwd("/uufs/chpc.utah.edu/common/home/u1123911/Linkage_RP1/ROH")
#install.packages("detectRUNS")
library("detectRUNS")
#install.packages("ggplot2")
#install.packages("stargazer")
library(stargazer)


genotypeFilePath <- system.file(
  "extdata", "Kijas2016_Sheep_subset.ped", package="detectRUNS")

genotypeFilePath <- paste0(getwd(),"/","ROH_data.ped")

mapFilePath <- paste0(getwd(),"/","ROH_data.map")

slidingRuns <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  windowSize = 15, 
  threshold = 0.05,
  minSNP = 20, 
  ROHet = FALSE, 
  maxOppWindow = 1, 
  maxMissWindow = 1,
  maxGap = 10^6, 
  minLengthBps = 250000, 
  minDensity = 1/10^3, # SNP/kbps
  maxOppRun = NULL,
  maxMissRun = NULL
) 

consecutiveRuns <- consecutiveRUNS.run(
  genotypeFile =genotypeFilePath,
  mapFile = mapFilePath,
  minSNP = 20,
  ROHet = FALSE,
  maxGap = 10^6,
  minLengthBps = 250000,
  maxOppRun = 1,
  maxMissRun = 1
)

summaryList <- summaryRuns(
  runs = slidingRuns, mapFile = mapFilePath, genotypeFile = genotypeFilePath, 
  Class = 6, snpInRuns = TRUE)

pdf("ROH_analysis.pdf", width = 15, height = 10)
plot_Runs(runs = slidingRuns)

plot_SnpsInRuns(
  runs = slidingRuns[slidingRuns$chrom==8,], genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath)

plot_SnpsInRuns(
  runs = slidingRuns[slidingRuns$chrom==9,], genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath)

dev.off()

pdf("ROH_manhattan.pdf", width = 15, height = 10)
plot_manhattanRuns(
  runs = consecutiveRuns, 
  genotypeFile = genotypeFilePath, mapFile = mapFilePath)

dev.off()
Inbreeding_total <- Froh_inbreeding(runs = slidingRuns,mapFile =mapFilePath, genome_wide=TRUE)
Inbreeding_per_Chr <- Froh_inbreeding(runs = slidingRuns,mapFile =mapFilePath, genome_wide=FALSE)

Froh_inbreedingClass(runs = slidingRuns,mapFile =mapFilePath, Class = 2)


##Reports
stargazer(Inbreeding_total, summary=FALSE, rownames=FALSE,
          type = "text", title="Frequency of ROH per individual RP1", digits=4, out="fROH_RP1_individual.txt")

stargazer(Inbreeding_per_Chr, summary=TRUE, rownames=TRUE,
          type = "text", title="Frequency of ROH per CHR RP1", digits=4, out="fROH_RP1_CHR.txt")
