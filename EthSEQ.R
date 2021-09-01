#install.packages('EthSEQ', repos='http://cran.us.r-project.org')
#BiocManager::install("SNPRelate")
#models https://github.com/cibiobcg/EthSEQ_Data/tree/master/EthSEQ_Models
setwd("/uufs/chpc.utah.edu/common/home/u1123911/EthSEQ")
library(EthSEQ)

models.dir = file.path("Models/")
vcf.dir = file.path("VCF/")
out.dir = file.path("Results_RP1/")
Sample = "RP1.rs.filtered.final.vcf"


## Download genotype data in VCF format
#dir.create(data.dir)
#download.file("https://github.com/aromanel/EthSEQ_Data/raw/master/Sample_SS2.vcf",
#              destfile = file.path(data.dir,"Sample_SS2.vcf"))

## Run the analysis
ethseq.Analysis(
  target.vcf =  file.path(vcf.dir,Sample),
  out.dir = out.dir,
  model.available = "SS2.All",
  model.folder = models.dir,
  verbose=TRUE,
  composite.model.call.rate = 1,
  space = "2D") # Default space is 2D

# Load and display computed ethnicity annotations
ethseq.annotations = read.delim(file.path(out.dir,"Report.txt"),sep="\t",as.is=TRUE,header=TRUE)
head(ethseq.annotations)

