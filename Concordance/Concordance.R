# 1075  bcftools merge E1_WES.vcf.gz E1_XLC.vcf.gz > E1_Combined.vcf
# 1082  bcftools view -Oz --max-alleles 2 --exclude-types indels E1_Combined.vcf.gz -o Intermediate.vcf.gz
# 1084  bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t\QUAL\t\FILTER\t\INFO\t\GT[\t%GT]\n' Intermediate.vcf.gz > Intermediate2.vcf
# 1091  bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t\QUAL\t\FILTER\t\INFO\t\GT[\t%GT]\n' E1_XLC.vcf.gz > E1_XLC.vcf
# 1093  bcftools view -Oz --max-alleles 2 --exclude-types indels E1_XLC.vcf.gz -o Intermediate.vcf.gz
# 1095  bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t\QUAL\t\FILTER\t\INFO\t\GT[\t%GT]\n' Intermediate.vcf.gz > E1_XLC.tab
# 1100  bcftools view -Oz --max-alleles 2 --exclude-types indels E1_WES.vcf.gz -o Intermediate.vcf.gz
# 1102  bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t\QUAL\t\FILTER\t\INFO\t\GT[\t%GT]\n' Intermediate.vcf.gz > E1_WES.tab


library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(grid)

WESname= "E1_WES.RP1"
XLCname= "E1_XLC.RP1"

XLC_File = read.table(file = XLCname,
                          header=TRUE, sep="\t", stringsAsFactors = FALSE, fill = T,
                          colClasses = c("character","character","character","character","character",
                                         "character","character","character","character","character") )
colnames(XLC_File) <- c("CHR","POS", "ID","REF","ALT","QUAL","FILTER","INFO","GT","GENOTYPE")
XLC_File$CHR <- gsub("chr","",XLC_File$CHR)

XLC_File$Key <-paste0(XLC_File$CHR,":",XLC_File$POS,XLC_File$REF,":",XLC_File$ALT) # make new column

WES_File = read.table(file = WESname,
                        header=TRUE, sep="\t", stringsAsFactors = FALSE, fill = T,
                        colClasses = c("character","character","character","character","character",
                                     "character","character","character","character","character") )
colnames(WES_File) <- c("CHR","POS", "ID","REF","ALT","QUAL","FILTER","INFO","GT","GENOTYPE")
WES_File$Key <-paste0(WES_File$CHR,":",WES_File$POS,WES_File$REF,":",WES_File$ALT) # make new column

merged_File <- merge(WES_File, XLC_File, by.x ="Key", by.y ="Key", all.x=TRUE)

merged_File <-merged_File[,c("Key","CHR.x","POS.x","REF.x","ALT.x","GENOTYPE.x","GENOTYPE.y")]
colnames(merged_File) <- c("Identifier","Chr","Pos","Ref","Alt","WES_GT","XLC_GT")

merged_File <- filter(merged_File, is.na(XLC_GT)==FALSE)

merged_File$Chr <- as.numeric(merged_File$Chr)
merged_File$Pos <- as.numeric(merged_File$Pos)

merged_File$XLC_GT <- ifelse(is.na(merged_File$XLC_GT)==TRUE,"0",merged_File$XLC_GT)
merged_File$XLC_GT <- ifelse(merged_File$XLC_GT=="0/1","1",merged_File$XLC_GT)
merged_File$XLC_GT <- ifelse(merged_File$XLC_GT=="1/1","2",merged_File$XLC_GT)

#merged_File$WES_GT <- ifelse(merged_File$WES_GT=="./.","0",merged_File$WES_GT)
merged_File <- filter(merged_File,merged_File$WES_GT!="./.")
merged_File$WES_GT <- ifelse(merged_File$WES_GT=="0/1","1",merged_File$WES_GT)
merged_File$WES_GT <- ifelse(merged_File$WES_GT=="1/1","2",merged_File$WES_GT)

merged_File$XLC_GT <- as.numeric(merged_File$XLC_GT)
merged_File$WES_GT <- as.numeric(merged_File$WES_GT)

merged_File$PosMB <- round(merged_File$Pos/1000000,digits = 4)
##ConcordanceMatrix
merged_File$Concordance <- (merged_File$WES_GT-merged_File$XLC_GT)*-1
merged_File$Concordance <- as.factor(merged_File$Concordance)

i=22
Chr<- filter(merged_File, Chr==i)
Chr <- group_by(Chr,Concordance)
Chr_summary <- summarize(Chr, Total=n())
Chr_summary$Percentage <- round(Chr_summary$Total/ sum(Chr_summary$Total)*100, digits=2)

print(paste0("Summary for Chromosome: ",i))  

plot1 <-  ggplot(Chr, aes(PosMB, PosMB)) + 
    geom_point(aes(colour = Concordance)) + scale_color_manual(values=c("green", "blue", "red","orange")) +
    ggtitle(paste0("Chr", i ," Correlation Plot"))
  
plot2 <-  ggplot(data=Chr, aes(x=Concordance, colour = Concordance)) + 
    geom_bar(aes(y = (..count..)/sum(..count..))) + scale_color_manual(values=c("green", "blue", "red","orange")) +
    scale_y_continuous(labels = percent) + theme_minimal() + ggtitle(paste0("Chr", i ," Concordance Count"))

pdf(paste0("Chr",i,".pdf"))
plot1
plot2
grid.newpage()
grid.table(Chr_summary)
dev.off()
