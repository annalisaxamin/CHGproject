#Project question 10.
library(data.table)
install.packages("SPIAssay")
library(SPIAssay)
setwd("~/Documents/HumanGenomics/Project_2022/")

normal = fread("Control.csv",data.table=F)
normal$af = normal$altCount/normal$totalCount
tumor = fread("Tumor.csv",data.table=F)
tumor$af = tumor$altCount/tumor$totalCount
#AF = alternative base / total counts
#How do I determine the genotype from the AF, 0/1/2?
pileup.normal = normal[,c(1,2,4,5,14,8)]
colnames(pileup.normal) = c("chr","pos","ref","alt","af","cov")

pileup.tumor = tumor[,c(1,2,4,5,14,8)]
colnames(pileup.tumor) = c("chr","pos","ref","alt","af","cov")
                           

                           
# CODING: 0 for homozygous ref, 1 for heterozygous, 2 for homozygous alt
for (j in 1:nrow(pileup.normal) )  
{
  if(pileup.normal$af[j]>0.85) {pileup.normal$genotype[j] <- 2}
  else if(pileup.normal$af[j]<0.25) {pileup.normal$genotype[j]<-0} 
  else if(pileup.normal$af[j]>0.4  & pileup.normal$af[j]<0.6) {pileup.normal$genotype[j] <- 1} 
  else {pileup.normal$genotype[j] <- NA}
}

for (j in 1:nrow(pileup.tumor) )  
{
  if(pileup.tumor$af[j]>0.85) {pileup.tumor$genotype[j] <- 2}
  else if(pileup.tumor$af[j]<0.25) {pileup.tumor$genotype[j]<-0} 
  else if(pileup.tumor$af[j]>0.4  & pileup.tumor$af[j]<0.6) {pileup.tumor$genotype[j] <- 1} 
  else {pileup.tumor$genotype[j] <- NA}
}



