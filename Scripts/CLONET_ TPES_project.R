BiocManager::install("CLONETv2")
BiocManager::install("TPES")

library(data.table)
library(CLONETv2)
library(TPES)

setwd("C:/Users/OMISTAJA/OneDrive - TUNI.fi/Kurssit kesken/Computational human genomics/LAB sessions/10_Purity_ploidy_estimation/")

#For the tumor purity estimation, SNPs included are those identified as heterozygous around #6000 SNPs included
control = fread("Control_readcounts_for_CLONET.csv",data.table=F)
control$af = control$altCount/control$totalCount
tumor = fread("Tumor_readcounts_for_CLONET.csv",data.table=F)
tumor$af = tumor$altCount/tumor$totalCount

pileup.control = control[,c(1,2,4,5,14,8)]
colnames(pileup.control) = c("chr","pos","ref","alt","af","cov")

pileup.tumor = tumor[,c(1,2,4,5,14,8)]
colnames(pileup.tumor) = c("chr","pos","ref","alt","af","cov")

seg.tb <- fread("C:/Users/OMISTAJA/OneDrive - TUNI.fi/Kurssit kesken/Computational human genomics/LAB sessions/09_Somatic_copynumber_calling/SCNA.copynumber.called.seg",data.table=F)

bt <- compute_beta_table(seg.tb, pileup.tumor, pileup.control)

## Compute ploidy table with default parameters
pl.table <- compute_ploidy(bt)

adm.table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl.table)

allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
                                                                ploidy_table = pl.table, 
                                                                admixture_table = adm.table)


check.plot <- check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl.table,
                                         admixture_table = adm.table)
print(check.plot)


# TPES
snv.reads = fread("C:/Users/OMISTAJA/OneDrive - TUNI.fi/Kurssit kesken/Computational human genomics/LAB sessions/08_Somatic_Variant_calling/somatic.pm",data.table=F)
snv.reads = snv.reads[which(snv.reads$somatic_status=="Somatic"),]
snv.reads = snv.reads[,c("chrom","position","position","tumor_reads1","tumor_reads2")]
colnames(snv.reads) = c("chr","start","end","ref.count","alt.count")
snv.reads$sample = "Sample.1"

TPES_purity(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)

TPES_report(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)

