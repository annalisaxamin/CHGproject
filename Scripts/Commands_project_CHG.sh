#!/bin/bash

echo “//////// Control ////////” > Control_stats.txt
echo “//////// Tumor ////////” > Tumor_stats.txt

echo “1. Properly paired reads ***************************************************” >> Control_stats.txt
echo “1. Properly paired reads ***************************************************” >> Tumor_stats.txt

samtools sort Control.bam > Control.sorted.bam
samtools index Control.sorted.bam

samtools sort Tumor.bam > Tumor.sorted.bam
samtools index Tumor.sorted.bam

echo Control Statistics > Control_stats.txt
echo Samtools sorted >> Control_stats.txt
samtools flagstat Control.sorted.bam >> Control_stats.txt

echo Tumor Statistics > Tumor_stats.txt
echo Samtools sorted >> Tumor_stats.txt
samtools flagstat Tumor.sorted.bam >> Tumor_stats.txt 

echo “2. Realign and recalibrate  ***************************************************” >> Control_stats.txt
echo “2. Realign and recalibrate ***************************************************” >> Tumor_stats.txt

echo +++++++ Realignment +++++++ >> Control_stats.txt
echo +++++++ Realignment +++++++ >> Tumor_stats.txt

java -jar ../Tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.bam -o realigner.intervals -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T IndelRealigner -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.bam -targetIntervals realigner.intervals -o Control.sorted.realigned.bam -L Captured_Regions.bed

samtools index Control.sorted.realigned.bam

samtools view Control.sorted.realigned.bam |grep OC | wc -l >> Control_stats.txt
samtools flagstat Control.sorted.realigned.bam >> Control_stats.txt

java -jar ../Tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.bam -o realigner.intervals -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T IndelRealigner -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.bam -targetIntervals realigner.intervals -o Tumor.sorted.realigned.bam -L Captured_Regions.bed

samtools index Tumor.sorted.realigned.bam

samtools view Tumor.sorted.realigned.bam | grep OC | wc -l >> Tumor_stats.txt
samtools flagstat Tumor.sorted.realigned.bam >> Tumor_stats.txt

echo +++++++ Recalibration +++++++ >> Control_stats.txt
echo +++++++ Recalibration +++++++ >> Tumor_stats.txt

java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -o recal.table -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T PrintReads -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.realigned.bam -BQSR recal.table -o Control.sorted.realigned.recalibrated.bam -L Captured_Regions.bed --emit_original_quals

#samtools index Control.sorted.realigned.recalibrated.bam

java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.realigned.recalibrated.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -BQSR recal.table -o after_recal.table -L Captured_Regions.bed


java -jar ../Tools/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ../Annotations/human_g1k_v37.fasta -before recal.table -after after_recal.table -csv recal.csv -plots Control_recal.pdf


samtools view Control.sorted.realigned.recalibrated.bam | grep OQ | wc -l >> Control_stats.txt
samtools flagstat Control.sorted.realigned.recalibrated.bam >> Control_stats.txt

java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -o recal.table -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T PrintReads -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.bam -BQSR recal.table -o Tumor.sorted.realigned.recalibrated.bam -L Captured_Regions.bed --emit_original_quals

#samtools index Tumor.sorted.realigned.recalibrated.bam

java -jar ../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.recalibrated.bam -knownSites ../Annotations/hapmap_3.3.b37.vcf -BQSR recal.table -o after_recal.table -L Captured_Regions.bed

java -jar ../Tools/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ../Annotations/human_g1k_v37.fasta -before recal.table -after after_recal.table -csv recal.csv -plots Tumor_recal.pdf

samtools view Tumor.sorted.realigned.recalibrated.bam | grep OQ | wc -l >> Tumor_stats.txt
samtools flagstat Tumor.sorted.realigned.recalibrated.bam >> Tumor_stats.txt

echo “3. heterozygous SNPs - GATK  ***************************************************” >> Control_stats.txt
echo “3. heterozygous SNPs - GATK  ***************************************************” >> Tumor_stats.txt

java -jar ../Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../Annotations/human_g1k_v37.fasta -I Control.sorted.realigned.recalibrated.bam -o Control.GATK.vcf -L Captured_Regions.bed

# Done. There were 1 WARN messages, the first 1 are repeated below. 
# WARN  12:50:44,486 InbreedingCoeff - Annotation will not be calculated. InbreedingCoeff # requires at least 10 unrelated samples. 
##### ^^^^ can safely ignore it

vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf Control.GATK.vcf --out Control.GATK --recode --recode-INFO-all

java -jar ../Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../Annotations/human_g1k_v37.fasta -I Tumor.sorted.realigned.recalibrated.bam -o Tumor.GATK.vcf -L Captured_Regions.bed

vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf Tumor.GATK.vcf --out Tumor.GATK --recode --recode-INFO-all

echo taking only the heterozygous SNPs
grep -E "(^#|0/1)" Control.GATK.vcf > Control.het.vcf
grep -E "(^#|0/1)" Tumor.GATK.vcf > Tumor.het.vcf

echo Annotating SNPs with ClinVar
java -Xmx4g -jar ../Tools/snpEff/snpEff.jar -v hg19kg Control.het.vcf -s Control.het.html > Control.het.ann.vcf 		# SNP database

java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar Annotate ../Annotations/hapmap_3.3.b37.vcf  Control.het.ann.vcf > Control.het.ann_2.vcf		# ancestry

java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar Annotate ../Annotations/clinvar_Pathogenic.vcf Control.het.ann_2.vcf > Control.het.ann_3.vcf

echo SNP in Control.het.ann_3.vcf >> Control_stats.txt
grep -v ‘^#’ Control.het.ann_3.vcf | wc -l  >> Control_stats.txt

grep -E "(^#|CLNSIG=Pathogenic)" Control.het.ann_3.vcf > Control.het.ann_3.CVPath.vcf
# we get 1 SNP

echo Tumor het to clinvar path
java -Xmx4g -jar ../Tools/snpEff/snpEff.jar -v hg19kg Tumor.het.vcf -s Tumor.het.html > Tumor.het.ann.vcf

java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar Annotate ../Annotations/hapmap_3.3.b37.vcf Tumor.het.ann.vcf > Tumor.het.ann_2.vcf

java -Xmx4g -jar ../Tools/snpEff/SnpSift.jar Annotate ../Annotations/clinvar_Pathogenic.vcf Tumor.het.ann_2.vcf > Tumor.het.ann_3.vcf

echo SNP in Tumor.het.ann_3.vcf >> Tumor_stats.txt
grep -v ‘^#’ Tumor.het.ann_3.vcf | wc -l  >> Tumor_stats.txt 

grep -E "(^#|CLNSIG=Pathogenic)" Tumor.het.ann_3.vcf > Tumor.het.ann_3.CVPath.vcf
# we get 1 SNP

echo Ethseq start
echo Control.sorted.realigned.recalibrated.bam > BAMs_List.txt
Rscript RunEthSEQ.R
echo Ethseq end

echo copy number calling
samtools mpileup -q 1 -f ../Annotations/human_g1k_v37.fasta Control.sorted.realigned.recalibrated.bam Tumor.sorted.realigned.recalibrated.bam | java -jar ../Tools/VarScan.v2.3.9.jar copynumber --output-file SCNA --mpileup 1

java -jar ../Tools/VarScan.v2.3.9.jar copyCaller SCNA.copynumber --output-file SCNA.copynumber.called

echo SPM calling
samtools mpileup -q 1 -f ../Annotations/human_g1k_v37.fasta Control.sorted.realigned.recalibrated.bam > Control.sorted.realigned.recalibrated.pileup

samtools mpileup -q 1 -f ../Annotations/human_g1k_v37.fasta Tumor.sorted.realigned.recalibrated.bam > Tumor.sorted.realigned.recalibrated.pileup

java -jar ../Tools/VarScan.v2.3.9.jar somatic Control.sorted.realigned.recalibrated.pileup Tumor.sorted.realigned.recalibrated.pileup --output-snp somatic.pm --output-indel somatic.indel --output-vcf 1

java -jar ../Tools/VarScan.v2.3.9.jar somatic Control.sorted.realigned.recalibrated.pileup Tumor.sorted.realigned.recalibrated.pileup --output-snp somatic.pm --output-indel somatic.indel

##### 7. Determine which DNA repair genes (use file DNA_Repair_Genes.bed) overlap both...
### 7.1. ...heterozygous deletions...
grep -E "(^chrom|del)" SCNA.copynumber.called > SCNA.copynumber.called.del
bedtools intersect -a DNA_Repair_Genes.bed -b SCNA.copynumber.called.del > bedtools/repair..copynumber.del.bed
bedtools intersect -a SCNA.copynumber.called.del -b DNA_Repair_Genes.bed > bedtools/copynumber.del..repair.bed

### 7.2. ...and heterozygous SNPs of the patient that are in Clinvar.
bedtools intersect -a DNA_Repair_Genes.bed -b Control.het.ann_3.CVPath.vcf > bedtools/repair..Control.het.clnvar.bed

### 7.3. [both]
bedtools intersect -a bedtools/repair..copynumber.del.bed -b bedtools/repair..Control.het.clnvar.bed > bedtools/repair..copynumber.del..Control.het.clnvar.bed

##### 8. Determine which DNA repair genes overlap both...
### 8.1 ...heterozygous deletions...
# already performed in 7.1

### 8.2 ...and somatic point mutations of the patient...
grep -E "(^#|SOMATIC)" somatic.pm.vcf > somatic.somatic.pm.vcf
bedtools intersect -a DNA_Repair_Genes.bed -b somatic.somatic.pm.vcf > bedtools/repair..somatic.somatic.pm.bed

### 8.3. [both]
bedtools intersect -a bedtools/repair..copynumber.del.bed -b bedtools/repair..somatic.somatic.pm.bed > bedtools/repair..copynumber.del..somatic.somatic.pm.bed
