library(vcfR)
library(tidyr)
library(dplyr)

##########################################################################################################################################################################

calls <- function(x){
  
  # Function for extracting info from VCF file
  
  vcf <- read.vcfR(x)
  vcf.tidy <- vcfR2tidy(vcf)
  vcf.tidy$fix$Locus <- paste0(vcf.tidy$fix$CHROM,":",vcf.tidy$fix$POS)
  vcf.tidy$fix$AF <- vcf.tidy$gt$gt_AF
  return(vcf.tidy$fix)
  
}

##########################################################################################################################################################################
# Matthews Correlation Coefficient

mcc <- function(TP,TN,FP,FN) {
  MCC <- ( TP*TN - FP*FN ) /
    sqrt( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) )
}

##########################################################################################################################################################################
# Function to select positions that have a coverage of at least x and compare to Gold Standard

# GS      = Gold Standard SNPs (heterozygous; from exome)
# GS.pos  = set of callable exon positions (coverage>=20) in Gold Standard
# cov.pos = coverage table for all exon positions in RNA-seq (from GATK DepthOfCoverage)
# x       = minimum coverage in RNA-seq to use as threshold
# vars    = variants called on RNA-seq data (pass all MuTect filters)

compare.to.GS <- function(x, GS, GS.pos, cov.pos, vars){
  
  # Select exon positions callable in exome and in RNA (coverage>=threshold)
  callable <- cov.pos %>% filter(Total_Depth>=x) %>% filter(Locus %in% GS.pos$Locus)
  
  # Select variants at positions callable in both
  callable.GS <- GS %>% filter(Locus %in% callable$Locus)
  callable.vars <- vars %>% filter(Locus %in% callable$Locus)
  
  # Put all variant calls into one data frame (one row per exon position)
  callable.all <- full_join(callable.GS,callable.vars, by=c("CHROM","POS"))
  
  # True Positive   = variant called in both exome and RNA-seq
  TP <- nrow(callable.all[!is.na(callable.all$ChromKey.x) &
                           !is.na(callable.all$ChromKey.y) &
                           callable.all$REF.x==callable.all$REF.y &
                           callable.all$ALT.x==callable.all$ALT.y ,])
  
  # False Positive  = variant called in RNA-seq but not on exome
  # or variant (REF/ALT) differs between the two
  FP <- nrow(callable.all[is.na(callable.all$ChromKey.x) &
                            !is.na(callable.all$ChromKey.y) ,]) +
    nrow(callable.all[!is.na(callable.all$ChromKey.x) &
                        !is.na(callable.all$ChromKey.y) &
                        callable.all$REF.x!=callable.all$REF.y |
                        callable.all$ALT.x!=callable.all$ALT.y ,])
  
  # False Negative  = Gold-Standard SNP not detected in RNA-seq
  FN <- nrow(callable.all[!is.na(callable.all$ChromKey.x) &
                            is.na(callable.all$ChromKey.y) ,])
  
  out <- data.frame("min.Coverage"=x,
                   "Total"=nrow(callable),
                   "GS.pos"=nrow(callable.GS),
                   "TP"=TP,
                   "FP"=FP,
                   "FN"=FN,
                   "TN"=nrow(callable)-TP-FP-FN,
                   "MCC"=mcc(TP,nrow(callable)-TP-FP-FN,FP,FN))
  
  return(out)
  
}

##########################################################################################################################################################################
# Make table of comparison to Gold Standard for a range of minimum consensus sequence coverages

comp.table <- function(min,max, GS, GS.pos, cov.pos, vars) {
  
  comp <- lapply(min:max, compare.to.GS, GS, GS.pos, cov.pos, vars)
  comp.table <- comp[[1]]
  for (n in 2:length(comp)){comp.table <- rbind(comp.table,comp[[n]])}
  
  return(comp.table)
  
}

##########################################################################################################################################################################
# Gold Standard coverage and SNPs

# Coverage of exome data on all exon positions
exome.cov.P1 <- read.table("exome.all_exons.DoC.P1.recal_reads", header=TRUE,sep="\t")
exome.cov.P2 <- read.table("exome.all_exons.DoC.P2.recal_reads", header=TRUE,sep="\t")
exome.cov.P4 <- read.table("exome.all_exons.DoC.P4.recal_reads", header=TRUE,sep="\t")

# Combine all coverage data
exome.cov <- data.frame(exome.cov.P1[c(1,4)],
                        Depth_for_P2=exome.cov.P2$Depth_for_P2,
                        Depth_for_P4=exome.cov.P4$Depth_for_P4)
exome.cov$min.Depth <- apply(exome.cov[2:4], 1, min)

# Select positions with coverage of at least 20 in all three samples and save/load as RDS
pos.GS <- exome.cov[exome.cov$min.Depth>=20,]
saveRDS(pos.GS, file="exome_exonPos_minCoverage20.rds")
#pos.GS <- readRDS("exome_exonPos_minCoverage20.rds")

# Get SNP information from exome data
snp <- readRDS("gold_het_snps.rds")
var.snp <- snp$fix
var.snp$Locus <- paste0(var.snp$CHROM,":",var.snp$POS)

##########################################################################################################################################################################
# MuTect2 without UMIs

# Load coverage info for all exon positions (from GATK DepthOfCoverage)
cov.MT <- data.frame(read.table("/coverage/recalibrated",
                                header=TRUE, sep="\t"))
# Load all variants called by MuTect2 on pre-processed reads; select only "heterozygous" ones that pass all filters
var.MT <- calls("MuTect.vcf")
var.MT <- var %>% filter(FILTER=="PASS")

# Make table of comparison to Gold Standard for a range of minimum RNA coverages
ct.MuTect <- comp.table(4,30, var.snp, pos.GS, cov.MT, var.MT)

##########################################################################################################################################################################
# MuTect2 without UMIs (dedupped)

cov.dedupped <- data.frame(read.table("/coverage/recalibrated.dedupped",
                                      header=TRUE, sep="\t"))

var.dd <- calls("MuTect_dedupped.vcf")
var.dd <- var.dd %>% filter(FILTER=="PASS")

# Make table of comparison to Gold Standard for a range of minimum RNA coverages
ct.dedupped <- comp.table(4,30, var.snp, pos.GS, cov.dedupped, var.dd)

##########################################################################################################################################################################
# MuTect2 on consensus sequences from UMI-tools grouped reads PER CELL taking out low MAPQ earlier

# Load coverage info for all exon positions (from GATK DepthOfCoverage)
cov.UMItools <- data.frame(read.table("/coverage/consensus.Aligned.recalibrated",
                                        header=TRUE, sep="\t"))
var.UMItools <- calls("consensus.MuTect.vcf")
var.UMItools <- var.UMItools %>% filter(FILTER=="PASS")

# Make table of comparison to Gold Standard for a range of minimum RNA coverages
ct.UMItools <- comp.table(4,30, var.snp, pos.GS, cov.UMItools, var.UMItools)

##########################################################################################################################################################################
# Plot False Discovery Rate

fdr.plot <- ggplot(NULL, aes(min.Coverage, FP/(TP+FP))) +
  geom_line(data=ct.MuTect, color="black", linetype=2) +
  geom_line(data=ct.dedupped, color="black", linetype=1) +
  geom_line(data=ct.UMItools, color="chartreuse3") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16),
        plot.title=element_text(size=20), legend.box.background=element_rect()) +
  labs(x="Minimum coverage", y="FDR")

ggsave(filename="FDR.png",fdr.plot, width=8,height=10)

##########################################################################################################################################################################
# Plot sensitivity

sens.plot <- ggplot(NULL, aes(min.Coverage, TP/GS.pos)) +
  geom_line(data=ct.MuTect, color="black", linetype=2) +
  geom_line(data=ct.dedupped, color="black", linetype=1) +
  geom_line(data=ct.UMItools, color="chartreuse3") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16),
        plot.title=element_text(size=20), legend.box.background=element_rect()) +
  labs(x="Minimum coverage", y="Sensitivity")

ggsave(filename="Sensitivity.png",sens.plot, width=8,height=8)

