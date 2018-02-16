#!/bin/bash

###################################################################################
#                                                                                 #
#              **** UMI-based scRNA-seq Consensus Variant Calls ****              #
#                                                                                 #
#  Calls consensus sequences using fgbio CallMolecularConsensusReads              #
#  and variants using GATK MuTect2.                                               #
#                                                                                 #
#  Takes as input FASTA file (-i) with cell barcode and UMI in read names.        #
#                                                                                 #
#  All other input files are hard-coded for now:                                  #
#  /data/ngs/genomes/Human/hg19/for_gatk/hg19_ucsc.fa                             #
#  /data/ngs/genomes/Human/hg19/for_gatk/STARidx_hg19_ucsc                        #
#  /data/ngs/genomes/Human/hg19/Homo_sapiens_onlyChr.chrM.ERCC.gfp.GRCh37.75.gtf  #
#  /data/ngs/genomes/Human/hg19/Homo_sapiens_onlyChr.ERCC.GRCh37.75.refFlat       #
#  /data/share/htp/Charlotte_UMI_SNV/contaminant_list.txt                         #
#  /data/ngs/genomes/Human/hg19/Human_polymorphism/ExAC/release0.3.1/NonTCGA      #
#                                         /ExAC_nonTCGA.r0.3.1.sites.vep.chr.vcf  #
#                                                                                 #
#  Output files have sample name (-s) in file names.                              #
#                                                                                 #
###################################################################################
#                                                                                 #
#  Tools used:                                                                    #
#                                                                                 #
#  UMI-tools 0.5.1                                                                #
#  STAR 2.5.3a                                                                    #
#  Samtools 1.1                                                                   #
#  fgbio 0.3.0                                                                    #
#  Genome Analysis ToolKit 3.6                                                    #
#                                                                                 #
###################################################################################


### 1) READ ARGUMENTS

# Get options
while getopts s:i: option
do
  case "${option}"
  in
  s) SAMPLE=${OPTARG};;
  i) FASTA=${OPTARG};;
  esac
done


### 2) MAP READS

# /data/share/htp/Charlotte_UMI_SNV/AML491/SC_McSCRB_ChristianeMouseman/GATK/alignment/get_STARidx.sh must have been run first

# Map reads to reference genome with STAR
STAR --genomeDir /data/ngs/genomes/Human/hg19/for_gatk/STARidx_hg19_ucsc \
  --runThreadN 10 \
  --readFilesCommand zcat \
  --sjdbGTFfile /data/ngs/genomes/Human/hg19/Homo_sapiens_onlyChr.chrM.ERCC.gfp.GRCh37.75.gtf \
  --outFileNamePrefix $SAMPLE. \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMmultNmax 1 \
  --outFilterMultimapNmax 50 \
  --outSAMunmapped Within \
  --sjdbOverhang 49 \
  --twopassMode Basic \
  --readFilesIn $FASTA


### 3) CONSENSUS

# Filter out reads with mapping quality <30
samtools view -q 30 -b -o $SAMPLE.Aligned.fitered.bam $SAMPLE.Aligned.sortedByCoord.out.bam

# Index mapped reads
samtools index $SAMPLE.Aligned.filtered.bam

# Sort/group by coordinate and UMI (only if they also have same cell barcode)
umi_tools group \
  -I $SAMPLE.Aligned.filtered.bam \
  --group-out=groups_perCell.tsv \
  --output-bam -S $SAMPLE.grouped.bam \
  --method=directional \
  --edit-distance-threshold=1 \
  --per-cell

# Call consensus reads
java -jar /data/share/htp/Charlotte_UMI_SNV/fgbio/fgbio-0.3.0.jar CallMolecularConsensusReads \
  -i $SAMPLE.grouped.bam \
  -o $SAMPLE.consensus.bam \
  -r $SAMPLE.consensus_rejects.bam --tag BX --min-reads 1


### 4) FASTQ

# Convert BAM to FASTQ format
samtools bam2fq $SAMPLE.consensus.bam > $SAMPLE.consensus.fq

# Run fastQC (optional)
mkdir -p fastq_consensus_$SAMPLE
fastqc -o fastq_consensus_$SAMPLE \
  -f fastq \
  -contaminants /data/share/htp/Charlotte_UMI_SNV/contaminant_list.txt \
  $SAMPLE.consensus.fq


### 5) MAP CONSENSUS SEQS

# STAR alignment
STAR --genomeDir /data/ngs/genomes/Human/hg19/for_gatk/STARidx_hg19_ucsc \
  --runThreadN 10 --sjdbGTFfile /data/ngs/genomes/Human/hg19/Homo_sapiens_onlyChr.chrM.ERCC.gfp.GRCh37.75.gtf \
  --outFileNamePrefix $SAMPLE.consensus. \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMmultNmax 1 \
  --outFilterMultimapNmax 50 \
  --outSAMunmapped Within \
  --sjdbOverhang 49 \
  --twopassMode Basic \
  --readFilesIn $SAMPLE.consensus.fq


### 6A) PICARD (optional)

# Get Picard metrics
picard-tools CollectMultipleMetrics I=$SAMPLE.consensus.Aligned.sortedByCoord.out.bam \
  O=$SAMPLE.consensus_alignment.multiple_metrics R=/data/ngs/genomes/Human/hg19/for_gatk/hg19_ucsc.fa
picard-tools CollectRnaSeqMetrics I=$SAMPLE.consensus.Aligned.sortedByCoord.out.bam \
  O=$SAMPLE.consensus_alingnment.RNAseq_Metrics \
  REF_FLAT=/data/ngs/genomes/Human/hg19/Homo_sapiens_onlyChr.ERCC.GRCh37.75.refFlat \
  STRAND=FIRST_READ_TRANSCRIPTION_STRAND #RIBOSOMAL_INTERVALS=/data/share/htp/Charlotte_UMI_SNV/hg38.rRNA.interval_list


### 6B) PREPARE FOR VARIANT CALLING

# Add read groups
picard-tools AddOrReplaceReadGroups I=$SAMPLE.consensus.Aligned.sortedByCoord.out.bam \
  O=$SAMPLE.consensus.Aligned.RG.bam \
  SO=coordinate \
  RGLB=McSCRB \
  RGPL=illumina \
  RGPU=hiseq \
  RGSM=AML491

# Index
samtools index $SAMPLE.consensus.Aligned.RG.bam

# Split into exon segments and hard clip parts in intronic regions, reassign mapping quality
java -jar /opt/bin/GenomeAnalysisTK.jar -T SplitNCigarReads \
  -R /data/ngs/genomes/Human/hg19/for_gatk/hg19_ucsc.fa \
  -I $SAMPLE.consensus.Aligned.RG.bam \
  -o $SAMPLE.consensus.Aligned.split.bam \
  -rf ReassignOneMappingQuality \
  -RMQF 255 -RMQT 60 \
  -U ALLOW_N_CIGAR_READS

# Decompress VCF file (GATK tools do not take piped input!)
#gunzip -k /data/ngs/genomes/Human/hg19/Human_polymorphism/ExAC/release0.3.1/NonTCGA/ExAC_nonTCGA.r0.3.1.sites.vep.chr.vcf.gz

## BQSR

# Analyze covariation before recalibration
java -jar /opt/bin/GenomeAnalysisTK.jar -T BaseRecalibrator \
  -R /data/ngs/genomes/Human/hg19/for_gatk/hg19_ucsc.fa \
  -I $SAMPLE.consensus.Aligned.split.bam \
  -knownSites /data/ngs/genomes/Human/hg19/Human_polymorphism/ExAC/release0.3.1/NonTCGA/ExAC_nonTCGA.r0.3.1.sites.vep.chr.vcf \
  -o recal_data.table

# Analyze covariation after recalibration
java -jar /opt/bin/GenomeAnalysisTK.jar -T BaseRecalibrator \
  -R /data/ngs/genomes/Human/hg19/for_gatk/hg19_ucsc.fa \
  -I $SAMPLE.consensus.Aligned.split.bam \
  -knownSites /data/ngs/genomes/Human/hg19/Human_polymorphism/ExAC/release0.3.1/NonTCGA/ExAC_nonTCGA.r0.3.1.sites.vep.chr.vcf \
  -BQSR recal_data.table \
  -o post_recal_data.table

# Plot covariation before and after recalibration
java -jar /opt/bin/GenomeAnalysisTK.jar -T AnalyzeCovariates \
-R /data/ngs/genomes/Human/hg19/for_gatk/hg19_ucsc.fa \
-before recal_data.table \
-after post_recal_data.table \
-plots recalibration_plots.pdf

# Apply recalibration
java -jar /opt/bin/GenomeAnalysisTK.jar -T PrintReads \
  -R /data/ngs/genomes/Human/hg19/for_gatk/hg19_ucsc.fa \
  -I $SAMPLE.consensus.Aligned.split.bam \
  -BQSR recal_data.table \
  -o $SAMPLE.consensus.Aligned.recalibrated.bam

# Delete uncompressed VCF file
#rm /data/ngs/genomes/Human/hg19/Human_polymorphism/ExAC/release0.3.1/NonTCGA/ExAC_nonTCGA.r0.3.1.sites.vep.chr.vcf


### 7) CALL VARIANTS

java -jar /opt/bin/GenomeAnalysisTK.jar -T MuTect2 \
  -R /data/ngs/genomes/Human/hg19/for_gatk/hg19_ucsc.fa \
  -I:tumor $SAMPLE.consensus.Aligned.recalibrated.bam \
  -o $SAMPLE.consensus.MuTect.vcf
