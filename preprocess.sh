#!/bin/bash

# Pre-process AML scRNA-seq data for UMI-consensus-based variant calling

# Trim adapter sequences
cutadapt -b file:/data/share/ngs/genomes/artificial_sequences/illumina_adapters.fa \
  -b file:/data/share/ngs/genomes/artificial_sequences/illumina_adapters_Rev.fa \
  -b file:/data/share/ngs/genomes/artificial_sequences/illumina_adapters_RevComp.fa \
  -b file:/data/share/ngs/genomes/artificial_sequences/SCRB.fa \
  -b file:/data/share/ngs/genomes/artificial_sequences/SCRB_Rev.fa \
  -b file:/data/share/ngs/genomes/artificial_sequences/SCRB_RevComp.fa \
  -a "A{100}" \
  -m 25 \
  -n 2 \
  -q 10,10 \
  -o cdnareads_cut.fq.gz \
  cdnareads.fq.gz

# Match reads with barcodes
perl fqcheck.pl barcodes.matched.fq.gz \
  cdnareads_cut.matched.fq.gz \
  sample_SC \
  sc_combined

# Paste barcodes to beginning of reads
#Unzip
gunzip barcodes.matched.fq.gz
gunzip cdnareads_cut.matched.fq.gz
#Paste
paste barcodes.matched.fq cdnareads_cut.matched.fq | awk '{if(NR%2==0){print $1$2;}else{print $1" "$2;}}' | sed 's/+ +/+/' | pigz -p 15 > reads_with_barcodes_cut.fq.gz
#Re-zip
pigz -p 15 barcodes.matched.fq barcodes.matched.fq.gz
pigz -p 15 cdnareads_cut.matched.fq cdnareads_cut.matched.fq.gz

# Move UMIs and cell barcodes (including i7) to read names
umi_tools extract --extract-method=string \
  --bc-pattern=CCCCCCCCCCCCCCNNNNNNNNNN \
  --log=processed.log \
  --stdin=reads_with_barcodes_cut.fq.gz \
  --stdout=reads_with_barcodes.processed.fq.gz
