#! /usr/bin/env bash

curl -s https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel | grep YRI | cut -f1,1 > yri.ids

chip=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
seq=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
region=20:10000000-10100000

bcftools view -S yri.ids -r $region $seq -Ou | bcftools view -Ob -o truth.bcf
bcftools index truth.bcf
bcftools view -r $region $chip -Ob -o guess.bcf
bcftools index guess.bcf
