
Compile:

```
make
bin/r2
```

A small test:

```
chip=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
seq=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
region=20:10000000-12000000

bcftools view -r $region $seq -Ob -o truth.bcf
bcftools index truth.bcf
bcftools view -r $region $chip -Ob -o guess.bcf
bcftools index guess.bcf

bin/r2 guess.bcf truth.bcf

bin/r2 --missing-as-homref guess.bcf truth.bcf ##Note there are no indels on the chip so VAR[GUESS]==0 hence covariance is undefined
```
