# ibp_pgs_pipelines

## Install a conda environment

```
# install nextflow using mamba (requires conda/mamba)
mamba create -n ibp_pgs_pipelines --channel bioconda \
  nextflow==20.10.0 \
  bcftools=1.9 \
  tabix

# Activate environment
conda activate ibp_pgs_pipelines

```

## This is a nextflow pipeline designed to calculate polygenic scores at IBP

The pipeline can be cloned along with all its dependencies using the following command:

`git clone --recursive https://github.com/vaqm2/ibp_pgs_pipelines.git`

The pipeline can be run as follows to display the parameter options and help message:

`nextflow run main.nf --help`

## The inputs required

* A reference GWAS in VCF format, indexed using tabix.
These can be obtained from the following resource: <https://gwas.mrcieu.ac.uk/datasets/>
* A json file containing paths to PLINK format `.bed/.bim/.fam` files.
These can be obtained for iPSYCH datasets from the `JSON/` folder that comes with the pipeline
* Linkage disequilibrium matrices for the SNPs used to calculate PGS
These are hardcoded in the `nextflow.config` file
* A phenotype file with the columns: `FID IID [Phenotype_Name]`
* A covariate file with the columns: `FID IID Age Sex PC1 PC2...`
This is hardcoded in the `nextflow.config` file for iPSYCH

The pipeline also requires the R packages `ggplot2, dplyr, fmsb, data.table` present in the environment
PGS computation using PRS-CS requires python packages `scipy, h5py` installed in the environment

## The outputs produced

* A file containing PGS from different methods.
For example, if your trait is named BMI, the score file will be named: `iPSYCH2012_All_Imputed_2021_QCed_BMI_Scores.txt`
with the columns:

```
    IID [ Sample identifer ]
    PT_5E8 [ PRSice pruning and thresholding score with SNPs P < 5E-8 ]
    PT_1E6 [ PRSice pruning and thresholding score with SNPs P < 1E-6 ]
    PT_0.05 [ PRSice pruning and thresholding score with SNPs P < 0.05 ]
    PT_1 [ PRSice pruning and thresholding score with all SNPs ]
    sBayesR_UKBB_2.8M  [ sBayesR scores using UKBB as reference LD and 2.8M high quality SNPs ]
    sBayesR_UKBB_HM3 [ sBayesR scores using UKBB as reference LD and HapMap3 SNPs ]
    PRSCS_UKBB_HM3 [ PRS-CS scores using UKBB as reference LD and HapMap3 SNPs ]
    PRSCS_1KG_HM3` [ PRS-CS scores using 1000 Genomes as reference LD and HapMap3 SNPs ]
```

* A file containing the variance explained in phenotype by each score with columns

```
    Method
    r2 [ Nagelkerke R2 ]
    P
```

* A bar plot for variance explained


## Generating example data from 1KGP


### Fetch source data
Align everything to GRCh37

```
# Download 1000 Genomes phase 3 release
mkdir -p 1000_Genomes_phase_3_release
cd 1000_Genomes_phase_3_release

Accessible from 1000 Genome data portal:
'https://www.internationalgenome.org/data-portal/data-collection/phase-3'

# Download required files
for chr in {1..22};do
  wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
done

# And their tabix index
for chr in {1..22};do
  wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
done

# Select the same example sumstats and ld files as in the IBP `ldsc_regression_rG` project
# This set is already GRCh37
cd ..
git clone git@github.com:BioPsyk/ldsc_regression_rG.git
mkdir example_data
cp -r ldsc_regression_rG/tests/example_data/test1/folder1/* example_data/
cp -r ldsc_regression_rG/tests/example_data/test1/folder2/* example_data/
cp -r ldsc_regression_rG/tests/example_data/test1/eur_w_ld_chr example_data/

```


### Find overlap between sumstats and 1KGP genotypes

```
# Extract location
mkdir extr_loc
for chr in $(seq 1 22); do \
 ( \
 echo "$chr starting ..."; \
  zcat 1000_Genomes_phase_3_release/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | awk '!/^#/{print $1":"$2, NR}' | LC_ALL=C sort -k 1,1 > extr_loc/chr${chr}_extracted_location; \
 echo "$chr done ..."; \
 ) & \
done; wait

# Prep one of the sumstats with location and index
# Reduce data to 10000 variants ( can be increased if needed )
seedval=1337
random_seed_file_source="random_seed_file_source"
openssl enc -aes-256-ctr -pass pass:"$seedval" -nosalt </dev/zero 2>/dev/null | head -10000 > ${random_seed_file_source}
zcat example_data/sumstat_2/sumstat_cleaned.gz | awk '{print $1":"$2, NR}' | shuf --random-source=${random_seed_file_source} -n10000 - | LC_ALL=C sort -k 1,1 > extr_loc/sumstat_2_extracted_location

# Find overlap between sumstats and 1KGP genotypes for each chromosome
mkdir overlap_loc
for chr in $(seq 1 22); do \
 ( \
 echo "$chr starting ..."; \
   LC_ALL=C join -1 1 -2 1 extr_loc/chr${chr}_extracted_location extr_loc/sumstat_2_extracted_location > overlap_loc/chr${chr}_joined_location; \
 echo "$chr done ..."; \
 ) & \
done; wait

# Select this set of variants from all source data
mkdir subset_1kgp
for chr in $(seq 1 22); do \
 ( \
 echo "$chr starting ..."; \
 head -n 3000 <(zcat 1000_Genomes_phase_3_release/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz) |  grep "#" > subset_1kgp/chr_${chr}_subset_1kgp;
 awk '
 NR==FNR{sel[$2]++;next}
 FNR in sel{print $0}
 ' overlap_loc/chr${chr}_joined_location <(zcat 1000_Genomes_phase_3_release/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz) >> subset_1kgp/chr_${chr}_subset_1kgp ; \
 bgzip -c subset_1kgp/chr_${chr}_subset_1kgp > subset_1kgp/chr_${chr}_subset_1kgp.gz; \
 tabix -p vcf subset_1kgp/chr_${chr}_subset_1kgp.gz; \
 rm subset_1kgp/chr_${chr}_subset_1kgp; \
 echo "$chr done ..."; \
 ) & \
done; wait

mkdir subset_eur_w_ld_chr
for i in $(seq 1 22);do
  echo $i
  awk -vOFS="\t" 'NR==FNR{a[$1];next}; FNR==1{print $0}; ($1":"$3 in a)' extr_loc/sumstat_2_extracted_location <(zcat example_data/eur_w_ld_chr/${i}.l2.ldscore.gz) | gzip -c > subset_eur_w_ld_chr/${i}.l2.ldscore.gz
# calculate new *.l2.M_5_50 files
  zcat subset_eur_w_ld_chr/${i}.l2.ldscore.gz | awk '$5>=0.05{print $0}' | wc -l > subset_eur_w_ld_chr/${i}.l2.M_5_50
done

mkdir subset_sumstats
for i in 2 3 5 6 7;do
  echo $i
  awk -vOFS="\t" 'NR==FNR{a[$1];next}; FNR==1{print $0}; ($1":"$2 in a)' extr_loc/sumstat_2_extracted_location <(zcat example_data/sumstat_$i/sumstat_cleaned.gz) | gzip -c > subset_sumstats/${i}.l2.ldscore.gz
done


```

# convert to the appropriate formats



