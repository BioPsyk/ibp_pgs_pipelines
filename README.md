# ibp_pgs_pipelines

## This is a nextflow pipeline designed to calculate polygenic scores at IBP

The pipeline can be cloned along with all its dependencies using the following command:

`git clone --recursive https://github.com/vaqm2/ibp_pgs_pipelines.git`

The pipeline can be run as follows to display the parameter options and help message:

`nextflow run main.nf --help`

## The inputs required

* A reference GWAS in VCF format, indexed using tabix.
These can be obtained from the output of the IBP sum stats cleaning pipeline: `https://github.com/BioPsyk/cleansumstats` \
Please make sure that the VCF has standard errors in the format column as the older summary stats might be missing these.
The code ignores SNPs without standard errors. \
The file to use from the output directory is labeled: `cleaned_GRCh37_full_altbased.vcf.gz`
* A json file containing paths to PLINK format `.bed/.bim/.fam` files.
These can be obtained for iPSYCH datasets on iPSYCH-GDK and a 1000 Genomes test dataset for GDK open from the `JSON/` folder that comes with the pipeline
* Linkage disequilibrium matrices for the SNPs used to calculate PGS
These are hardcoded in the `nextflow.config` file for iPSYCH-GDK and `gdk_open.config` for GDK-open
* A phenotype file with the columns: `FID IID [Phenotype_Name]`
* A covariate file with the columns: `FID IID Age Sex PC1 PC2...`
This is hardcoded in the `nextflow.config` file for iPSYCH and `gdk_open.config` for GDK-open. 

* The pipeline requires the R packages `ggplot2, dplyr, fmsb, data.table, argparser`
* The pipeline requires the python packages `scipy, h5py, cyvcf2`

The file rquirements.txt produced by conda can be used to reproduce the developer environment to run this pipeline with the following command:
`conda create --name <env> --file requirements.txt`

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
    r2_L [ Liability adjusted R2 for case-control outcomes ]
    P
```

* A bar plot for variance explained
