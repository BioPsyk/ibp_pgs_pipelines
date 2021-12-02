# ibp_pgs_pipelines

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
    PT_5E8 [ Pruning and thresholding score with SNPs P \< 5E-8 ]
    PT_1E6 [ Pruning and thresholding score with SNPs P \< 1E-6 ]
    PT_0.05 [ Pruning and thresholding score with SNPs P \< 0.05 ]
    PT_1 [ Pruning and thresholding score with all SNPs ]
    sBayesR_UKBB_2.8M  [ sBayesR scores using UKBB as reference LD and 2.8M high quality SNPs ]
    sBayesR_UKBB_HM3 [ sBayesR scores using UKBB as reference LD and HapMap3 SNPs ]
    PRSCS_UKBB_HM3 [ PRSCS scores using UKBB as reference LD and HapMap3 SNPs ]
    PRSCS_1KG_HM3` [ PRSCS scores using 1000 Genomes as reference LD and HapMap3 SNPs ]
```

* A file containing the variance explained in phenotype by each score with columns

    Method
    r2 [ Nagelkerke R2 ]
    P

* A bar plot for variance explained
