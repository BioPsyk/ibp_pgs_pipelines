#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    input:
        tuple val(chr),
            path(gwas),
            val(N),
            path(ld_file),
            val(plink_prefix),       
            path(bed), 
            path(bim), 
            path(fam), 
            val(traitName)
    
    output:
    tuple val(chr),
        path("${traitName}_pst_eff_a1_b0.5_phiauto_chr22.txt")

    script:
        """
        echo "#Header\npython ${projectDir}/bin/PRScs.py --ref_dir=$workDir \
                --sst_file=$gwas \
                --bim_prefix=$plink_prefix \
                --n_gwas=$N \
                --chrom=$chr \
                --out_dir=$projectDir/$traitName" > ${traitName}_pst_eff_a1_b0.5_phiauto_chr22.txt
        """ 
}