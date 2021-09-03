#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    label 'big_mem'
    
    input:
        tuple val(chr),
            path(gwas),
            val(N),
            path(ld_bin),
            path(ld_info),
            val(cohort),
            val(plink_prefix),       
            path(bed), 
            path(bim), 
            path(fam), 
            val(traitName),
            path(prscs)
    
    output:
        path("${traitName}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt")

    script:
        """
        mkdir ${cohort}
        mv ${ld_bin} ${ld_info} ${cohort}/
        python ${prscs} --ref_dir=\$PWD/${cohort} \
            --sst_file=$gwas \
            --bim_prefix=$plink_prefix \
            --n_gwas=$N \
            --chrom=$chr \
            --out_dir=${traitName}
        """ 
}