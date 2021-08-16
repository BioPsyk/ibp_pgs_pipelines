#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    input:
        val chr
        path gwas_chr
        val N
        path ld_mat
        tuple val(bfile), path(plink_files)
    
    output:
        path "${out_prefix}_prscs_chr${chr}.snpRes"

    script:
        """
        echo "python /bin/PRScs.py --ref_dir=$ld_mat \
                --sst_file=$gwas_chr \
                --bim_prefix=$bfile \
                --n_gwas=$N \
                --chrom=$chr \
                --out_dir=$workDir"
        """ 
}