#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_sbayesr {
    input:
        val chr
        path gwas_chr
        tuple val ld_path, path ld_files
        val out_prefix
    
    output:
        path "${out_prefix}_sBayesR_chr${chr}.snpRes"

    script:
        """
        echo "$gctb --sbayes R \
                --gwas-summary $gwas_chr \
                --ldm $ld_path \
                --gamma 0.0,0.01,0.1,1 \
                --pi 0.95.0.02,0.02,0.01 \
                --burn-in 2000 \
                --out-freq 10 \
                --out $out_prefix \
                --exclude-mhc \
                --impute-n"
             """
}