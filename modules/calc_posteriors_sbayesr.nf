#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_sbayesr {
    input:
        tuple val(chr),
            path(gwas_chr),
            val(ld_prefix),
            path(ld_bin), 
            path(ld_info), 
            val(out_prefix)
    
    output:
        path("${out_prefix}_sbayesr_chr${chr}.snpRes")

    script:
        """
        ${projectDir}/bin/gctb --sbayes R \
            --gwas-summary ${gwas_chr} \
            --ldm ${ld_prefix} \
            --gamma 0.0,0.01,0.1,1 \
            --pi 0.95.0.02,0.02,0.01 \
            --burn-in 5000 \
            --chain-length 25000 \
            --out ${out_prefix}_sbayesr_chr${chr} \
            --exclude-mhc \
            --no-mcmc-bin \
            --thread 8 \
            --seed 80851 \
            --impute-n
        """
}