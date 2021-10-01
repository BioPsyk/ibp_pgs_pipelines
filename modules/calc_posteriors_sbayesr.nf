#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_sbayesr {
    label 'big_mem'
    
    input:
        tuple val(chr),
            path(gwas_chr),
            val(ld_prefix),
            path(ld_bin), 
            path(ld_info), 
            val(ld_cohort),
            val(ld_population),
            val(ld_snp_set),
            val(ld_format),
            val(trait),
            path(gctb)
    
    output:
        tuple val(chr),
            path("${trait}_${ld_cohort}_${ld_population}_${ld_snp_set}_sbayesr_chr${chr}.snpRes")

    script:
        """
        ./gctb --sbayes R \
            --gwas-summary ${gwas_chr} \
            --ldm ${ld_prefix} \
            --gamma 0.0,0.01,0.1,1 \
            --pi 0.95,0.02,0.02,0.01 \
            --burn-in 5000 \
            --chain-length 25000 \
            --out ${trait}_${ld_cohort}_${ld_population}_${ld_snp_set}_sbayesr_chr${chr} \
            --exclude-mhc \
            --no-mcmc-bin \
            --thread 8 \
            --seed 80851 \
            --impute-n
        """
}