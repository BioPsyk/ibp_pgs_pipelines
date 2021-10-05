#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process calc_score {
    label 'mod_mem'
    
    input:
        tuple val(chr),
        path(snp_posteriors),
        val(col_nums),
        val(out_prefix),
        val(bfile),
        path(bed),
        path(bim),
        path(fam),
        val(target_cohort),
        val(target_population),
        val(target_snp_set),
        path(plink)

    output:
        path "${target_cohort}_${target_population}_${target_snp_set}_${out_prefix}_chr${chr}.sscore"
    
    script:
        """
        ./plink2 --bfile ${bfile} \
        --out ${target_cohort}_${target_population}_${target_snp_set}_${out_prefix}_chr${chr} \
        --score ${snp_posteriors} ${col_nums} header cols=+scoresums ignore-dup-ids

        awk '{gsub(/^\#/, ""); print}' ${target_cohort}_${target_population}_${target_snp_set}_${out_prefix}_chr${chr}.sscore > tmp.score
        mv tmp.score ${target_cohort}_${target_population}_${target_snp_set}_${out_prefix}_chr${chr}.sscore
        """
}