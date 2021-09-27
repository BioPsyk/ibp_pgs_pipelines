#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process calc_score {
    publishDir launchDir
    label 'mod_mem'
    
    input:
        tuple val(chr),
        path(snp_posteriors),
        val(col_nums),
        val(trait),
        val(method),
        val(bfile),
        path(bed),
        path(bim),
        path(fam),
        path(plink)
    output:
        path "${trait}_${method}_chr${chr}.sscore"
    script:
        """
        ./plink2 --bfile ${bfile} \
        --out ${trait}_${method}_chr${chr} \
        --score ${snp_posteriors} ${col_nums} header cols=+scoresums ignore-dup-ids
        """
}