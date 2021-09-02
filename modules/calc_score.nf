#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process calc_score {
    input:
        tuple path(snp_posteriors),
        val(rsId_col_num),
        val(eff_allele_col_num),
        val(eff_size_col_num),
        val(trait),
        val(method),
        val(plink_prefix),
        path(bed),
        path(bim),
        path(fam),
        path(plink)
    output:
        path "${trait}_${method}.sscore"
    script:
        """
        ${plink} --bfile ${plink_prefix} \
        --out ${trait}_${method} \
        --score ${snp_posteriors} ${rsId_col_num} ${eff_allele_col_num} ${eff_size_col_num}
        """
}