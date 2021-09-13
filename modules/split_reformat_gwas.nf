#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Nextflow process to split a gwas VCF into chunks
// Generates a space delimited file for prs-cs input
// Generates a file in mt-cojo format for sBayesR

process split_reformat_gwas {
    label 'mod_mem'

    input:
        tuple val(chr),
            val(traitName),
            path(vcf),
            path(vcf_idx),
            val(N),
            val(method),
            path(split_gwas)

    output:
        if (method == "prsice") {
            path("${traitName}_${method}_chr${chr}.txt")
        }
        else {
            tuple val(chr), 
                path("${traitName}_${method}_chr${chr}.txt")
        }

    script:
        """
        python ${split_gwas} --vcf $vcf --chromosome $chr --format $method
        """
}