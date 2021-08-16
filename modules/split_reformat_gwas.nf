#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Nextflow process to split a gwas VCF into chunks
// Generates a space delimited file for prs-cs input
// Generates a file in mt-cojo format for sBayesR

process split_reformat_gwas {
    input:
        val chr
        val traitName
        path gwas
        val N
        each method

    output:
        val chr
        path "${traitName}_${method}_chr${chr}.txt"

    script:
        """
        python /bin/split_gwas_vcf.py --vcf $gwas --chromosome $chr --format $method
        """
}