#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for sBayesR

process sbayesR {
    input:
        val chr, file gwas_chr from sBayesR_chunk_ch
        path sBayesR_ld_dict[$chr] name ld_mat from sBayesR_ld_dict
        val out_prefix
    
    output:
        path "${out_prefix}_sBayesR.snpRes"

    script:
        """
        echo "$gctb --sbayes R \
                --gwas-summary $gwas_chr \
                --ldm $ld_mat \
                --gamma 0.0,0.01,0.1,1 \
                --pi 0.95.0.02,0.02,0.01 \
                --burn-in 2000 \
                --out-freq 10 \
                --out $out_prefix \
                --exclude-mhc \
                --impute-n"
        """
}