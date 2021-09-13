#!/usr/bin/env nextflow

process run_prsice {
    publishDir launchDir
    label 'big_mem'

    input: 
        tuple path(source),
            val(bfile),
            path(bed),
            path(bim),
            path(fam),
            val(trait),
            val(p_vals),
            val(binary),
            path(pheno),
            path(prsice)
    output:
        path "${trait}.all_score"

    script:
    """
        ./PRSice_linux \
            --base $source \
            --target $bfile \
            --thread 8 \
            --binary-target $binary \
            --pheno $pheno \
            --out $trait \
            --all-score \
            --fastscore \
            --bar-levels $p_vals \
            --no-regress
    """
}