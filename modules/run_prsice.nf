#!/usr/bin/env nextflow

process run_prsice {
    publishDir launchDir
    label 'big_mem'

    input: 
        tuple path(source),
            path(bed),
            path(bim),
            path(fam),
            val(trait),
            val(p_value_thresholds),
            val(binary),
            path(pheno),
            path(prsice),
    output:
        "$trait.all_score"

    script:
    """
        ./PRsice_linux \
            --base $source \
            --target $bed.getBaseName() \
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