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
            path(prsice),
            val(binary),
            path(pheno),
            val(p_vals)
    output:

    script:
    """
        Rscript $prsice --dir . \
            --prsice ./PRsice \
            --base $source \
            --target $bed.getBaseName() \
            --thread 8 \
            --stat BETA \
            --beta \
            --binary-target $binary \
            --snp SNP \
            --chr CHR \
            --bp BP \
            --A1 A1 \
            --A2 A2 \
            --stat BETA \
            --pvalue P \
            --pheno $pheno \
            --pheno-col 3 \
            --out $trait \
            --all-score \
            --fastscore \
            --bar-levels $p_vals \
            --no-regress
    """
}