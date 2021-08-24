#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process score {
    input: 
        tuple val(chr),
            path(posteriors),
            path(bed),
            path(bim),
            path(fam),
            val(trait),
            val(method),
            val(out_dir)
    output:
        val(chr),
        path("${trait}_${method}_chr${chr}.profile")
    script:
        echo "plink --bfile $bed.getBaseName() \
                --score $posteriors 2 4 header sum \
                --out "${trait}_${method}_chr${chr}"
}