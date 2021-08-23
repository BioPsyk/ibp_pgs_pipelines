#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process merge_posteriors {
    publishDir "$out_dir/$method/posteriors/"

    input:
        tuple path(posteriors),
            val(method),
            val(trait),
            path(out_dir)

    output:
        path("${trait}_${method}_snpPosteriorEffects.txt")

    script:
        """
        echo "
        head -n 1 $posteriors[0] > ${trait}_${method}_snpPosteriorEffects.txt
        for i in {1..22};do cat $posteriors[\$i] >> ${trait}_${method}_snpPosteriorEffects.txt"
        """
}