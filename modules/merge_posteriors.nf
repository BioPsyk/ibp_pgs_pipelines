#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

merge_posteriors {
    input:
        path posteriors
        val method
    output:
        path ${method}_snpPosteriorEffects.txt
    script:
        """
        echo "
        head -n 1 $posteriors[0] > ${method}_snpPosteriorEffects.txt
        for i in {1..22};do cat $posteriors[\$i] >> ${method}_snpPosteriorEffects.txt"
        """
}