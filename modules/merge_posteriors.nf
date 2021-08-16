#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process merge_posteriors {
    publishDir "${out_dir}/${method}/posteriors/", mode: "copy"

    input:
        path posteriors
        val method
        path out_dir

    output:
        path ${method}_snpPosteriorEffects.txt

    script:
        """
        echo "
        head -n 1 $posteriors[0] > ${method}_snpPosteriorEffects.txt
        for i in {1..22};do cat $posteriors[\$i] >> ${method}_snpPosteriorEffects.txt"
        """
}