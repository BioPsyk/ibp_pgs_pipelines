#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process eval_scores {
    input:
        path prscs_scores,
        path sbayesr_scores,
        val traitName
    output:
        path "${traitName}.VarianceExplained.txt"
    script:
        """
        eval_prs.R ${prscs_scores} ${sbayesr_scores} ${traitName}
        """
}