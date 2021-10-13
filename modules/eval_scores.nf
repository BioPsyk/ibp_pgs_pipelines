#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process eval_scores {
    label 'mod_mem'
    publishDir launchDir
    
    input:
        tuple path(prsice_scores),
            path(sbayesr_ukbb_big_scores),
            path(sbayesr_ukbb_hm3_scores),
            path(prscs_ukbb_hm3_scores),
            path(prscs_1kg_hm3_scores),
            path(pheno_file),
            path(covs_file),
            val(binary),
            val(prevalence),
            val(case_proportion),
            val(out_prefix),
            path(eval_scores_script)

    output:
        path "${out_prefix}_VarianceExplained.txt"
        path "${out_prefix}_VarianceExplained.png"
        path "${out_prefix}_Scores.txt"
        
    script:
    if(binary == "T")
        """
        Rscript ./evaluate_pgs.R $prscs_scores \
            $sbayesr_ukbb_big_scores \
            $sbayesr_ukbb_hm3_scores \
            $prscs_ukbb_hm3_score \
            $prscs_1kg_hm3_score \
            $out_prefix \
            $pheno_file \
            $covs_file \
            --binary \
            --prevalence $prevalence \
            --case_pct $case_proportion
        """
    else
        """
        Rscript ./evaluate_pgs.R $prscs_scores \
            $sbayesr_ukbb_big_scores \
            $sbayesr_ukbb_hm3_scores \
            $prscs_ukbb_hm3_scores \
            $prscs_1kg_hm3_score \
            $out_prefix \
            $pheno_file \
            $covs_file
        """
}