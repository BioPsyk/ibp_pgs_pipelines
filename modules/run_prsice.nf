#!/usr/bin/env nextflow

process run_prsice {
    label 'big_mem'

    input: 
        tuple val(chr),
            path(source),
            val(bfile),
            path(bed),
            path(bim),
            path(fam),
            val(target_cohort),
            val(target_population),
            val(target_snp_set),
            val(trait),
            val(p_vals),
            val(binary),
            path(pheno),
            path(prsice),
            path(col_check_script)
    
    output:
        path "${trait}_chr${chr}_colCheck.txt"

    script:
    """
        ./PRSice_linux \
            --base $source \
            --target $bfile \
            --thread 8 \
            --binary-target $binary \
            --pheno $pheno \
            --out ${trait}_chr${chr} \
            --all-score \
            --fastscore \
            --bar-levels $p_vals \
            --no-regress \
            --score sum

        Rscript ./fill_missing_cols_prsice.R ${trait}_chr${chr}.all_score ${trait}_${chr}
    """
}