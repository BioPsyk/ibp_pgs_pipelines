#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    input:
        tuple val(chr),
            path(gwas),
            val(N),
            path(ld_file),        
            path(bed), 
            path(bim), 
            path(fam), 
            path(out_dir),
            val(traitName)
    
    output:
        path "${out_prefix}_prscs_chr${chr}.snpRes"

    script:
        """
        echo "python ${projectDir}/bin/PRScs.py --ref_dir=${ld.file.getParent()}/${ld.file.getBaseName()} \
                --sst_file=$gwas \
                --bim_prefix=${bim.getParent()}/${bim.getBaseName()} \
                --n_gwas=$N \
                --chrom=$chr \
                --out_dir=$out_dir/$traitName"
        """ 
}