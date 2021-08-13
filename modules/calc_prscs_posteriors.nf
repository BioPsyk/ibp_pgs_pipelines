#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// calculate per chromosome posterior SNP effects for PRS-CS

process prs_cs {
    input:
        val chr, path gwas_chr
        path ld_mat
        val N
        path out_dir
        val bfile, path plink_files

    output:
        val chr, path posteriors into prscs_posterior_chunks

    script:
        """
        echo "python /bin/PRScs.py --ref_dir=$ld_mat \
                                --sst_file=$gwas_chr \
                                --bim_prefix=$bfile \
                                --n_gwas=$N \
                                --chrom=$chr \
                                --out_dir=$out_dir"
                                """ 
}