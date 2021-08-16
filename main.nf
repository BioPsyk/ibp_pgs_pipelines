#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

include { split_reformat_gwas as split_for_prscs } from './modules/split_reformat_gwas.nf'
include { split_reformat_gwas as split_for_sbayesr } from './modules/split_reformat_gwas.nf'
include { calc_posteriors_sbayesr } from './modules/calc_posteriors_sbayesr.nf'
include { calc_posteriors_prscs } from './modules/calc_posteriors_prscs.nf'
include { merge_posteriors as merge_sbayesr } from './modules/merge_posteriors.nf'
include { merge_posteriors as merge_prscs } from './modules/merge_posteriors.nf'

def help_msg() {
    log.info """
    A nextflow pipeline for performing polygenic score analysis at IBP
    Author: Vivek Appadurai | vivek.appadurai@regionh.dk

    Usage: nextflow run main.nf 
    
    Options

    --ref <gwas.vcf.gz> [A reference file of gwas summary stats in VCF.gz format] (Default: PGC SCZ 2014)
    --trait <ipsych_scz_2014> [A prefix for output files, preferably containing name of trait and dataset] (Default: Simple name of reference file)
    --N <10000> [Sample size of the reference dataset] (Default: 77096)
    --dir </path/to/results/> [Place where the results from the pipeline are to be stored] (Default: Launch directory)
    --prscs_ld_files <prscs_eur_ld.json> [A .json format file with paths to chromosome wise LD matrices in hdf5 format] (Default: prscs_eur_ld.json)
    --sbayesR_ld_files <sbayesR_eur_ld.json> [A .json format file with paths to chromosome wise LD matrices in sparse bin format] (Default: sbayesR_eur_ld.json)
    --help <prints this message>
    """
}

params.ref              = "/test/ieu-b42.vcf.gz"
params.bfile            = "/test/Genomes1k_Phase1_CEU_chr1_22"
params.N                = "77096"
params.dir              = launchDir
params.trait            = file(params.ref).getSimpleName()
params.prscs_ld_files   = "prscs_eur_ld.json"
params.sbayesr_ld_files = "sbayesR_eur_ld.json"
params.help             = false

if(params.help)
{
    help_msg()
    exit 0
}

log.info """
======================================================
I B P - P R S -  P I P E L I N E _ v. 1.0 - N F
======================================================
Reference GWAS              : $params.ref
Trait Name                  : $params.trait
Reference GWAS Sample Size  : $params.N
Target dataset              : $params.bfile
Prs-CS LD Files             : $params.prscs_ld_files
sBayesR-LD-Files            : $params.sbayesr_ld_files
Output Directory            : $params.dir
======================================================
"""

def jsonSlurper_prscs   = new JsonSlurper()
def jsonSlurper_sBayesR = new JsonSlurper()
def prscs_ld_json       = new File(params.prscs_ld_files)
def sBayesR_ld_json     = new File(params.sbayesr_ld_files)
String prscs_ld_files   = prscs_ld_json.text
String sBayesR_ld_files = sBayesR_ld_json.text 
def prscs_ld_dict       = jsonSlurper_prscs.parseText(prscs_ld_files) 
def sBayesR_ld_dict     = jsonSlurper_sBayesR.parseText(sBayesR_ld_files)

workflow {
    split_for_prscs(Channel.of(1..22), 
        params.trait, 
        params.ref, 
        params.N,
        "prscs")

    split_for_sbayesr(Channel.of(1..22), 
        params.trait, 
        params.ref, 
        params.N,
        "sbayesr")

    calc_posteriors_prscs(split_for_prscs.out[0], 
        split_for_prscs.out[1],
        params.N,  
        Channel.fromPath(prscs_ld_dict.get("${split_for_prscs.out[0]}")),
        Channel.fromFilePairs("${params.bfile}.{bed, bim, fam}", checkIfExists: true))

    calc_posteriors_sbayesr(split_for_sbayesr.out[0],
        split_for_sbayesr.out[1],
        Channel.fromFilePairs("${sBayesR_ld_dict.get("${split_for_sbayesr.out[0]}").get("bin").getBaseName}.{bin, info}", checkIfExists: true),
        params.trait)

    merge_prscs(calc_posteriors_prscs.out.collect(), "prscs", params.dir)
    merge_sbayesr(calc_posteriors_sbayesr.out.collect(), "sBayesR", params.dir)
}