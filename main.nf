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
    --prscs_ld_files <"ldblk_1kg_chr*.hdf5"> [Paths to chromosome wise LD matrices in hdf5 format] (Default: "ldblk_1kg_chr*.hdf5")
    --sbayesR_ld_files <"ukbEURu_hm3_chr*_v3_50k.ldm.sparse.*"> [Paths to chromosome wise LD matrices in .bin/.info format] (Default: "ukbEURu_hm3_chr*_v3_50k.ldm.sparse.*")
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
============================================================================================================
I B P - P R S -  P I P E L I N E _ v. 1.0 - N F
============================================================================================================
Reference GWAS              : $params.ref
Trait Name                  : $params.trait
Reference GWAS Sample Size  : $params.N
Target dataset              : $params.bfile
Prs-CS LD Files             : $params.prscs_ld_files
sBayesR-LD-Files            : $params.sbayesr_ld_files
Output Directory            : $params.dir
============================================================================================================
"""

def jsonSlurper_prscs   = new JsonSlurper()
def jsonSlurper_sBayesR = new JsonSlurper()
def prscs_ld_json       = new File(params.prscs_ld_files)
def sBayesR_ld_json     = new File(params.sbayesr_ld_files)
String prscs_ld_paths   = prscs_ld_json.text
String sBayesR_ld_paths = sBayesR_ld_json.text 
def prscs_ld_dict       = jsonSlurper_prscs.parseText(prscs_ld_paths) 
def sBayesR_ld_dict     = jsonSlurper_sBayesR.parseText(sBayesR_ld_paths)

// Channel.fromFilePairs("${ref_vcf_gz.getParent()}/${ref_vcf_gz.getBaseName()}.{gz,gz.tbi}", size : 2, checkIfExists : true).set{ref_ch}
prscs_ld_ch = Channel.of(1..22) | map {a -> [a, prscs_ld_dict.get(a.toString())]}
sbayesr_ld_ch = Channel.of(1..22) | map {a -> [a, sBayesR_ld_dict.get(a.toString())]}
Channel.fromFilePairs("${params.bfile}.{bed,bim,fam}", size : 3, checkIfExists : true)
    { file -> file.getSimpleName().replaceAll(/.+chr/,'') }
    .set{geno_ch}
ref_ch = Channel.of(1..22) | map {a -> [a, file("${params.ref}"), file("${params.ref}.tbi")]}

workflow {
    Channel.of(1..22) \
    | combine(Channel.from(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.from(params.N)) \
    | combine(Channel.from('prscs')) \
    | split_for_prscs \
    | combine(Channel.of(params.N)) \
    | combine(prscs_ld_ch, by: 0) \
    | combine(geno_ch, by: 0) \
    | combine(Channel.of(params.dir)) \
    | combine(Channel.of(params.trait)) \
    | calc_posteriors_prscs.collect() \
    | combine(Channel.of('prscs')) \
    | combine(Channel.of(params.trait)) \
    | merge_prscs

    Channel.of(1..22) \
    | combine(Channel.from(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.from(params.N)) \
    | combine(Channel.from('sbayesr')) \
    | split_for_sbayesr \
    | combine(sbayesr_ld_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | calc_posteriors_sbayesr.collect() \
    | combine(Channel.of('sbayesr'))  \
    | combine(Channel.of(params.trait)) \ 
    | merge_sbayesr
} 