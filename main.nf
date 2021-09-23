#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

include { split_reformat_gwas as split_for_prscs } from './modules/split_reformat_gwas.nf'
include { split_reformat_gwas as split_for_sbayesr } from './modules/split_reformat_gwas.nf'
include { split_reformat_gwas as split_for_prsice } from './modules/split_reformat_gwas.nf'
include { calc_posteriors_sbayesr } from './modules/calc_posteriors_sbayesr.nf'
include { calc_posteriors_prscs } from './modules/calc_posteriors_prscs.nf'
include { calc_score as calc_score_prscs} from './modules/calc_score.nf'
include { calc_score as calc_score_sbayesr } from './modules/calc_score.nf'
include { run_prsice as run_prsice} from './modules/run_prsice.nf'
//include { eval_scores as eval_prs }

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
    --prscs_ld <prscs_eur_ld.json> [Paths to chromosome wise LD matrices in hdf5 format] (Default: "ldblk_1kg_chr*.hdf5")
    --sbayesR_ld <"sbayesr_eur_ld.json> [Paths to chromosome wise LD matrices in .bin/.info format] (Default: "ukbEURu_hm3_chr*_v3_50k.ldm.sparse.*")
    --covs <file.covs> [Path to covariates you might want to include in the NULL PGS model, such as age, gender, 10 PCs]
    --pheno <trait.pheno> [Path to the true phenotype file to evaluate the PGS performance]
    --binary <T/F> [Is the outcome binary or continuous?] (Default: T)
    --p_value_thresholds <comma separated list of p-values> [p-value thresholds to be used for pruning and thresholding method] (Default: 5e-8,1e-6,0.05,1)
    --help <prints this message>
    """
}

params.ref                = "/test/ieu-b42.vcf.gz"
params.bfile              = "/test/Genomes1k_Phase1_CEU_chr1_22"
params.N                  = "77096"
params.dir                = launchDir
params.trait              = file(params.ref).getSimpleName()
params.prscs_ld           = ""
params.sbayesr_ld         = "sbayesR_eur_ld.json"
params.help               = false
params.covs               = ""
params.pheno              = ""
params.binary             = "T"
params.p_value_thresholds = "5e-8,1e-6,0.05,1"
split_gwas_path           = "$projectDir/bin/split_gwas_vcf.py"
prscs_path                = "$projectDir/bin/PRScs/PRScs.py"
sbayesr_path              = "$projectDir/bin/gctb_2.03beta_Linux/gctb"
plink_path                = "$projectDir/bin/plink2"
prsice_path               = "$projectDir/bin/PRSice_linux"

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
Prs-CS LD Directory         : $params.prscs_ld
sBayesR LD File Paths       : $params.sbayesr_ld
Output Directory            : $params.dir
PGS Covariates              : $params.covs
Phenotype                   : $params.pheno
p-value thresholds for P&T  : $params.p_value_thresholds
============================================================================================================
"""

def jsonSlurper_sBayesR = new JsonSlurper()
def jsonSlurper_prscs   = new JsonSlurper()
def sBayesR_ld_json     = new File(params.sbayesr_ld)
def prscs_ld_json       = new File(params.prscs_ld)
String sBayesR_ld_paths = sBayesR_ld_json.text 
String prscs_ld_paths   = prscs_ld_json.text
def sBayesR_ld_dict     = jsonSlurper_sBayesR.parseText(sBayesR_ld_paths)
def prscs_ld_dict       = jsonSlurper_prscs.parseText(prscs_ld_paths)

sbayesr_ld_ch = Channel.of(1..22) 
    | map {a -> [a, file(sBayesR_ld_dict.get(a.toString())["bin"]).getBaseName(), 
        sBayesR_ld_dict.get(a.toString())["bin"], 
        sBayesR_ld_dict.get(a.toString())["info"]]}
prscs_ld_ch = Channel.of(1..22) 
    | map {a -> [a, prscs_ld_dict.get(a.toString())["bin"],
        prscs_ld_dict.get(a.toString())["info"],
        prscs_ld_dict.get(a.toString())["pop"]]}
ref_ch = Channel.of(1..22) 
    | map {a -> [a, params.ref, "${params.ref}.tbi"]}
plink_ch = Channel.of(1..22) 
    | map {a -> [a, file("${params.bfile}.bed").getSimpleName(), 
        "${params.bfile}.bed", 
        "${params.bfile}.bim", 
        "${params.bfile}.fam"]}

workflow {

    //Run PRS-CS

    Channel.of(1..22) \
    | combine(Channel.of(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.of(params.N)) \
    | combine(Channel.of('prscs')) \
    | combine(Channel.of(split_gwas_path)) \
    | split_for_prscs \
    | combine(Channel.of(params.N)) \
    | combine(prscs_ld_ch, by: 0) \
    | combine(plink_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of(prscs_path)) \
    | calc_posteriors_prscs \
    | combine(Channel.of("2 4 6")) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of('prscs')) \
    | combine(Channel.of(params.bfile)) \
    | combine(Channel.of(plink_ch, by: 0)) \
    | combine(Channel.of(plink_path)) \
    | calc_score_prscs \
    | collectFile(name: "${params.trait}_prscs.sscore", 
        keepHeader: true,
        skip: 1)

    //Run SBayesR

    Channel.of(1..22) \
    | combine(Channel.of(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.of(params.N)) \
    | combine(Channel.of('sbayesr')) \
    | combine(Channel.of(split_gwas_path)) \
    | split_for_sbayesr \
    | combine(sbayesr_ld_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of(sbayesr_path)) \
    | calc_posteriors_sbayesr \
    | combine(Channel.of("2 5 8")) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of('sbayesr')) \
    | combine(plink_ch, by: 0) \
    | combine(Channel.of(plink_path)) \
    | calc_score_sbayesr \
    | collectFile(name: "${params.trait}_sbayesr.sscore", 
        keepHeader: true,
        skip: 1)

    // Run PRSice

    Channel.of(1..22) \
    | combine(Channel.of(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.of(params.N)) \
    | combine(Channel.of('prsice')) \
    | combine(Channel.of(split_gwas_path)) \
    | split_for_prsice \
    | combine(plink_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of(params.p_value_thresholds)) \
    | combine(Channel.of(params.binary)) \
    | combine(Channel.of(params.pheno)) \
    | combine(Channel.of(prsice_path)) \
    | run_prsice \
    | collectFile(name: "${params.trait}_prsice.all_score")

    //eval_prs(calc_score_prscs.out, calc_score_sbayesr.out, $params.covs, $params.trait, $params.pheno) 
} 