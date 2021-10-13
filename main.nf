#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

include { split_reformat_gwas as split_for_prscs } from './modules/split_reformat_gwas.nf'
include { split_reformat_gwas as split_for_sbayesr } from './modules/split_reformat_gwas.nf'
include { split_reformat_gwas as split_for_prsice } from './modules/split_reformat_gwas.nf'
include { calc_posteriors_sbayesr as calc_posteriors_sbayesr_ukbb_eur_big } from './modules/calc_posteriors_sbayesr.nf'
include { calc_posteriors_sbayesr as calc_posteriors_sbayesr_ukbb_eur_hm3 } from './modules/calc_posteriors_sbayesr.nf'
include { calc_posteriors_prscs as calc_posteriors_prscs_ukbb_eur_hm3 } from './modules/calc_posteriors_prscs.nf'
include { calc_posteriors_prscs as calc_posteriors_prscs_1kg_eur_hm3 } from './modules/calc_posteriors_prscs.nf'
include { calc_score as calc_score_prscs_ukbb_eur_hm3 } from './modules/calc_score.nf'
include { calc_score as calc_score_prscs_1kg_eur_hm3 } from './modules/calc_score.nf'
include { calc_score as calc_score_sbayesr_ukbb_eur_big } from './modules/calc_score.nf'
include { calc_score as calc_score_sbayesr_ukbb_eur_hm3 } from './modules/calc_score.nf'
include { run_prsice } from './modules/run_prsice.nf'
include { eval_scores } from './modules/eval_scores.nf'

def help_msg() {
    log.info """
    A nextflow pipeline for performing polygenic score analysis at IBP
    Author: Vivek Appadurai | vivek.appadurai@regionh.dk

    Usage: nextflow run main.nf 
    
    Options

    --ref <gwas.vcf.gz> [A reference file of gwas summary stats in VCF.gz format] (Default: PGC SCZ 2014)
    --n <77096> [Reference GWAS Sample Size] (Default: 77096)
    --n_cases <Number of cases> [If the outcome is binary, case proportion is used to calculate liability transformed r2] (Default: 33640)
    --prevalence <Prevalence in decimals> [If the outcome is binary, population prevalence is used to calculate liability transformed r2] (Default: 0.1)
    --target <iPSYCH2012_All_Imputed_2021_QCed.json> [JSON file of target genotypes to score] (Default: iPSYCH2012 Imputed in 2021)
    --trait <ipsych_scz_2014> [A prefix for output files] (Default: Simple name of reference file)
    --covs <file.covs> [Path to covariates you might want to include in the NULL PGS model, such as age, gender, 10 PCs]
    --pheno <trait.pheno> [Path to the true phenotype file to evaluate the PGS performance]
    --binary <T/F> [Is the outcome binary or continuous?] (Default: T)
    --p_vals <comma separated list of p-values> [thresholds to be use for pruning and thresholding] (Default: 5e-8,1e-6,0.05,1)
    --help <prints this message>
    """
}

params.ref        = "/test/ieu-b42.vcf.gz"
params.n          = 77096
params.n_cases    = 33640
params.prevalence = 0.1
params.trait      = file(params.ref).getSimpleName()
params.covs       = ""
params.pheno      = ""
params.binary     = "T"
params.p_vals     = "5e-8,1e-6,0.05,1"
params.help       = false
target_prefix     = file(params.target).getSimpleName()
case_proportion   = params.n_cases/params.n 

if(params.help)
{
    help_msg()
    exit 0
}

log.info """
============================================================================================================
I B P - P R S -  P I P E L I N E _ v. 1.0 - N F
============================================================================================================
Reference GWAS               : $params.ref
Trait Name                   : $params.trait
Reference GWAS Sample Size   : $params.n
Target Genotypes             : $params.target
Target Prefix to use         : $target_prefix
PRSCS 1000G hapmap3 SNPs LD  : $params.prscs_1000G_hm3_eur_ld
PRSCS UKBB hapmap3 LD        : $params.prscs_ukbb_hm3_eur_ld
sBayesR UKBB hapmap3 SNPs LD : $params.sbayesr_ukbb_hm3_eur_ld
sBayesR UKBB 2.5M SNPs LD    : $params.sbayesr_ukbb_big_eur_ld
Output Directory             : $launchDir
Covariates                   : $params.covs
Phenotype                    : $params.pheno
p-value thresholds for P&T   : $params.p_vals
PLINK Path                   : $params.plink_path
Split GWAS Path              : $params.split_gwas_path
sBayesR Path                 : $params.sbayesr_path
PRS_CS Path                  : $params.prscs_path
PRSICE Path                  : $params.prsice_path
============================================================================================================
"""

//Parse the .json inputs for LD file paths, target PLINK genotypes

String sbayesr_ukbb_hm3_eur_ld   = new File(params.sbayesr_ukbb_hm3_eur_ld).text
String sbayesr_ukbb_big_eur_ld   = new File(params.sbayesr_ukbb_big_eur_ld).text
String prscs_1kg_hm3_eur_ld      = new File(params.prscs_1000G_hm3_eur_ld).text
String prscs_ukbb_hm3_eur_ld     = new File(params.prscs_ukbb_hm3_eur_ld).text
String genotypes                 = new File(params.target).text
def sbayesr_ukbb_hm3_eur_ld_dict = new JsonSlurper().parseText(sbayesr_ukbb_hm3_eur_ld)
def sbayesr_ukbb_big_eur_ld_dict = new JsonSlurper().parseText(sbayesr_ukbb_big_eur_ld)
def prscs_1kg_hm3_eur_ld_dict    = new JsonSlurper().parseText(prscs_1kg_hm3_eur_ld)
def prscs_ukbb_hm3_eur_ld_dict   = new JsonSlurper().parseText(prscs_ukbb_hm3_eur_ld)
def genotypes_dict               = new JsonSlurper().parseText(genotypes)

sbayesr_ukbb_hm3_eur_ld_ch = Channel.of(1..22) 
    | map {a -> [
        a,
        file(sbayesr_ukbb_hm3_eur_ld_dict[a.toString()]."bin").getBaseName(),
        sbayesr_ukbb_hm3_eur_ld_dict[a.toString()]."bin", 
        sbayesr_ukbb_hm3_eur_ld_dict[a.toString()]."info",
        sbayesr_ukbb_hm3_eur_ld_dict."meta"."cohort",
        sbayesr_ukbb_hm3_eur_ld_dict."meta"."population",
        sbayesr_ukbb_hm3_eur_ld_dict."meta"."snps",
        sbayesr_ukbb_hm3_eur_ld_dict."meta"."format"
        ]
    }

sbayesr_ukbb_big_eur_ld_ch = Channel.of(1..22) 
    | map {a -> [
        a,
        file(sbayesr_ukbb_big_eur_ld_dict[a.toString()]."bin").getBaseName(),
        sbayesr_ukbb_big_eur_ld_dict[a.toString()]."bin", 
        sbayesr_ukbb_big_eur_ld_dict[a.toString()]."info",
        sbayesr_ukbb_big_eur_ld_dict."meta"."cohort",
        sbayesr_ukbb_big_eur_ld_dict."meta"."population",
        sbayesr_ukbb_big_eur_ld_dict."meta"."snps",
        sbayesr_ukbb_big_eur_ld_dict."meta"."format"
        ]
    }

prscs_ukbb_hm3_eur_ld_ch = Channel.of(1..22) 
    | map {a -> [
        a,
        prscs_ukbb_hm3_eur_ld_dict[a.toString()]."bin", 
        prscs_ukbb_hm3_eur_ld_dict."info",
        prscs_ukbb_hm3_eur_ld_dict."meta"."cohort",
        prscs_ukbb_hm3_eur_ld_dict."meta"."population",
        prscs_ukbb_hm3_eur_ld_dict."meta"."snps",
        prscs_ukbb_hm3_eur_ld_dict."meta"."format"
        ]
    }

prscs_1kg_hm3_eur_ld_ch = Channel.of(1..22) 
    | map {a -> [
        a,
        prscs_1kg_hm3_eur_ld_dict[a.toString()]."bin", 
        prscs_1kg_hm3_eur_ld_dict."info",
        prscs_1kg_hm3_eur_ld_dict."meta"."cohort",
        prscs_1kg_hm3_eur_ld_dict."meta"."population",
        prscs_1kg_hm3_eur_ld_dict."meta"."snps",
        prscs_1kg_hm3_eur_ld_dict."meta"."format"
        ]
    }

ref_ch = Channel.of(1..22) 
    | map {a -> [a, params.ref, "${params.ref}.tbi"]}

genotypes_ch = Channel.of(1..22)
    | map {a -> [
        a,
        file(genotypes_dict[a.toString()]."bed").getBaseName(),
        genotypes_dict[a.toString()]."bed",
        genotypes_dict[a.toString()]."bim",
        genotypes_dict[a.toString()]."fam",
        genotypes_dict."meta"."cohort",
        genotypes_dict."meta"."population",
        genotypes_dict."meta"."snps"
        ]
    }

workflow {
    Channel.of(1..22) \
    | combine(Channel.of(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.of('prscs')) \
    | combine(Channel.of(params.split_gwas_path)) \
    | combine(Channel.of(n)) \
    | split_for_prscs \
    | set { prscs_input_ch } 

    Channel.of(1..22) \
    | combine(Channel.of(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.of('sbayesr')) \
    | combine(Channel.of(params.split_gwas_path)) \
    | combine(Channel.of(n)) \
    | split_for_sbayesr \
    | set { sbayesr_input_ch } 

    //Run PRS-CS

    Channel.of(1..22) \
    | combine(prscs_input_ch, by: 0) \
    | combine(Channel.of(params.n)) \
    | combine(prscs_ukbb_hm3_eur_ld_ch, by: 0) \
    | combine(genotypes_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of(params.prscs_path)) \
    | calc_posteriors_prscs_ukbb_eur_hm3 \
    | combine(Channel.of("2 4 6")) \
    | combine(Channel.of("${params.trait}_prscs_ukbb_eur_hm3")) \
    | combine(genotypes_ch, by: 0) \
    | combine(Channel.of(params.plink_path)) \
    | calc_score_prscs_ukbb_eur_hm3 \
    | collectFile(name: "${target_prefix}_${params.trait}_prscs_ukbb_eur_hm3.score",
        keepHeader: true,
        skip: 1) \
    | set { prscs_ukbb_eur_hm3_score_ch } 

    Channel.of(1..22) \
    | combine(prscs_input_ch, by: 0) \
    | combine(Channel.of(params.n)) \
    | combine(prscs_1kg_hm3_eur_ld_ch, by: 0) \
    | combine(genotypes_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of(params.prscs_path)) \
    | calc_posteriors_prscs_1kg_eur_hm3 \
    | combine(Channel.of("2 4 6")) \
    | combine(Channel.of("${params.trait}_prscs_1kg_hm3_eur")) \
    | combine(genotypes_ch, by: 0) \
    | combine(Channel.of(params.plink_path)) \
    | calc_score_prscs_1kg_eur_hm3 \
    | collectFile(name: "${target_prefix}_${params.trait}_prscs_1kg_eur_hm3.score",
        keepHeader: true,
        skip: 1) \
    | set { prscs_1kg_eur_hm3_score_ch }

    //Run SBayesR

    Channel.of(1..22) \
    | combine(sbayesr_input_ch, by: 0) \
    | combine(sbayesr_ukbb_big_eur_ld_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of(params.sbayesr_path)) \
    | calc_posteriors_sbayesr_ukbb_eur_big \
    | combine(Channel.of("2 5 8")) \
    | combine(Channel.of("${params.trait}_sbayesr_ukbb_eur_2.8M")) \
    | combine(genotypes_ch, by: 0) \
    | combine(Channel.of(params.plink_path)) \
    | calc_score_sbayesr_ukbb_eur_big \
    | collectFile(name: "${target_prefix}_${params.trait}_sbayesr_ukbb_eur_2.8M.score", 
        keepHeader: true,
        skip: 1) \
    | set { sbayesr_ukbb_eur_big_score_ch }

    Channel.of(1..22) \
    | combine(sbayesr_input_ch, by: 0) \
    | combine(sbayesr_ukbb_hm3_eur_ld_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of(params.sbayesr_path)) \
    | calc_posteriors_sbayesr_ukbb_eur_hm3 \
    | combine(Channel.of("2 5 8")) \
    | combine(Channel.of("${params.trait}_sbayesr_ukbb_eur_hm3")) \
    | combine(genotypes_ch, by: 0) \
    | combine(Channel.of(params.plink_path)) \
    | calc_score_sbayesr_ukbb_eur_hm3 \
    | collectFile(name: "${target_prefix}_${params.trait}_sbayesr_ukbb_eur_hm3.score", 
        keepHeader: true,
        skip: 1) \
    | set { sbayesr_ukbb_eur_hm3_score_ch }

    // Run PRSice

    Channel.of(1..22) \
    | combine(Channel.of(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.of('prsice')) \
    | combine(Channel.of(params.split_gwas_path)) \
    | combine(Channel.of(n)) \
    | split_for_prsice \
    | combine(genotypes_ch, by: 0) \
    | combine(Channel.of(params.trait)) \
    | combine(Channel.of(params.p_vals)) \
    | combine(Channel.of(params.binary)) \
    | combine(Channel.of(params.pheno)) \
    | combine(Channel.of(params.prsice_path)) \
    | run_prsice \
    | collectFile(name: "${target_prefix}_${params.trait}_prsice.score", 
        keepHeader: true,
        skip: 1) \
    | set { prsice_score_ch }

    // Evaluate the different PGS methods

    prsice_score_ch \
    | combine(sbayesr_ukbb_eur_big_score_ch) \ 
    | combine(sbayesr_ukbb_eur_hm3_score_ch) \
    | combine(prscs_ukbb_eur_hm3_score_ch) \
    | combine(prscs_1kg_eur_hm3_score_ch) \ 
    | combine(Channel.of($params.pheno)) \
    | combine(Channel.of($params.covs)) \
    | combine(Channel.of($params.binary)) \
    | combine(Channel.of($params.prevalence)) \
    | combine(Channel.of(case_proportion)) \
    | combine(Channel.of("${target_prefix}_${params.trait}")) \
    | combine(Channel.of($params.eval_scores_path)) \
    | eval_scores
} 