#!/usr/bin/env Rscript

require(dplyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(fmsb, quietly = TRUE)
require(argparser, quietly = TRUE)

parser = arg_parser("Evaluate performance of Polygenic Scores", 
                    hide.opts = TRUE)
parser = add_argument(parser, 
                      "prsice", 
                      help = "Path  to PRSice scores")
parser = add_argument(parser, 
                      "sbayesr_ukbb_big", 
                      help = "Path to sBayesR scores with 2.5M UKBB SNP set")
parser = add_argument(parser, 
                      "sbayesr_ukbb_hm3", 
                      help = "Path to sBayesR scores with Hapmap3 SNP set")
parser = add_argument(parser, 
                      "prscs_ukbb_hm3", 
                      help = "Path to PRS-CS  scores with Hapmap3 SNP set")
parser = add_argument(parser, 
                      "prscs_1kg_hm3", 
                      help = "Path to PRS-CS scores with 1000 Genomes SNP set")
parser = add_argument(parser, 
                      "out",
                      help = "Output Prefix",
                      default = "PGS_Eval")
parser = add_argument(parser,
                      "pheno",
                      help = "Phenotype file to evaluate, expected cols: FID IID Pheno")
parser = add_argument(parser,
                      "covar",
                      help = "File of covariates to use in tge PGS null model, expected cols: FID IID...")
parser = add_argument(parser,
                      "--binary",
                      flag = TRUE,
                      help = "Flag if the outcome is binary")
parser = add_argument(parser,
                      "--prevalence",
                      help = "Prevalence in case of a binary trait for liability transformation")
parser = add_argument(parser,
                      "--case_pct",
                      help = "Case proportion in case of a binary trait for liability transformation")
options = parse_args(parser)

if(isTRUE(options$binary) && (!exists(options$prevalence) || !exists(options$case_pct))) {
    print("ERROR: --binary flag implies binary trait and requires --prevalence & --case_pct!")
    stop()
}

scaled_unique_scores = function (x_df) {
    colnames(x_df) = c("FID", "IID", "N_MISS_ALLELE_CT", 
                       "NAMED_ALLELE_DOSAGE_SUM", 
                       "SCORE_AVG", "SCORE_SUM")
    x_df = x_df %>% 
        group_by(IID) %>%
        mutate(score_genome_sum = sum(SCORE_SUM)) %>%
        ungroup() %>%
        as.data.frame() %>%
        mutate(SCORE_SCALED = scale(score_sum)) %>%
        select(IID, SCORE_SCALED) %>% 
        unique()
    return (x_df)
}

pheno = read.table(options$pheno, header = TRUE)
pheno = pheno %>% select(-FID)
covar = read.table(options$covar, header = TRUE)
covar = covar %>% select(-FID)

# Read PRSice scores, assumes there are scores at 4 different p-value thresholds
# 5e-08, 1e-06, 0.05, 1 and the scores are summed but not averaged

prsice_scores = read.table(options$prsice, header = TRUE)
colnames(prsice_scores) = c("FID", "IID", 
                            "PT_5E8", "PT_1E6", "PT_0.05", "PT_1")
prsice_scores = prsice_scores %>% 
    group_by(IID) %>% 
    mutate(PT_5E8_sum = sum(PT_5E8)) %>% 
    mutate(PT_1E6_sum = sum(PT_1E6)) %>% 
    mutate(PT_0.05_sum = sum(PT_0.05)) %>%
    mutate(PT_1_sum = sum(PT_1)) %>%
    ungroup() %>%
    as.data.frame() %>%
    mutate(PT_5E8_scaled = scale(PT_5E8_sum)) %>%
    mutate(PT_1E6_scaled = scale(PT_1E6_sum)) %>%
    mutate(PT_0.05_scaled = scale(PT_0.05_sum)) %>%
    mutate(PT_1_scaled = scale(PT_1_sum)) %>%
    select(-PT_5E8, -PT_1E6, -PT_0.05, -PT_1, 
           -PT_1eE8_sum, -PT_1E6_sum, -PT_0.05_sum, -PT_1_sum, -FID) %>%
    unique()

# Read sBayesR, PRS-CS scores, assumes the scores are summed but not averaged

sbayesr_ukbb_big_scores = read.table(options$sbayesr_ukbb_big, header = TRUE)
sbayesr_ukbb_big_scores = sbayesr_ukbb_big_scores %>% scaled_unique_scores()
colnames(sbayesr_ukbb_big_scores) = c("IID", "sBayesR_UKBB_2.5M")

sbayesr_ukbb_hm3_scores = read.table(options$sbayesr_ukbb_hm3, header = TRUE)
sbayesr_ukbb_hm3_scores = sbayesr_ukbb_hm3_scores %>% scaled_unique_scores()
colnames(sbayesr_ukbb_hm3_scores) = c("IID", "sBayesR_UKBB_Hapmap3")

prscs_ukbb_hm3_scores   = read.table(options$prscs_ukbb_hm3, header = TRUE)
prscs_ukbb_hm3_scores   = prscs_ukbb_hm3_scores %>% scaled_unique_scores()
colnames(prscs_ukbb_hm3_scores) = c("IID", "PRSCS_UKBB_Hapmap3")

prscs_1kg_hm3_scores    = read.table(options$prscs_1kg_hm3, header = TRUE)
prscs_1kg_hm3_scores    = prscs_1kg_hm3_scores %>% scaled_unique_scores()
colnames(prscs_1kg_hm3_scores) = c("IID", "PRSCS_1KG_Hapmap3")

pgs = inner_join(prsice_scores, sbayesr_ukbb_hm3_scores, by = c("IID"))
pgs = inner_join(pgs, sbayesr_ukbb_big_scores, by = c("IID"))
pgs = inner_join(pgs, prscs_ukbb_hm3_scores, by = c("IID"))
pgs = inner_join(pgs, prscs_1kg_hm3_scores, by = c("IID"))
pgs = inner_join(pgs, pheno, by = c("IID"))
pgs = inner_join(pgs, covar, by = c("IID"))

null_model = glm(data = pgs, pheno ~ . -PT_5E8_scaled, -PT_1E6_scaled,
                 -PT_0.05_scaled, -PT_1_scaled, 
                 -sBayesR_UKBB_2.5M, -sBayesR_UKBB_Hapmap3, 
                 -PRSCS_UKBB_Hapmap3, -PRSCS_1KG_Hapmap3, -IID)