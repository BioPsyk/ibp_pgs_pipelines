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
                      help = "Path to sBayesR scores with 2.8M UKBB SNP set")
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
                      help = "Phenotype file to evaluate, 
                      expected cols: FID IID Pheno")
parser = add_argument(parser,
                      "covar",
                      help = "File of covariates to use in the PGS null model, 
                      expected cols: FID IID...")
parser = add_argument(parser,
                      "--binary",
                      flag = TRUE,
                      help = "Flag if the outcome is binary")
parser = add_argument(parser,
                      "--prevalence",
                      help = "Prevalence in case of a binary trait
                      for liability transformation")
options = parse_args(parser)

if(isTRUE(options$binary) & !exists(options$prevalence)) {
    print("ERROR: --binary flag implies binary trait and requires --prevalence")
    stop()
}

##################### Functions used in the script #############################

# Sum PGS across chromosomes and scale to zero mean, unit variance

sum_scale_scores = function (x_df, score_col_name) {
    colnames(x_df) = c("IID", "SCORE_CHR")
    x_df = x_df %>% 
        group_by(IID) %>%
        mutate(SCORE_SUM = sum(SCORE_CHR)) %>%
        ungroup() %>%
        as.data.frame() %>%
        mutate(SCORE_SCALED = scale(SCORE_SUM)) %>%
        select(IID, SCORE_SCALED) %>% 
        unique()
    colnames(x_df) = c("IID", score_col_name)
    return (x_df)
}

# Calculate r2 and p-value of PGS

calculate_r2_p = function(x_df, y_df, binary) {
    x_df = semi_join(x_df, y_df, by = c("IID"))
    pgs  = inner_join(x_df, y_df, by = c("IID"))
    r2   = 0
    p    = 0
    
    if(isTRUE(binary)) {
        fit_null = glm(data = x_df, Pheno ~ . -IID)
        fit_bin  = glm(data = pgs, Pheno ~ . -IID)
        r2       = NagelkerkeR2(fit_bin) - NagelkerkeR2(fit_null)
        p        = pchisq(deviance(fit_null) - deviance(fit_bin),
                          df.residual(fit_null) - df.residual(fit_bin),
                          lower.tail = F) 
    } else {
        fit_null = lm(data = x_df, Pheno ~ . -IID)
        fit_con  = lm(data = pgs, Pheno ~ . -IID)
        r2       = summary(fit_con)$r.squared - summary(fit_null)$r.squared
        p        = pchisq(deviance(fit_null) - deviance(fit_con),
                          df.residual(fit_null) - df.residual(fit_con),
                          lower.tail = F) 
    }
    return (c(r2, p))
}

# Transform variance explained to liability scale for binary traits

liability_transform = function(r2, k, p) {
    x     = qnorm(1 - k)
    z     = dnorm(x)
    i     = z / k
    cc    = k * (1 - k) * k * (1 - k) / (z * z * p * (1 - p))
    theta = i * ((p - k) / (1 - k)) * (i * ((p - k) / (1- k)) - x)
    e     = 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
    r2_L  = cc * e * r2 / (1 + cc * e * theta * r2)
    return (r2_L)
}

################################################################################

pheno           = read.table(options$pheno, header = TRUE)
colnames(pheno) = c("FID", "IID", "Pheno")
pheno           = pheno %>% select(-FID)
covar           = read.table(options$covar, header = TRUE)
covar           = covar %>% select(-FID)

# Read PRSice scores, assumes there are scores at 4 different p-value thresholds
# 5e-08, 1e-06, 0.05, 1 and that the scores are summed but not averaged

prsice_scores = read.table(options$prsice, header = TRUE)
colnames(prsice_scores) = c("FID", "IID", 
                            "PT_5E8", "PT_1E6", "PT_0.05", "PT_1")
prsice_5E8_scores  = prsice_scores %>% 
    select(IID, PT_5E8) %>% 
    sum_scale_scores("PT_5E8")
prsice_1E6_scores  = prsice_scores %>% 
    select(IID, PT_1E6) %>% 
    sum_scale_scores("PT_1E6")
prsice_0.05_scores = prsice_scores %>% 
    select(IID, PT_0.05) %>% 
    sum_scale_scores("PT_0.05")
prsice_1_scores = prsice_scores %>% 
    select(IID, PT_1) %>% 
    sum_scale_scores("PT_1")

# Read sBayesR, PRS-CS scores, assumes the scores are summed but not averaged

sbayesr_ukbb_2.8m_scores = read.table(options$sbayesr_ukbb_big, header = TRUE)
sbayesr_ukbb_2.8m_scores = sbayesr_ukbb_2.8m_scores %>% 
    select(IID, SCORE1_SUM) %>% 
    sum_scale_scores("sBayesR_UKBB_2.8M")

sbayesr_ukbb_hm3_scores = read.table(options$sbayesr_ukbb_hm3, header = TRUE)
sbayesr_ukbb_hm3_scores = sbayesr_ukbb_hm3_scores %>% 
    select(IID, SCORE1_SUM) %>%
    sum_scale_scores("sBayesR_UKBB_HM3")

prscs_ukbb_hm3_scores = read.table(options$prscs_ukbb_hm3, header = TRUE)
prscs_ukbb_hm3_scores = prscs_ukbb_hm3_scores %>% 
    select(IID, SCORE1_SUM) %>% 
    sum_scale_scores("PRSCS_UKBB_HM3")

prscs_1kg_hm3_scores = read.table(options$prscs_1kg_hm3, header = TRUE)
prscs_1kg_hm3_scores = prscs_1kg_hm3_scores %>% 
    select(IID, SCORE1_SUM) %>% 
    sum_scale_scores("PRSCS_1KG_HM3")

pheno_cov = inner_join(pheno, covar, by = c("IID"))

# Calculate Nagelkerke R2 for binary traits

prsice_5E8_eval = calculate_r2_p(pheno_cov,
                                 prsice_5E8_scores,
                                 options$binary)
prsice_1E6_eval = calculate_r2_p(pheno_cov,
                                 prsice_1E6_scores,
                                 options$binary)
prsice_0.05_eval = calculate_r2_p(pheno_cov,
                                  prsice_0.05_scores,
                                  options$binary)
prsice_1_eval = calculate_r2_p(pheno_cov,
                               prsice_1_scores,
                               options$binary)
sbayesr_ukbb_2.8m_eval = calculate_r2_p(pheno_cov,
                                        sbayesr_ukbb_2.8m_scores,
                                        options$binary)
sbayesr_ukbb_hm3_eval = calculate_r2_p(pheno_cov,
                                       sbayesr_ukbb_hm3_scores,
                                       options$binary)
prscs_1kg_hm3_eval = calculate_r2_p(pheno_cov,
                                    prscs_1kg_hm3_scores,
                                    options$binary)
prscs_ukbb_hm3_eval = calculate_r2_p(pheno_cov,
                                     prscs_ukbb_hm3_scores,
                                     options$binary)

pgs = inner_join(prsice_5E8_scores, prsice_1E6_scores, by = c("IID"))
pgs = inner_join(pgs, prsice_0.05_scores, by = c("IID"))
pgs = inner_join(pgs, prsice_1_scores, by = c("IID"))
pgs = inner_join(pgs, sbayesr_ukbb_2.8m_scores, by = c("IID"))
pgs = inner_join(pgs, sbayesr_ukbb_hm3_scores, by = c("IID"))
pgs = inner_join(pgs, prscs_ukbb_hm3_scores, by = c("IID"))
pgs = inner_join(pgs, prscs_1kg_hm3_scores, by = c("IID"))

all_methods_eval = calculate_r2_p(pheno_cov,
                                  pgs,
                                  options$binary)

r2_out = data.frame(Method = c("PRsice_5E8",
                               "PRSice_1E6",
                               "PRSice_0.05",
                               "PRSice_1",
                               "sBayesR_UKBB_2.8M",
                               "sBayesR_UKBB_HM3",
                               "PRSCS_UKBB_HM3",
                               "PRSCS_1KG_HM3",
                               "ALL"),
                    r2 = c(prsice_5E8_eval[1],
                           prsice_1E6_eval[1],
                           prsice_0.05_eval[1],
                           prsice_1_eval[1],
                           sbayesr_ukbb_2.8m_eval[1],
                           sbayesr_ukbb_hm3_eval[1],
                           prscs_ukbb_hm3_eval[1],
                           prscs_1kg_hm3_eval[1],
                           all_methods_eval[1]),
                    P = c(prsice_5E8_eval[2],
                          prsice_1E6_eval[2],
                          prsice_0.05_eval[2],
                          prsice_1_eval[2],
                          sbayesr_ukbb_2.8m_eval[2],
                          sbayesr_ukbb_hm3_eval[2],
                          prscs_ukbb_hm3_eval[2],
                          prscs_1kg_hm3_eval[2],
                          all_methods_eval[2]))

if(isTRUE(options$binary)) {
    n_cases = sum(pheno$Pheno)
    n_samples = nrow(pheno$Pheno)
    if(n_cases > n_samples) {
        n_cases = n_samples - n_cases
    }
    case_pct = n_cases/n_samples

    prsice_5E8_r2_L        = liability_transform(prsice_5E8_eval[1], 
                                                 options$case_pct, 
                                                 options$prevalence)
    prsice_1E6_r2_L        = liability_transform(prsice_1E6_eval[1], 
                                                 options$case_pct, 
                                                 options$prevalence)
    prsice_0.05_r2_L       = liability_transform(prsice_0.05_eval[1], 
                                                 options$case_pct, 
                                                 options$prevalence)
    prsice_1_r2_L          = liability_transform(prsice_1_eval[1], 
                                                 options$case_pct, 
                                                 options$prevalence)
    sBayesR_ukbb_2.8m_r2_L = liability_transform(sbayesr_ukbb_2.8m_eval[1], 
                                                 options$case_pct, 
                                                 options$prevalence)
    sBayesR_ukbb_hm3_r2_L  = liability_transform(sbayesr_ukbb_hm3_eval[1], 
                                                 options$case_pct, 
                                                 options$prevalence)
    prscs_ukbb_hm3_r2_L    = liability_transform(prscs_ukbb_hm3_eval[1], 
                                                 options$case_pct, 
                                                 options$prevalence)
    prscs_1kg_hm3_r2_L     = liability_transform(prscs_1kg_hm3_eval[1], 
                                                 options$case_pct, 
                                                 options$prevalence)
    
    r2_L = data.frame("r2_L"= c(prsice_5E8_r2_L, 
                                prsice_1E6_r2_L,
                                prsice_0.05_r2_L,
                                prsice_1_r2_L,
                                sBayesR_ukbb_2.8m_r2_L,
                                sBayesR_ukbb_hm3_r2_L,
                                prscs_ukbb_hm3_r2_L,
                                prscs_1kg_hm3_r2_L))
    r2_out = cbind(r2_out, r2_L)
    
    # Make plot
    
    png(paste0(options$out, "_", "VarianceExplainedLiabilityScale.png"), 
        width = 5, 
        height = 5, 
        units = "in", 
        res = 300)
    
    ggplot(r2_out, aes(x = Method, y = r2_L)) + 
        geom_bar(stat = "identity") + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab("")
    
    dev.off()
}

#------------------------------ Make plots -------------------------------------

png(paste0(options$out, "_", "VarianceExplained.png"), 
    width = 5, 
    height = 5, 
    units = "in", 
    res = 300)

ggplot(r2_out, aes(x = Method, y = r2)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("")

dev.off()

#--------------------------- Write outputs -------------------------------------

write.table(r2_out, paste0(options$out, "_", "VarianceExplained.txt"), 
                           row.names = F, 
                           quote = F, 
                           sep = "\t")

write.table(pgs, paste0(options$out, "_", "Scores.txt"), 
            row.names = F, 
            quote = F, 
            sep = "\t")

######################### --- End of Script --- ################################