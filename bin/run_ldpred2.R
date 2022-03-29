#!/usr/bin/env Rscript

require(data.table, quietly = TRUE)
require(magrittr, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(bigsnpr, quietly = TRUE)
require(bigreadr, quietly = TRUE)
require(argparser, quietly = TRUE)

parser = arg_parser("Run LDPred2 to calculate Polygenic Scores", 
                    hide.opts = TRUE)
parser = add_argument(parser, 
                      "sumstats", 
                      help = "Sumstats in LdPred2 format")
parser = add_argument(parser, 
                      "genotypes", 
                      help = "Genotype file in BigSnpR format")
parser = add_argument(parser,
                      "map",
                      help = "map .RDS file for LD reference")
parser = add_argument(parser,
                      "ldm",
                      help = "LDpred2 Sparse LD matrix .RDS file")
parser = add_argument(parser, 
                      "out",
                      help = "Output Prefix",
                      default = "PGS_Eval")
parser = add_argument(parser,
                      "chromosome",
                      help = "Chromosome to analyze")
options = parse_args(parser)

# Reading input summary stats

sumstats = bigreadr::fread2(options$sumstats)
sumstats = sumstats %>% rename(n_eff = N)

# Get maximum amount of cores

NCORES = 1 # change to nb_cores() when R pkg bigparallelr is fixed for bugs

# Information for variants provided in the LD reference

map_ldref = readRDS(options$map)
map_ldref = map_ldref %>% filter(chr == options$chromosome)
info_snp  = snp_match(sumstats, map_ldref)
info_snp  = tidyr::drop_na(tibble::as_tibble(info_snp))
sd_ldref  = with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss     = with(info_snp, 2 / sqrt(n_eff * beta_se^2))
is_bad    = sd_ss < (0.5 * sd_ldref) | 
    sd_ss > (sd_ldref + 0.1) | 
    sd_ss < 0.1 | 
    sd_ldref < 0.05
df_beta = info_snp[!is_bad, ]

tmp      = tempfile(tmpdir = "tmp-data")
ind_chr  = which(df_beta$chr == options$chromosome)
ind.chr2 = df_beta$`_NUM_ID_`[ind.chr]
ind.chr3 = match(ind.chr2, which(map_ldref$chr == chr))
corr_chr = readRDS(options$ldm)[ind.chr3, ind.chr3]
corr     = as_SFBM(corr_chr, tmp)

# Heritability estimation of LD score regression to be used as a starting value 
# in LDpred2-auto

ldsc = with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref), 
                              chi2 = (beta / beta_se)^2,
                              sample_size = n_eff,
                              ncores = NCORES))
h2_est = ldsc[["h2"]]

# Obtain posteriors - INFINITESIMAL MODEL

beta_inf = snp_ldpred2_inf(corr, df_beta, h2 = h2_est)

# GRID MODEL

h2_seq = round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq  = signif(seq_log(1e-5, 1, length.out = 21), 2)

grid.param = expand.grid(p = p_seq, 
                         h2 = h2_seq, 
                         sparse = c(FALSE, TRUE))

beta_grid = snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)

# AUTO MODEL

multi_auto <- snp_ldpred2_auto(corr, 
                               df_beta,
                               h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, length.out = 30),
                               allow_jump_sign = FALSE,
                               shrink_corr = 0.95,
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

# LASSOSUM2 GRID

beta_lassosum2 = snp_lassosum2(corr, df_beta, ncores = NCORES)

# Write SNP posteriors to file

write.table(beta_inf, 
            paste0(options$out, "_LdPred2_Beta_Infinitesimal.txt"), 
            sep = " ", 
            row.names = F, 
            quote = F)

write.table(beta_grid, 
            paste0(options$out, "_LdPred2_Beta_Grid.txt"), 
            sep = " ", 
            row.names = F, 
            quote = F)

write.table(beta_auto, 
            paste0(options$out, "_LdPred2_Beta_Auto.txt"), 
            sep = " ", 
            row.names = F, 
            quote = F)

write.table(beta_lassosum2, 
            paste0(options$out, "_LassuSum2_Beta_Grid.txt"), 
            sep = " ", 
            row.names = F, 
            quote = F)

# cleanup

file.remove(paste0(tmp, ".sbk"))

# ---------------- END OF CODE --------------#
