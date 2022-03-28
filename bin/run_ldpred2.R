#!/usr/bin/env Rscript

require(data.table, quietly = TRUE)
require(magrittr, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(bigsnpr, quietly = TRUE)
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
                      "out",
                      help = "Output Prefix",
                      default = "PGS_Eval")
options = parse_args(parser)

sumstats = bigreadr::fread2(options$sumstats)

# Calculate LD matrix from genotypes

# Get maximum amount of cores
NCORES = nb_cores()
# Open a temporary file
tmp = tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

# Initialize variables for storing the LD score and LD matrix

corr       = NULL
ld         = NULL
info_snp   = NULL
fam.order  = NULL
obj.bigSNP = snp_attach(options$genotypes)
map        = obj.bigSNP$map
names(map) = c("chr", "rsid", "dist", "pos", "a1", "a0")
    
# Intersect sumstats with genotype map
    
info_snp  = snp_match(sumstats, map)

# Assign the genotype to a variable for easier downstream analysis
    
genotype = obj.bigSNP$genotypes
    
# Rename the data structures

CHR  = map$chr
POS  = map$pos
POS2 = map$dist

for (chr in 1:22) {
    #Extract SNPs in chromosome
    
    ind.chr = which(info_snp$chr == chr)
    ind.chr2 = info_snp$`_NUM_ID_`[ind.chr]
    
    # Calculate LD
    
    corr0 = snp_cor(
        genotype,
        ind.col = ind.chr2,
        ncores = NCORES,
        info.pos = POS2[ind.chr2],
        size = 3/1000
    )
    
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } 
    else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }   
}

fam.order <- as.data.table(obj.bigSNP$fam)

# Rename fam order

setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

# Perform LD Score regression 

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(ld, 
                 length(ld), 
                 chi2 = (df_beta$beta / df_beta$beta_se)^2,
                 sample_size = df_beta$n_eff, 
                 blocks = NULL)
h2_est <- ldsc[["h2"]]

# Obtain posteriors - INFINITESIMAL MODEL

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)

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
                               vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
                               alllow_jump_sign = FALSE,
                               shrink_corr = 0.95,
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto)
    auto$beta_est)

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

# ---------------- END OF CODE --------------#
