#!/usr/bin/env Rscript

# This script fills the missing columns in prsice output
# Missing columns occur if there are no SNPs in GWAS passing a certain threshold

require(dplyr, quietly = TRUE)
require(data.table, quietly = TRUE)

args       = commandArgs(trailingOnly = TRUE)
score      = fread(args[1])
out_prefix = args[2]

if(!"Pt_5e.08" %in% colnames(score)) {
    score = score %>% mutate(Pt_5e.08 = 0)
}

if(!"Pt_1e.06" %in% colnames(score)) {
    score = score %>% mutate(Pt_1e.06 = 0)
}

if(!"Pt_0.05" %in% colnames(score)) {
    score = score %>% mutate(Pt_0.05 = 0)
}

if(!"Pt_1" %in% colnames(score)) {
    score = score %>% mutate(Pt_1 = 0)
}

score = score %>% select(FID, IID, Pt_5e.08, Pt_1e.06, Pt_0.05, Pt_1)

write.table(score, 
            paste0(out_prefix, "_", "colCheck.txt"),
            row.names = F, 
            quote = F, 
            sep = " ")