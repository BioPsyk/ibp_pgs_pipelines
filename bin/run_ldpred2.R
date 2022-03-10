#!/usr/bin/env Rscript

require(data.table)
require(magrittr)
require(dplyr)
require(ggplot2)
require(bigsnpr)

args = commandArgs(trailingOnly = TRUE)
phenotype = args[1]
covariates = args[2]