#!/usr/bin/env python

import sys
from cyvcf2 import VCF
import argparse
from os import path
import re
import numpy


def main():
    parser = argparse.ArgumentParser(description = """Splits a summary stats file in VCF 
    format by chromosome.
    Reformats to:
    mt-cojo format for sBayesR,
    tab delimited format format for prs-cs""")
    parser.add_argument('--chromosome', type = str, help = "Chromosome to split from VCF", required = True)
    parser.add_argument('--vcf', type = str, help = "Path to a gzipped VCF", required = True)
    parser.add_argument('--format', type = str, help = "prscs | sbayesr | prsice | ldpred2", required = True)
    parser.add_argument('--out',  type = str, help = "Output prefix", required = True)
    parser.add_argument('--n',  type = int, help = "GWAS Sample Size", required = True)

    args        = parser.parse_args()
    args.format = args.format.lower()
    out         = args.out + "_" + args.format + "_" + "chr" + str(args.chromosome) + ".txt"
    out_fh      = open(out, "w")
    snp_dict    = {}

    if(args.format == "prscs"):
        out_fh.write("SNP A1 A2 BETA P\n")
    elif(args.format == "sbayesr"):
        out_fh.write("SNP A1 A2 freq b se p N\n")
    elif(args.format == "prsice"):
        out_fh.write("SNP CHR BP A1 A2 BETA SE P\n")
    elif(args.format == "ldpred2"):
        # a1 is effect allele and a0 is non-effect allele for LDPred2
        out_fh.write("rsid chr pos a1 a0 beta beta_se p N\n")
    else:
        sys.exit("FATAL: --format is invalid")

    if path.exists(args.vcf):
        gwas_vcf = VCF(args.vcf)
        for variant in gwas_vcf(args.chromosome):
            alt_allele = ''.join(variant.ALT)
            ref_allele = variant.REF
            if(len(alt_allele) > 1 or len(ref_allele) > 1): # Ignores INDELs in the input GWAS dataset
                continue
            if (variant.ID):
                if variant.ID in snp_dict.keys():
                    continue
                else:
                    out_fh.write(variant.ID)
                    out_fh.write(" ")
                    snp_dict[variant.ID] = 1
            else: # Create a new ID in chr_pos_a1_a2 format, if variant ID does not exist in input VCF 
                out_fh.write(variant.CHROM)
                out_fh.write("_")
                out_fh.write(str(variant.start + 1)) # Positions returned by CyVCF2 API are 0-based
                out_fh.write("_")
                out_fh.write(ref_allele) 
                out_fh.write("_")
                out_fh.write(alt_allele)
                out_fh.write(" ")
            if(args.format == "prsice" or args.format == "ldpred2"):
                out_fh.write(variant.CHROM)
                out_fh.write(" ")
                out_fh.write(str(variant.start + 1))
                out_fh.write(" ")
            out_fh.write(alt_allele)
            out_fh.write(" ")
            out_fh.write(ref_allele)
            out_fh.write(" ")
            if(args.format == "sbayesr"):
                out_fh.write(str(variant.INFO.get('AF')))
                out_fh.write(" ")
            out_fh.write(str(variant.format('ES').flat[0])) # Expects effects to be BETAs
            out_fh.write(" ")
            if (args.format != "prscs"):
                out_fh.write(str(variant.format('SE').flat[0]))
                out_fh.write(" ")
            out_fh.write(str(pow(10, -1 * variant.format('LP').flat[0])))
            out_fh.write(" ")
            if (args.format == "sbayesr" or args.format == "ldpred2"):
                out_fh.write(str(args.n))
            out_fh.write("\n")
    else:
        sys.exit("FATAL: File not found: " + args.vcf + "!!")

    out_fh.close()

if __name__ == '__main__':
    main()