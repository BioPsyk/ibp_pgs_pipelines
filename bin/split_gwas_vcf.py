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
    parser.add_argument('--vcf', type = str,help = "Path to a gzipped VCF", required = True)
    parser.add_argument('--format', type = str, help = "prscs or sbayesr or prsice", required = True)
    parser.add_argument('--out',  type = str, help = "Output prefix", required = True)

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
            else:
                out_fh.write(variant.CHROM)
                out_fh.write("_")
                out_fh.write(str(variant.start + 1)) # Positions returned by CyVCF2 API are 0-based
                out_fh.write("_")
                out_fh.write(ref_allele) 
                out_fh.write("_")
                out_fh.write(alt_allele)
                out_fh.write(" ")
            if(args.format == 'prsice'):
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
            out_fh.write(str(variant.format('ES').flat[0]))
            out_fh.write(" ")
            if (args.format == "sbayesr" or args.format == "prsice"):
                out_fh.write(str(variant.format('SE').flat[0]))
                out_fh.write(" ")
            out_fh.write(str(pow(10, -1 * variant.format('LP').flat[0])))
            out_fh.write(" ")
            if (args.format == "sbayesr"):
                out_fh.write(str(variant.format('SS').flat[0]))
            out_fh.write("\n")
    else:
        sys.exit("FATAL: File not found: " + args.vcf + "!!")

    out_fh.close()

if __name__ == '__main__':
    main()