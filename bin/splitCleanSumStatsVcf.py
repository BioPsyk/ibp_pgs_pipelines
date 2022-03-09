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
    tab delimited format format for prs-c or prsice""")
    parser.add_argument('--chromosome', type = str, help = "Chromosome to split from VCF", required = True)
    parser.add_argument('--vcf', type = str, help = "Path to a gzipped VCF", required = True)
    parser.add_argument('--format', type = str, help = "prscs or sbayesr or prsice", required = True)
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
    else:
        sys.exit("ERROR: Unsupported --format specified. Select either prscs or sbayesr or prsice!")

    if path.exists(args.vcf):
        gwas_vcf = VCF(args.vcf)
        for variant in gwas_vcf(args.chromosome):
            chromosome = variant.chromosome
            position   = variant.start + 1
            ref_allele = variant.REF
            alt_allele = ''.join(variant.ALT)
            effect     = str(variant.format('ES').flat[0])
            std_error  = str(variant.format('SE').flat[0])
            p_value    = str(variant.format('EP').flat[0])
            maf        = variant.format('AF').flat[0]
            N          = str(args.N)

            if(len(alt_allele) > 1 or len(ref_allele) > 1): # Ignores INDELs in the input GWAS dataset
                print("INFO: Ignoring INDEL/MULTI ALLELIC VARIANT at :" + snp + ":" + chromosome + " " + position)
                continue

            if (maf > 0.5):
                maf = 1 - 0.5

            maf = str(maf)

            if (variant.ID):
                snp = variant.ID
                if snp in snp_dict.keys():
                    print("INFO: Ignoring DUPLICATE VARIANT :" + snp + "\n")
                    continue
            else:
                snp = chromosome + "_" + position + "_" + ref_allele + "_" + alt_allele

            out_fh.write(snp + " ")
            if(args.format == 'prsice'):
                out_fh.write(chromosome + " ")
                out_fh.write(position + " ")
            out_fh.write(alt_allele + " ")
            out_fh.write(ref_allele + " ")
            if(args.format == "sbayesr"):
                out_fh.write(maf + " ")
            out_fh.write(effect + " ")
            if (args.format == "sbayesr" or args.format == "prsice"):
                out_fh.write(std_error + " ")
            out_fh.write(p_value + " ")
            if (args.format == "sbayesr"):
                out_fh.write(N)
            out_fh.write("\n")
    else:
        sys.exit("ERROR: File not found: " + args.vcf + "!!")

    out_fh.close()

if __name__ == '__main__':
    main()