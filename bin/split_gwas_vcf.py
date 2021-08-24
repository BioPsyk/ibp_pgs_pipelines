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
    parser.add_argument('--format', type = str, help = "prscs or sbayesR", required = True)

    args        = parser.parse_args()
    out         = re.sub('^.*\/', '', args.vcf)
    out         = re.sub('\.vcf.gz$', '', out)
    args.format = args.format.lower()
    out         = out + "_" + args.format + "_" + "chr" + str(args.chromosome) + ".txt"
    out_fh      = open(out, "w")

    if(args.format == "prscs"):
        out_fh.write("SNP A1 A2 BETA P\n")
    elif(args.format == "sbayesr"):
        out_fh.write("SNP A1 A2 freq b se p N\n")
    else:
        sys.exit("FATAL: --format is either prscs or sbayesR")

    if path.exists(args.vcf):
        gwas_vcf = VCF(args.vcf)
        for variant in gwas_vcf(args.chromosome):
            if (variant.ID):
                out_fh.write(variant.ID)
                out_fh.write(" ")
            else:
                out_fh.write(variant.CHROM)
                out_fh.write("_")
                out_fh.write(str(variant.start))
                out_fh.write("_")
                out_fh.write(variant.REF) 
                out_fh.write("_")
                out_fh.write(''.join(variant.ALT))
                out_fh.write(" ")
            out_fh.write(''.join(variant.ALT))
            out_fh.write(" ")
            out_fh.write(variant.REF)
            out_fh.write(" ")
            if(args.format == "sbayesr"):
                out_fh.write(str(variant.INFO.get('AF')))
                out_fh.write(" ")
            out_fh.write(str(variant.format('ES').flat[0]))
            out_fh.write(" ")
            if (args.format == "sbayesr"):
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