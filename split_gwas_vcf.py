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
    parser.add_argument('--chromosome', type = str, help = "Chromosome to split from VCF")
    parser.add_argument('--vcf', type = str,help = "Path to a gzipped VCF")

    args           = parser.parse_args()
    out            = re.sub('\.vcf.gz$', '', args.vcf)
    sbayesr_out    = out + "_" + "sBayesR" + "_" + str(args.chromosome) + ".txt"
    prscs_out      = out + "_" + "prscs" + "_" + str(args.chromosome) + ".txt"
    sbayesr_out_fh = open(sbayesr_out, "w")
    prscs_out_fh   = open(prscs_out, "w")

    sbayesr_out_fh.write("SNP A1 A2 freq b se p N\n")
    prscs_out_fh.write("SNP A1 A2 BETA P\n")
 
    if path.exists(args.vcf):
        gwas_vcf = VCF(args.vcf)
        for variant in gwas_vcf(args.chromosome):
            sbayesr_out_fh.write(variant.ID)
            sbayesr_out_fh.write(" ")
            prscs_out_fh.write(variant.ID)
            prscs_out_fh.write(" ")
            sbayesr_out_fh.write(''.join(variant.ALT))
            sbayesr_out_fh.write(" ")
            prscs_out_fh.write(''.join(variant.ALT))
            prscs_out_fh.write(" ")
            sbayesr_out_fh.write(variant.REF)
            sbayesr_out_fh.write(" ")
            prscs_out_fh.write(variant.REF)
            prscs_out_fh.write(" ")
            sbayesr_out_fh.write(str(variant.INFO.get('AF')))
            sbayesr_out_fh.write(" ")
            sbayesr_out_fh.write(str(variant.format('ES').flat[0]))
            sbayesr_out_fh.write(" ")
            prscs_out_fh.write(str(variant.format('ES').flat[0]))
            prscs_out_fh.write(" ")
            sbayesr_out_fh.write(str(variant.format('SE').flat[0]))
            sbayesr_out_fh.write(" ")
            sbayesr_out_fh.write(str(pow(10, -1 * variant.format('LP').flat[0])))
            sbayesr_out_fh.write(" ")
            prscs_out_fh.write(str(pow(10, -1 * variant.format('LP').flat[0])))
            prscs_out_fh.write("\n")
            sbayesr_out_fh.write(str(variant.format('SS').flat[0]))
            sbayesr_out_fh.write("\n")
    else:
        sys.exit("FATAL: File not found: " + args.vcf + "!!")

    sbayesr_out_fh.close()
    prscs_out_fh.close()

if __name__ == '__main__':
    main()