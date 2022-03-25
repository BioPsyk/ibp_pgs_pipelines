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
    tab delimited format format for prs-cs or prsice""")
    parser.add_argument('--chromosome', 
        type = str, 
        help = "Chromosome to split from VCF", 
        required = True)
    parser.add_argument('--vcf', 
        type = str, 
        help = "Path to a gzipped VCF", 
        required = True)
    parser.add_argument('--format', 
        type = str, 
        help = "prscs (or) sbayesr (or) prsice (or) ldpred", 
        required = True)
    parser.add_argument('--out',  
        type = str, 
        help = "Output prefix", 
        required = True)
    parser.add_argument('--n', 
        type = int, 
        help = "GWAS Sample Size", 
        required = True)

    args = parser.parse_args()
    args.format = args.format.lower()
    out = '_'.join([args.out, args.format, "chr"])
    out = ''.join([out, str(args.chromosome)])
    out = ".".join([out, "txt"])
    out_fh = open(out, "w")
    snp_dict = {}

    if(args.format == "prscs"):
        out_fh.write("SNP A1 A2 BETA P\n")
    elif(args.format == "sbayesr"):
        out_fh.write("SNP A1 A2 freq b se p N\n")
    elif(args.format == "prsice"):
        out_fh.write("SNP CHR BP A1 A2 BETA SE P\n")
    elif(args.format == "ldpred"):
        out_fh.write("rsid,chr,pos,a0,a1,beta,beta_se,N,p\n")
    else:
        sys.exit("ERROR: Unsupported --format specified: ", args.format)

    if path.exists(args.vcf):
        gwas_vcf = VCF(args.vcf)
        for variant in gwas_vcf(args.chromosome):
            chromosome = str(variant.CHROM)
            position = str(variant.start + 1)
            otherAllele = variant.REF
            effectAllele = ''.join(variant.ALT)

            if (variant.ID and variant.ID != "."):
                snp = variant.ID
                if snp in snp_dict.keys():
                    print("WARNING: DUPLICATE VARIANT : ", snp + ", Skipping..")
                    continue
            else:
                snp = '_'.join([chromosome, position, effectAllele, otherAllele])
            snp_dict[snp] = 1
            try:
                effect = str(variant.format('ES').flat[0])
            except:
                print("WARNING: Effect size not available at SNP: ", snp, ", Skipping..")
                continue
            try:
                std_error = str(variant.format('SE').flat[0])
            except:
                print("WARNING: Standard Error of effect not available at SNP: ", snp, ", Skipping..")
                continue
            try:
                p_value = str(variant.format('EP').flat[0])
            except:
                print("WARNING: p-value not available at SNP: ", snp, ", Skipping..")
                continue
            try:
                effectAlleleFreq = str(variant.format('AFKG').flat[0])
            except:
                print("WARNING: Allele Frequency missing at SNP: ", snp, ", Skipping..")
                continue
            N = str(args.n)

            if(len(effectAllele) > 1 or len(otherAllele) > 1): 
                print("WARNING: INDEL/MULTI ALLELIC VARIANT at SNP:", snp, ", Skipping..")
                continue

            if(args.format == "prscs"):
                out_line = ' '.join([snp, effectAllele, otherAllele, effect, p_value])
            elif(args.format == "sbayesr"):
                out_line = ' '.join([snp, effectAllele, otherAllele, effectAlleleFreq, effect, std_error, p_value, N])
            elif(args.format == "prsice"):
                out_line = ' '.join([snp, chromosome, position, effectAllele, otherAllele, effect, std_error, p_value])
            elif(args.format == "ldpred"):
                out_line = ','.join([snp, chromosome, position, otherAllele, effectAllele, effect, std_error, N, p_value])
            out_fh.write(out_line + "\n")
    else:
        sys.exit("ERROR: File not found: " + args.vcf + "!!")
    out_fh.close()

if __name__ == '__main__':
    main()