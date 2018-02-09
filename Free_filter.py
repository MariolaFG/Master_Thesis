#!/usr/bin/env python
"""This pipeline filters freebayes output to select true allelic SNPs

Usage: python Free_filter.py free_total.vcf free_filter.vcf"""

from __future__ import division
from sys import argv

def filtering(inpp, outp):
    """ Filters out homeologous SNPs
    inpp: input VCF file
    outp: output VCF file """
    
    for line in inpp:
    	if line.startswith("#"):
    		outp.write(line)
    	if not line.startswith("#"):
    		divis = line.split()
    		prob1 = divis[9:] #Contains a list with each genotype's probabilites
    		if len(prob1) < 2:
            # If the list of probabilities is smaller than the number of 
            # genotypes (two in this case) it is considered a homologous SNP
    			outp.write(line)
    
    
    inpp.close()
    outp.close()

if __name__ == "__main__":
    inpp = open(argv[1], "r")
    outp = open(argv[2], "w")
    
    filtering(inpp, outp)
