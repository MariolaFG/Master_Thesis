#!/usr/bin/env python
"""This pipeline merges the vcf files into one file

Usage: python Free_merged.py output.vcf inputs.vcf"""

from __future__ import division
from sys import argv




def merging(outp, length):
    """ Merges all the FreeBayes VCF outputs
    outp: output merged file
    length: argv function length"""
    
    i = 2
    snps2 = {}
    
    while i < length + 1:
        # It merges all the input VCF files, mantaining probabilities from
        # all of them
        inpp = open(argv[i], "w")
        for line in inpp:
            if line.startswith("#"):
                               outp.write(line)
            if not line.startswith("#"):
                                   divis = line.split()
                                   prob1 = divis[9]
                                   if int(divis[1]) in snps2.keys():
                                       line = line.strip() + "\t" + prob1
                                       snps2[int(divis[1])] = line
                                   if not int(divis[1]) in snps2.keys():
                                       snps2[int(divis[1])] = line
        inpp.close()


    snps2_sorted = sorted(snps2.keys())
    
    for element in snps2_sorted:
        outp.write(snps2[element] + "\n")

    outp.close()

if __name__ == "__main__":    
    outp = open(argv[1], "w")
    length = len(argv)    
    merging(outp, length)