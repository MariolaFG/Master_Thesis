#!/usr/bin/env python
"""This pipeline selects the SNPs detected by Bcftools call that are homozygous 
to the reference allele in at least one genotype
Usage: python Not_Anchors.py prefilter.vcf filter.vcf"""

from __future__ import division
from sys import argv

inpp = open(argv[1], "r")
outp = open(argv[2], "w")

snps_dic = {}
for line in inpp:

    if line.startswith("#"):
    	outp.write(line)
    	# This will write on the output file all the information contained on the Bcftools call file
    else:
        divis = line.split()
        position = divis[1]
        probs = divis[9:]
        pos = []

    for element in probs:
        pos.append(element.split(":")[1].split(","))



    i = 0
    a = range(0, (len(pos)-1))

    for element in pos:
        if element[0] == "0" and element[2] != "0":
            snps_dic[int(position)] = line



for snp in sorted(snps_dic.keys()):
    outp.write(snps_dic[snp] + "\n")


inpp.close()
outp.close()

