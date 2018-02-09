#!/usr/bin/env python

"""This pipeline selects the SNPs detected by Bcftools call that are homozygous 
to the reference allele in at least one genotype
Usage: python Not_Anchors.py prefilter.vcf filter.vcf"""

from __future__ import division
from sys import argv

inpp1 = open(argv[1], "r")
outp = open(argv[2], "w")


for line in inpp1:
	if line.startswith("#"):
                outp.write(line)
		# This will write on the output file all the information contained on the Bcftools call file
	if not line.startswith("#"):
		divis = line.split()
		prob1 = divis[9].split(":")[1].split(",")
		prob2 = divis[10].split(":")[1].split(",")

		if prob1[2] != "0" and prob2[2] != "0":
			if prob1[0] == "0":
				if prob2[0] != "0":
					outp.write(line)

			if prob2[0] == "0":
                        	if prob1[0] != "0":
                                	outp.write(line)

inpp1.close()
outp.close()
