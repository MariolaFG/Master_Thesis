#!/usr/bin/env python

"""This pipeline writes in a txt file the performance of the filtering tool detecting
SNPs
Usage: python Compare.py HaploGenerator_varianthaplos.txt(separated by commas) 
filtering_tool_output.vcf reference_genome Number_of_Organisms"""


from __future__ import division
from sys import argv
import subprocess


def reference(ref, genome, organisms):
	""" Creates a dictionary containing the SNPs generated by HaploGenerator

	ref = HaploGenerator varianthaplos.txt files
	genome = reference genome
	organisms = number of organisms"""

	ref_dic = {}
	ref_true = {}
	homo = {}
	total_list = []
	ref = ref.split(",")

	for reference in ref:
		reference = open(reference)
		for line in reference:
			line = line.strip().split("\t")
			if line[1] == "Chromosome_name":
				if line[2] in ref_dic.keys():
					ref_dic[line[2]].append(line[3:])
				if not line[2] in ref_dic.keys():
					ref_dic[line[2]] = [line[3:]]
	# For every HaploGenerator file, create a dictionary: Position = information.

	for element in ref_dic:
		if len(ref_dic[element]) == int(organisms):
			homo[element] = ref_dic[element][0]
			total_list.append(int(element))

		if len(ref_dic[element]) < int(organisms):
			if ref_dic[element][0].count("1") == ((int(organisms)-1) * 2):
				ref_true[element] = [ref_dic[element][0], "homo"]
				total_list.append(int(element))
			if ref_dic[element][0].count("1") != 0 and ref_dic[element][0].count("1") < int(organisms):
                                ref_true[element] = [ref_dic[element][0], "hetero"]
				total_list.append(int(element))
	# Checks every position in ref_dic and if it is present in all organisms is stored in the homeologous SNPs dictionary
	# if it is not present in all organisms it checks if it is a homozygous SNP or a heterozygous SNP and it adds it to
	# the true allelic dictionary

	genome = genome.read().replace(">Chromosome_name", "").replace("\n", "")
	true_per = len(ref_true)/len(genome)*100
	false_per = len(homo)/len(genome)*100
	# This two parameters represent the percentage of true allelic SNPs and the percentage of homeologous SNPs in the
	# genome

	return ref_true, homo, total_list, true_per, false_per

def filter_out(tool):
	"""Creates a dictionary containing the SNPs detected by the filtering tool

        tool = filter tool vcf output containing the SNP"""

	tool_dic = {}
	for line in tool:
    		if line.startswith("Chromosome_name"):
        		line = line.split("\t")
        		if line[1] in tool_dic:
            			tool_dic[line[1]].append(line[2:])
        		elif line[1] not in gatk_dic:
            			tool_dic[line[1]] = line[2:]

	return tool_dic

def true_positives(tool_dic, ref_true, homo):
	"""Calculates the true allelic SNPs and homeologous SNPs detected by the filtering tool

	tool_dic = dictionary containing the SNPs detected by the filtering tool
	ref_true = dictionary containing the true allelic SNPs generated by HaploGenerator
	homo = dictionary containing the homeologous SNPs generated by HaploGenerator"""

	snps = {}
	snps_called = ""
	snps_called_tool = ""
	homo_snps = []
	if len(tool_dic) != 0:

		for key in tool_dic:
    			for ref_key in ref_dic:
        			if key == ref_key:
					snps[str(key)] = ref_dic[key][1]

		homo_snps = []
		for key in tool_dic:
                	for ref_key in homo:
                       		if key == ref_key:
					if not key in homo_snps:
                               			homo_snps.append(key)

	# Checks every SNP detected by the filtering tool and and stores them in the dictionary of homeologous SNPs
	# or the dictionary of true allelic SNPs detected.

	if len(ref_dic) != 0:
		snps_called = "True positives: " + str(len(snps)/len(ref_dic)*100)[0:5] + " %"
	if len(ref_dic) == 0:
		snps_called = "No True SNPs in Genome"
	# Calculates the true allelic SNPs detected percentage

	return snps, snps_called, homo_snps

def false_positives(tool_dic, snps):
	"""Calculates the false positives detected by the filtering tool

        tool_dic = dictionary containing the SNPs detected by the filtering tool
        snps =  list of true allelic SNPs detected by the filtering tool """

	if len(tool_dic) != 0:

		posit = []
		number_false = 0
		for key in tool_dic:
			posit.append(key)
		for item in posit:
			if not item in snps:
				number_false = number_false + 1

		return number_false

def analysis(total_list, ref_true, snps, homo, homo_snps):
	"""Writes in a txt file the location of the SNPs generated by HaploGenerator and if it is 
	detected or not by the filtering tool

	total_list = list containing all the snps generated by HaploGenerator: homeologous and true allelic SNPs
	ref_true = dictionary containing the true allelic SNPs generated by HaploGenerator
	snps =  list of true allelic SNPs detected by the filtering tool
        homo = dictionary containing the homeologous SNPs generated by HaploGenerator
	homo_snps = list of homeologous SNPs detected by the filtering tool"""

	output = open("REPORT.txt", "a")

	analysis = ""
	i = 0
	total_true = sorted(ref_dic)
	homo = sorted(homo)
	sorted(total_list)
	recount = []

	while i < len(sorted(total_list)):
		if str(sorted(total_list)[i]) in homo_snps:
                        if i == 0:
                                analysis = analysis + "False SNP " + str(sorted(total_list)[i]) + " is not Anchored" + " \n"
				recount.append(str(sorted(total_list)[i]))

                        if i != 0 and sorted(total_list)[i] != sorted(total_list)[-1]:

                                if str(sorted(total_list)[i-1]) in homo and str(sorted(total_list)[i+1]) in homo:
                                        analysis = analysis + "False SNP " + str(sorted(total_list)[i]) +  " is Anchored (" + str(sorted(total_list)[i-1]) + "-" + str(sorted(total_list)[i+1]) + ")"  + " \n"
					recount.append(str(sorted(total_list)[i]))

				if str(sorted(total_list)[i-1]) in total_true:

					if len(total_true) > (i+1):
	                                	if str(sorted(total_list)[i-2]) in homo and str(sorted(total_list)[i+1]) in homo:
        	                                	analysis = analysis + "False SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i-1]) + ")" + "Anchors: (" + str(sorted(total_list)[i-2]) + "-" + str(sorted(total_list)[i+1]) + " 1" + " \n"
                	                                recount.append(str(sorted(total_list)[i]))
						else:
							analysis = analysis + "False SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i-1]) + ")" + " 2 \n"
                                                        recount.append(str(sorted(total_list)[i]))

					if len(total_list) < (i+2):
                                                analysis = analysis + "False SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ") " + "That is not Anchored 3 \n"
                                                recount.append(str(sorted(total_list)[i]))

                                if str(sorted(total_list)[i+1]) in total_true:

					if len(total_list) > (i+2):
	                                	if str(sorted(total_list)[i-1]) in homo and str(sorted(total_list)[i+2]) in homo:
        	                                	analysis = analysis + "False SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ") " + "Anchors: (" + str(sorted(total_list)[i-1]) + "-" + str(sorted(total_list)[i+2]) + " 1" + " \n"
                	                                recount.append(str(sorted(total_list)[i]))
						else:
							analysis = analysis + "False SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ")" + " 2 \n"
                                        	        recount.append(str(sorted(total_list)[i]))

					if len(total_list) < (i+2):
						analysis = analysis + "False SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ") " + "That is not Anchored 3 \n"
                                                recount.append(str(sorted(total_list)[i]))


                        if sorted(total_list)[i] == sorted(total_list)[-1]:
                                analysis = analysis + "False SNP " + str(sorted(total_list)[i]) + " is not Anchored" + " \n"
				recount.append(str(sorted(total_list)[i]))


		if str(sorted(total_list)[i]) in total_true:

			if not str(sorted(total_list)[i]) in snps.keys():
				if i == 0:
                               		analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) + " is not Anchored"+ " Type: " + ref_dic[str(sorted(total_list)[i])][1] + " \n"
					recount.append(str(sorted(total_list)[i]))

                        	if i != 0 and sorted(total_list)[i] != sorted(total_list)[-1]:

	                               	if str(sorted(total_list)[i-1]) in homo and str(sorted(total_list)[i+1]) in homo:
                                       		analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) +  " is Anchored (" + str(sorted(total_list)[i-1]) + "-" + str(sorted(total_list)[i+1]) + ") " + ref_dic[str(sorted(total_list)[i])][1] + "\n"
						recount.append(str(sorted(total_list)[i]))

                               		if str(sorted(total_list)[i-1]) in total_true:
						if len(total_list) > (i+1):
							if str(sorted(total_list)[i-2]) in homo and str(sorted(total_list)[i+1]) in homo:
		                                       		analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i-1]) + ")" + " Anchors: (" + str(sorted(total_list)[i-2]) + "-" + str(sorted(total_list)[i+1]) + ") 1 " + ref_dic[str(sorted(total_list)[i])][1] + " \n"
								recount.append(str(sorted(total_list)[i]))
							else:
                                	                        analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i-1]) + ") 2 " + ref_dic[str(sorted(total_list)[i])][1] + " \n"
                                        	                recount.append(str(sorted(total_list)[i]))

						if len(total_list) < (i + 1):
							analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i-1]) + ") " + "That is not Anchored 3 " + ref_dic[str(sorted(total_list)[i])][1] + "\n"
	                                                recount.append(str(sorted(total_list)[i]))

                               		if str(sorted(total_list)[i+1]) in total_true:
						if len(total_list) > (i+2):
							if str(sorted(total_list)[i-1]) in homo and str(sorted(total_list)[i+2]) in homo:
        	                               			analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ")" + " Anchors: (" + str(sorted(total_list)[i-1]) + "-" + str(sorted(total_list)[i+2]) + ") 1 " + ref_dic[str(sorted(total_list)[i])][1] + " \n"
								recount.append(str(sorted(total_list)[i]))
							else:
								analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ")" + " 2 " + ref_dic[str(sorted(total_list)[i])][1] + " \n"
                                                        	recount.append(str(sorted(total_list)[i]))
						if len(total_list) < (i+2):
							analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ") " + "That is not Anchored 3 " + ref_dic[str(sorted(total_list)[i])][1] + " \n"
	                                                recount.append(str(sorted(total_list)[i]))


				if sorted(total_list)[i] == sorted(total_list)[-1]:
                               		analysis = analysis + "True SNP not detected " + str(sorted(total_list)[i]) + " is not Anchored " + ref_dic[str(sorted(total_list)[i])][1] + " \n"
					recount.append(str(sorted(total_list)[i]))


			if str(sorted(total_list)[i]) in snps.keys():
				if i == 0:
                                	analysis = analysis + "True SNP " + str(sorted(total_list)[i]) + " is not Anchored" + " Type: " + snps[str(sorted(total_list)[i])]  + " \n"
					recount.append(str(sorted(total_list)[i]))

                        	if i != 0 and sorted(total_list)[i] != sorted(total_list)[-1]:
                                	if str(sorted(total_list)[i-1]) in homo and str(sorted(total_list)[i+1]) in homo:
                                	        analysis = analysis + "True SNP " + str(sorted(total_list)[i]) +  " is Anchored (" + str(sorted(total_list)[i-1]) + "-" + str(sorted(total_list)[i+1]) + ")" + " Type: " + snps[str(sorted(total_list)[i])] + " \n"
						recount.append(str(sorted(total_list)[i]))

				   	if str(sorted(total_list)[i-1]) in total_true:

						if len(total_list) > (i+1):
	                                                if str(sorted(total_list)[i-2]) in homo and str(sorted(total_list)[i+1]) in homo:
        	                                                analysis = analysis + "True SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i-1]) + ")" + " Anchors: (" + str(sorted(total_list)[i-2]) + "-" + str(sorted(total_list)[i+1]) + ") 1" + " Type: " + snps[str(sorted(total_list)[i])] + " \n"
                	                                        recount.append(str(sorted(total_list)[i]))
                        	                        else:
                                	                        analysis = analysis + "True SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i-1]) + ") 2" + " Type: " + snps[str(sorted(total_list)[i])] + " \n"
                                        	                recount.append(str(sorted(total_list)[i]))
						if len(total_list) < (i+1):
							analysis = analysis + "True SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i-1]) + ") " + "That is not Anchored 3" + " Type: " + snps[str(sorted(total_list)[i])] +  " \n"
                                                        recount.append(str(sorted(total_list)[i]))

                                	if str(sorted(total_list)[i+1]) in total_true:
						if len(total_list) > (i + 2):
							if str(sorted(total_list)[i-1]) in homo and str(sorted(total_list)[i+2]) in homo:
        	                                                analysis = analysis + "True SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ")" + " Anchors: (" + str(sorted(total_list)[i-1]) + "-" + str(sorted(total_list)[i+2]) + ") 1" + " Type: " + snps[str(sorted(total_list)[i])] + " \n"
                	                                        recount.append(str(sorted(total_list)[i]))
                        	                        else:
                                	                        analysis = analysis + "True SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ") 2"+ " Type: " + snps[str(sorted(total_list)[i])]  + " \n"
	                                       	                recount.append(str(sorted(total_list)[i]))
						if len(total_list) < (i + 2):
							analysis = analysis + "True SNP " + str(sorted(total_list)[i]) +  " is next to a True SNP (" + str(sorted(total_list)[i+1]) + ") " + "That is not Anchored 3" + " Type: " + snps[str(sorted(total_list)[i])] + " \n"
                                                        recount.append(str(sorted(total_list)[i]))


                       		if sorted(total_list)[i] == sorted(total_list)[-1]:
                               		analysis = analysis + "True SNP " + str(sorted(total_list)[i]) + " is not Anchored" + " Type: " + snps[str(sorted(total_list)[i])] + " \n"
					recount.append(str(sorted(total_list)[i]))

		i = i + 1



	output.write(analysis + "++++++++++++++++" + "\n"  + "\n"  + "\n")
	output.close()


def writer(name, ref_true, tool_dic, snps, snps_called, number_false, homo_snps, total_list, true_per, false_per):
	"""Writes in a txt file the filtering tool performance in SNP detection

	name = name of the VCF file
	total_list = list containing all the snps generated by HaploGenerator: homeologous and true allelic SNPs
        ref_true = dictionary containing the true allelic SNPs generated by HaploGenerator
	tool_dic = dictionary containing the SNPs detected by the filtering tool
        snps =  list of true allelic SNPs detected by the filtering tool
	snps_called = true allelic SNPs detected percentage
	number_false = false positives detected by the filtering tool
        homo_snps = list of homeologous SNPs detected by the filtering tool
	true_per = percentage of true allelic SNPs in the genome
	false_per = percentage of homeologous SNPs in the genome """

	output = open("REPORT.txt", "a")

	output.write(">>>" + name + "\n")

	if len(tool_dic) == 0:
		output.write("FILTERING TOOL DIDN'T FIND ANYTHING" + "\n")
	if len(tool_dic) != 0:
        	output.write("Homeologous % in Filtering Tool: " + str(len(homo_snps)/len(tool_dic)*100)[0:5] + " %" + "\n")
		output.write("False positives: " + str(((number_false-(len(homo_snps)))/len(tool_dic)*100)[0:5] + " %" + "\n")
		output.write("False + homeologous % in Filtering Tool: " + str(number_false/len(tool_dic)*100)[0:5] + " %" + "\n")

	output.write(snps_called + "\n")
        output.write("Homeologous SNPs in Filtering Tool: " + str(len(homo_snps)) + "\n")
	output.write("True SNPs: " + str(len(snps)) + "\n")
	output.write("Total SNPs in SWEEP: " + str(len(tool_dic)) + "\n")
	output.write("Total SNPs in Reference: " + str(len(ref_true)) + "\n")
	output.write("True + Homeologous SNPs: " + str(len(total_list)) + "\n")
	output.write("True %: " + str(true_per)[0:5] + " % \n" + "False %: " + str(false_per)[0:5] + " % \n \n")

	output.close()

if __name__ == "__main__":


	name = str(argv[2])
        ref_true, homo, total_list, true_per, false_per = reference(argv[1], open(argv[3]), argv[4])
	tool_dic = filter_out(open(argv[2]))

	if len(tool_dic) != 0:
		snps, snps_called, homo_snps = true_positives(tool_dic, ref_dic, homo)
		number_false = false_positives(tool_dic, snps)
		writer(name, ref_true, tool_dic, snps, snps_called, number_false, homo_snps, total_list, true_per, false_per)
		analysis(total_dic, ref_dic, snps, homo, homo_snps)

	if len(tool_dic) == 0:
		snps = []
		snps_called_gatk = []
		snps_called = "No SNPs detected"
		homo_snps = []
		number = 0
		writer(name, ref_true, tool_dic, snps, snps_called, number_false, homo_snps, total_list, true_per, false_per)
