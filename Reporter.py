#!/usr/bin/env python


"""This pipeline writes in a csv file the data obtained by different filtering tools
Usage: python Reporter.py number_file tool_ID"""

from __future__ import division
from sys import argv

inpp = open("REPORT.txt", "r")
input = inpp.read().split("\n")

sweep0 = {}
wrong0 = []

start0 = 0
end = 0

i = 0
n = 1

number_o = argv[1]
id = argv[2]

for line in input:
	# Divides with different names used for filtering tools
	if "SWEEP1" in line:
		start0 = i
		supername = line

	if "SWEEP2" in line:
                start0 = i
                supername = line

	if "Bcftools" in line:
		start0 = i
		supername = line
	if "Free" in line:
		start0 = i
		supername = line
	if "free" in line:
                start0 = i
                supername = line

	if "+++" in line:

		end = i

	if start0 != 0 and start0 < end:
#		print start0, end

		if "SWEEP DIDN'T FIND ANYTHING" in input[start0:end]:
			wrong0.append(input[start0:end])
		else:
			sweep0[n] = input[start0:end]
		n = n + 1
		start0 = end

	i = i + 1

name_o = "Results" + number_o + "_C" + id + ".csv"
output = open(name_o, "w")
# Stores in this file all the statistical data

for key in sweep0:
	out1 = str(sweep0[key])
	out1 = out1.replace("[", "").replace("]", "").replace("'", "")
#	print out1
#	out = out1.split("-p 2")[4]
	name = out1[0]
	output.write(name + str(key) + "," + out1 + "\n")

for element in wrong0:
	out = str(element)
        out = out.replace("[", "").replace("]", "").replace("'", "")
        output.write(name + "Wrong," + str(key) + "," + out + "\n")

output.write("\n \n \n")


# This following code, analyzes each SNP type

s0 = {}
s = []
for key in sweep0:
	u = 0
	for element in sweep0[key]:
		if "Total SNPs in Reference" in element:
			s.append(sweep0[key][u+2:])
#			if int(sweep0[key][13][25:])/int(sweep0[key][14][25:]) < 0.34:
#				s.append(sweep0[key][u+2:])
		u = u + 1


s0["True Anchored not detected"] = []
s0["True Anchored detected"] = []
s0["True not Anchored not detected"] = []
s0["True not Anchored detected"] = []
s0["True next to a True SNP not detected"] = []
s0["True next to a True SNP detected"] = []


count = 0

while count < len(s):
	for element in s[count]:
		if "True" in element and "is Anchored" in element:
			if "not detected" in element:
				element = element.split(" ")
				snp = element[4]
				anchors = element[7].replace("(", "").replace(")", "").split("-")
				s0["True Anchored not detected"].append([snp + "," + anchors[0] + "," + anchors[1] + "," + element[8]])

			else:
				element = element.split(" ")
		                snp = element[2]
		                anchors = element[5].replace("(", "").replace(")", "").split("-")
		                s0["True Anchored detected"].append([snp + "," + anchors[0] + "," + anchors[1] + "," + element[7]])

		if "True" in element and "is not Anchored" in element:
	                if "not detected" in element:
	                        element = element.split(" ")
	                        snp = element[4]
	                        s0["True not Anchored not detected"].append([snp + "," + element[8]])

	                else:
	                        element = element.split(" ")
	                        snp = element[2]
	                        s0["True not Anchored detected"].append([snp + "," + element[7]])

		if not "False" in element and "is next to a True SNP" in element:
			if "not detected" in element:
	                        element = element.split(" ")
	                        snp = element[4]

				if len(element) > 14 and "1" in element[14]:
		                        ts = element[11].replace("(", "").replace(")", "")
					anchors = element[13].replace("(", "").replace(")", "").replace("-", ",")
		                        s0["True next to a True SNP not detected"].append([snp + "," + ts + "," + anchors + "," + element[15]])
				if "2" in element[12]:
					ts = element[11].replace("(", "").replace(")", "")
                                        s0["True next to a True SNP not detected"].append([snp + "," + ts + "," + "NOT" + "," + element[13]])

	                else:
	                        element = element.split(" ")
	                        snp = element[2]
				if len(element) > 13 and "1" in element[12]:
                                        ts = element[9].replace("(", "").replace(")", "")
                                        anchors = element[11].replace("(", "").replace(")", "").replace("-", ",")
                                        s0["True next to a True SNP detected"].append([snp + "," + ts + "," + anchors + "," + element[14]])
                                if "2" in element[10]:
                                        ts = element[9].replace("(", "").replace(")", "")
                                        s0["True next to a True SNP detected"].append([snp + "," + ts + "," + "NOT" + "," + element[12]])


		if "False" in element and "is Anchored" in element:
                        if "not detected" in element:
                                element = element.split(" ")
                                snp = element[4]
                                anchors = element[7].replace("(", "").replace(")", "").split("-")
                                s0["True Anchored not detected"].append(["False," + snp + "," + anchors[0] + "," + anchors[1]])

                        else:
                                element = element.split(" ")
                                snp = element[2]
                                anchors = element[5].replace("(", "").replace(")", "").split("-")
                                s0["True Anchored detected"].append(["False," + snp + "," + anchors[0] + "," + anchors[1]])

                if "False" in element and "is not Anchored" in element:
                        if "not detected" in element:
                                element = element.split(" ")
                                snp = element[4]
                                s0["True not Anchored not detected"].append(["False," + snp])

                        else:
                                element = element.split(" ")
                                snp = element[2]
                                s0["True not Anchored detected"].append(["False," + snp])

                if "False" in element and "is next to a True SNP" in element:
                        if "not detected" in element:
                                element = element.split(" ")
                                snp = element[4]

                                if len(element) > 14 and "1" in element[14]:
                                        ts = element[11].replace("(", "").replace(")", "")
                                        anchors = element[13].replace("(", "").replace(")", "").replace("-", ",")
                                        s0["True next to a True SNP not detected"].append(["False," + snp + "," + ts + "," + anchors])
                                if "2" in element[12]:
                                        ts = element[11].replace("(", "").replace(")", "")
                                        s0["True next to a True SNP not detected"].append(["False," + snp + "," + ts + "," + "NOT"])

                        else:
                                element = element.split(" ")
                                snp = element[2]
                                if len(element) > 13 and "1" in element[12]:
                                        ts = element[9].replace("(", "").replace(")", "")
                                        anchors = element[11].replace("(", "").replace(")", "").replace("-", ",")
                                        s0["True next to a True SNP detected"].append(["False," + snp + "," + ts + "," + anchors])
                                if "2" in element[10]:
                                        ts = element[9].replace("(", "").replace(")", "")
					s0["True next to a True SNP detected"].append(["False," + snp + "," + ts + "," + "NOT"])

	count = count + 1



name_o = "Analysis" + number_o + "_C" + organisms + ".csv"
analizer = open(name_o, "w")
analizer.write(supername + " \n")


for key in s0:
	for element in s0[key]:
		line = str(element).replace("[", "").replace("]", "").replace("'", "")
		analizer.write(key + ","  + line + "\n")

analizer.write("\n")
# SNP type and detection are stored in this second file


analizer.close()
output.close()
