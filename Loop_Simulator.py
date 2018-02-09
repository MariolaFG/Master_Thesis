#!/usr/bin/env python
"""This pipeline runs the simulation of the genome and the SNP detection 
by different tools as many times as given. It is necesary to have a genome to
sample in path

Usage: python Loop_Simulator.py number_of_loops ploidy_level number_of_organisms"""

from __future__ import division
from Genome_Simulator import simulator_ploidy
import subprocess
from sys import argv
import random

def gen_break(start, end, genom):
    """ Picks part of a genome to create a smaller data set as big as desired
    
    start: position of the genome to be the first position in the new genome
    end: position of the genome to be the last position in the new genome
    genom: new genome's file"""
    
    
    little = open(genom, "w")
    breaker = open("Chromosome_file", "r")

    tot = ""
    i = 0

    for line in breaker:
        i = i + 1
        if start < i < end:
            line = line.strip()
            tot = tot + line

    little.write(">Chromosome_name" + "\n")

    e = 0
    e2 = 59
    while e < len(tot):
        if e < e2:
            little.write(tot[e])
            e = e + 1
            if e == e2:
                little.write(tot[e] + "\n")
                e = e + 1
                e2 = e2 + 60

    breaker.close()
    little.close()

def simulation_loop(number, ploidy, org):
    """ Runs the loop as many times as desired
    
    number: number of experiments desired
    ploidy: ploidy level
    org: number of organisms"""
    
    genom = 'Trial.fasta'
    o = 0
    sum_end = 350
    while o < int(number):
        start = random.randint(1, 10000)
        end = start + sum_end
        gen_break(start, end, genom)        
        #This produces the new genome to simulate SNPs in

        cmd = "samtools faidx Trial.fasta"
        e = subprocess.check_call(cmd, shell=True)

        simulator_ploidy(genom, ploidy, org)
        # Uses the simulation pipeline
        
        i = 5
        cove = 0
        while i < 100:
            
            cove = cove + 1
            cmd = "python Genome_reads.py Trial.fasta %s SWEEP %s 100 %s" % (ploidy, str(org), str(i))
            e = subprocess.check_call(cmd, shell=True)
            # Genome_reads is the same pipeline as Genome_Simulator, without the
            # simulation part to be able to do coverage screening.
            
            org_num = 1
            inpp = ""
            inpp2 = ""
            while org_num < int(org) + 1:
                inpp = inpp + "-b Alignment/sorted_merged" + str(org_num) + ".bam "
                if inpp2 == "":
                    inpp2 = inpp2 + "Simulated_Genome/genTotal" + str(org_num) + "_varianthaplos.txt"
                elif inpp2 != "":
                    inpp2 = inpp2 + ",Simulated_Genome/genTotal" + str(org_num) + "_varianthaplos.txt"
                org_num = org_num + 1
            # It gathers all BAM files and the variant reference files.    
                  
                        
            cmd = "perl ../../SWEEP.pl %s -g Trial.fasta -o SWEEP1_out.vcf" % (inpp)
            e = subprocess.check_call(cmd, shell=True)
                      
            cmd = "python ../../Compare.py %s SWEEP1_out.vcf Trial.fasta %s" % (inpp2, org)
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python Reporter.py %s %s" % (str(o+1), str(cove) + "_" + str(1))
            e = subprocess.check_call(cmd, shell=True)
            
            # Runs original SWEEP, compares with reference and stores data 
            # in a txt file
            
            cmd = "rm REPORT.txt"
            e = subprocess.check_call(cmd, shell=True)
            
            open_report = open("REPORT.txt", "w")
            open_report.write("\n")
            open_report.close()
            
            
            cmd = "perl ../../Modified_SWEEP.pl -b Alignment/sorted_merged1.bam -b Alignment/sorted_merged2.bam -g Trial.fasta -o SWEEP2_out.vcf --no_cleanup"
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python ../../Compare.py %s SWEEP2_out.vcf Trial.fasta %s" % (inpp2, org)
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python Reporter.py %s %s" % (str(o+1), str(cove) + "_" + str(2))
            e = subprocess.check_call(cmd, shell=True)
            
            # Runs Modified SWEEP, compares with reference and stores data 
            # in a txt file
            
            cmd = "rm REPORT.txt"
            e = subprocess.check_call(cmd, shell=True)
            
            open_report = open("REPORT.txt", "w")
            open_report.write("\n")
            open_report.close()
            
            
            
            cmd = "../../freebayes/bin/freebayes -f Trial.fasta -p %s --min-coverage %s --min-base-quality 9 --min-mapping-quality 13 -F 0.1  -b Alignment/sorted_merged1.bam > free_out.vcf" % (ploidy, ploidy)
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "../../freebayes/bin/freebayes -f Trial.fasta -p %s --min-coverage %s --min-base-quality 9 --min-mapping-quality 13 -F 0.1  -b Alignment/sorted_merged2.bam > free_out2.vcf" % (ploidy, ploidy)
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python Free_merged.py free_out.vcf free_out2.vcf free_total.vcf"
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python Free_filter.py free_total.vcf free_filter.vcf"
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python ../../Compare.py %s free_total.vcf Trial.fasta %s" % (inpp2, org)
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python Reporter.py %s %s" % (str(o+1), str(cove) + "_" + str(3))
            e = subprocess.check_call(cmd, shell=True)
            
            # Runs Filtered FreeBayes, compares with reference and stores data 
            # from total FreeBayes in a txt file
            
            cmd = "rm REPORT.txt"
            e = subprocess.check_call(cmd, shell=True)
            
            open_report = open("REPORT.txt", "w")
            open_report.write("\n")
            open_report.close()
            
            cmd = "python ../../Compare.py %s free_filter.vcf Trial.fasta %s" % (inpp2, org)
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python Reporter.py %s %s" % (str(o+1), str(cove) + "_" + str(6))
            e = subprocess.check_call(cmd, shell=True)
            
            # Stores data from filtered FreeBayes in a txt file
            
            cmd = "rm REPORT.txt"
            e = subprocess.check_call(cmd, shell=True)
            
            open_report = open("REPORT.txt", "w")
            open_report.write("\n")
            open_report.close()

            cmd = "python ../../Compare.py %s prefilter.vcf Trial.fasta %s" % (inpp2, org)
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "python Reporter.py %s %s" % (str(o+1), str(cove) + "_" + str(4))
            e = subprocess.check_call(cmd, shell=True)
            
            # Stores data from Bcftools in a txt file
            
            cmd = "rm REPORT.txt"
            e = subprocess.check_call(cmd, shell=True)
            
            open_report = open("REPORT.txt", "w")
            open_report.write("\n")
            open_report.close()
            
            cmd = "rm *.bcf *.vcf"
            e = subprocess.check_call(cmd, shell=True)
                
            cmd = "rm Alignment/*"
            subprocess.check_call(cmd, shell=True)
            
            cmd = "rm Reads/*"
            e = subprocess.check_call(cmd, shell=True)
            
            cmd = "rm Simulated_Genome/*"
            e = subprocess.check_call(cmd, shell=True)
            
            
            i = i + 100
        o = o + 1



if __name__ == "__main__":

	simulation_loop(argv[1], argv[2], argv[3])
