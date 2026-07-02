#!/usr/bin/env python3
"""
When using bing_2015 pipeline to call A/B compartment, +/- of eigenvalue is assigned using gene density. 
However, gene density is not a good predictor for A/B compartment since it will not change in different cell type or under different treatment.
So, if other datasets like ATAC-seq is available, it's a good idea to integrate with another dataset within the same cell type or under the same treatment.
This script is to flip the eigenvalue called based on gene density.

Usage: python3 compartment_flipper.py --input <input_compartment_file.bedgraph> --chrom <chrom#(23notX), seperated by ','> --arm <arm(p/q), separated by ','> --split <path to chromosome split file, 'chr# splie_site'>
"""

import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input the eigenvalue file in bedgraph format")
parser.add_argument("--chrom", help = "chromosome that needs to be flipped", type = str)
parser.add_argument("--arm", help = "specify which arm needs to be flipped", type = str)
parser.add_argument("--split", help = "path to chromosome split file", default = "/bar/yliang/genomes/private/hg38_chrom_split.txt", type = str)

args = parser.parse_args()

input_file = args.input
chrom_nums_input = args.chrom
chrom_arms_input = args.arm
split_file = args.split

output_file = os.path.splitext(os.path.basename(input_file))[0] + ".flipped.bedgraph"

###############################################################
# Part 1. Data Read in
###############################################################
split = {}

with open(split_file, 'r') as file:
	for line in file:
		entry = line.strip('\n').split('\t')
		split[entry[0]] = int(entry[1])

chrom_nums = chrom_nums_input.split(',')

chrom_arms = chrom_arms_input.split(',')

# Check if an argument was passed to the python script
if (len(chrom_nums) != len(chrom_arms)): 
    sys.exit("please enter chromosome number and its corresponding arm") 

eigen = {}
flag = 0

with open(input_file, 'r') as file:
	for line in file:
		entry = line.strip('\n').split('\t')
		eigen[flag] = [entry] # eigen_unflipped[flag] = [[chr, start, stop, eigen]]
		eigen[flag][0][3] = float(eigen[flag][0][3])
		flag += 1

###############################################################
# Part 2. Data Processing
###############################################################
for chr_num_num, chr_arm in zip(chrom_nums, chrom_arms):
	chr_num = "chr" + str(chr_num_num)
	chr_split = split[chr_num]
	for flag in eigen:
		if eigen[flag][-1][0] != chr_num:
			eigen[flag].append(eigen[flag][-1])
			continue
		if eigen[flag][-1][0] == chr_num:
			if chr_arm == "p" and int(eigen[flag][-1][2]) <= chr_split:
				eigen[flag].append([eigen[flag][-1][0], eigen[flag][-1][1], eigen[flag][-1][2], float(eigen[flag][-1][3])*-1])
				continue
			if chr_arm == "q" and int(eigen[flag][-1][1]) >= chr_split:
				eigen[flag].append([eigen[flag][-1][0], eigen[flag][-1][1], eigen[flag][-1][2], float(eigen[flag][-1][3])*-1])
				continue
			else:
				eigen[flag].append(eigen[flag][-1])
				continue

###############################################################
# Part 3. Data Output
###############################################################
#for flag in eigen:
#	print("%s\t%s\t%s\t%.9f" %(eigen[flag][-1][0], eigen[flag][-1][1], eigen[flag][-1][2], eigen[flag][-1][3]), file = open(output_file, 'a+'))

with open(output_file, 'w+') as file:
	for flag in eigen:
		file.write("%s\t%s\t%s\t%.9f\n" %(eigen[flag][-1][0], eigen[flag][-1][1], eigen[flag][-1][2], eigen[flag][-1][3]))
		#file.write('\t'.join([str(x) for x in [eigen[flag][-1][0], eigen[flag][-1][1], eigen[flag][-1][2], eigen[flag][-1][3]]]) + '\n')
