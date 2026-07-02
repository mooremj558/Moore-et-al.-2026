#!/usr/bin/env python3
"""
Normalize eigenvalue to [-1, 1]

Usage: python3 5_normalization.py --input <input eigenvalue bedgraph file>
"""

import argparse
import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns

# Get the path to input files
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "input eigenvalue bedgraph file", type = str)

args = parser.parse_args()
input_file = args.input

###############################################################
# Part 1. Data Read in
###############################################################
output_file = os.path.splitext(input_file)[0] + ".normalized.bedgraph"

data = {}
eigenvalue = []

with open(input_file, 'r') as file:
	index = -1
	for line in file:
		index += 1
		entry = line.strip('\n').split('\t')
		data[index] = [entry[0], entry[1], entry[2], float(entry[3])] #data[index] -> [chr, start, stop, eigenvalue]
		if entry[3] != 'nan' and entry[3] != 'NaN':
			eigenvalue.append(float(entry[3]))

max_eigenvalue = max(eigenvalue)

for index in data:
	data[index].append(data[index][3]/max_eigenvalue) #data[index] -> [chr, start, stop, eigenvalue, normalized eigenvalue]

print(input_file, max_eigenvalue)

for index in data:
	print('\t'.join(data[index][0:3] + [str(data[index][4])]), file = open(output_file, 'a+'))
