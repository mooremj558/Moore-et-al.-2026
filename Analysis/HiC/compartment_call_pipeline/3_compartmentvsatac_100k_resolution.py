import os

file_input = os.environ["input"]
file_peak = os.environ["peaknumber"]
file_split = os.environ["split"]
output = os.environ["output"]

chr_comp = []
start_comp = []
stop_comp = []
eigen_comp = []
peaknumber = []
chr_split = []
split = []

with open(file_input, 'r') as file1:
    for line in file1:
        chr, start, stop, eigen = line.strip('\n').split('\t')
        chr_comp.append(chr)
        start_comp.append(int(start))
        stop_comp.append(int(stop))
        eigen_comp.append(float(eigen))
with open(file_peak, 'r') as file2:
    for line in file2:
         number = line.strip('\n')
         peaknumber.append(int(number))
with open(file_split, 'r') as file3:
    for line in file3:
        chr, s = line.strip('\n').split('\t')
        chr_split.append(chr)
        split.append(int(s))

chr_flag = "chr1"
break_point = 125000000
p_plus_bin = 0
p_plus_peak = 0
p_minus_bin = 0
p_minus_peak = 0
p_zero_bin = 0
p_zero_peak = 0
q_plus_bin = 0
q_plus_peak = 0
q_minus_bin = 0
q_minus_peak = 0
q_zero_bin = 0
q_zero_peak = 0
flag = 0
for chr, start, stop, eigen, peak in zip(chr_comp, start_comp, stop_comp, eigen_comp, peaknumber):
    # divide into two arms
    flag += 1
    if chr != chr_flag:
        total_bin = p_plus_bin + p_zero_bin + p_minus_bin + q_plus_bin + q_zero_bin + q_minus_bin
        p_total_bin = p_plus_bin + p_zero_bin + p_minus_bin
        q_total_bin = q_plus_bin + q_zero_bin + q_minus_bin
        total_peak = p_plus_peak + p_zero_peak + p_minus_peak + q_plus_peak + q_zero_peak + q_minus_peak
        if break_point != 0:
            print("%s\t%d\t%d\t%.2f\t%d\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (
            chr_flag, total_bin, total_peak, (p_plus_bin+q_plus_bin)/total_bin*100,p_plus_bin+q_plus_bin, (p_plus_peak+q_plus_peak)/(p_plus_bin+q_plus_bin), (p_minus_bin+q_minus_bin)/total_bin*100, p_minus_bin+q_minus_bin,(p_minus_peak+q_minus_peak)/(p_minus_bin+q_minus_bin), p_plus_bin/p_total_bin*100, p_plus_peak/p_plus_bin, p_minus_bin/p_total_bin*100, p_minus_peak/p_minus_bin, p_plus_peak/p_plus_bin-p_minus_peak/p_minus_bin , q_plus_bin/q_total_bin*100, q_plus_peak/q_plus_bin, q_minus_bin/q_total_bin*100, q_minus_peak/q_minus_bin, q_plus_peak/q_plus_bin-q_minus_peak/q_minus_bin), file=open(output, 'a+'))
        if break_point == 0:
            print("%s\t%d\t%d\t%.2f\t%d\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (
            chr_flag, total_bin, total_peak, (p_plus_bin+q_plus_bin)/total_bin*100,p_plus_bin+q_plus_bin, (p_plus_peak+q_plus_peak)/(p_plus_bin+q_plus_bin), (p_minus_bin+q_minus_bin)/total_bin*100, p_minus_bin+q_minus_bin,(p_minus_peak+q_minus_peak)/(p_minus_bin+q_minus_bin), 0, 0, 0, 0, 0 , q_plus_bin/q_total_bin*100, q_plus_peak/q_plus_bin, q_minus_bin/q_total_bin*100, q_minus_peak/q_minus_bin, q_plus_peak/q_plus_bin-q_minus_peak/q_minus_bin), file=open(output, 'a+'))         
        p_plus_bin = 0
        p_plus_peak = 0
        p_minus_bin = 0
        p_minus_peak = 0
        p_zero_bin = 0
        p_zero_peak = 0
        q_plus_bin = 0
        q_plus_peak = 0
        q_minus_bin = 0
        q_minus_peak = 0
        q_zero_bin = 0
        q_zero_peak = 0
        chr_flag = chr
        for chr_s, s in zip(chr_split, split):
            if chr == chr_s:
                break_point = s
    if stop <= break_point:
        if eigen < 0:
            p_minus_bin += 1
            p_minus_peak += peak
        if eigen == 0:
            p_zero_bin += 1
            p_zero_peak += peak
        if eigen > 0:
            p_plus_bin += 1
            p_plus_peak += peak
    if start >= break_point:
        if eigen < 0:
            q_minus_bin += 1
            q_minus_peak += peak
        if eigen == 0:
            q_zero_bin += 1
            q_zero_peak += peak
        if eigen > 0:
            q_plus_bin += 1
            q_plus_peak += peak
    if flag == 30321:
        total_bin = p_plus_bin + p_zero_bin + p_minus_bin + q_plus_bin + q_zero_bin + q_minus_bin
        p_total_bin = p_plus_bin + p_zero_bin + p_minus_bin
        q_total_bin = q_plus_bin + q_zero_bin + q_minus_bin
        total_peak = p_plus_peak + p_zero_peak + p_minus_peak + q_plus_peak + q_zero_peak + q_minus_peak
        print("%s\t%d\t%d\t%.2f\t%d\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" % (
        chr_flag, total_bin, total_peak, (p_plus_bin+q_plus_bin)/total_bin*100,p_plus_bin+q_plus_bin, (p_plus_peak+q_plus_peak)/(p_plus_bin+q_plus_bin), (p_minus_bin+q_minus_bin)/total_bin*100, p_minus_bin+q_minus_bin,(p_minus_peak+q_minus_peak)/(p_minus_bin+q_minus_bin), 0, 0, 0, 0, 0 , q_plus_bin/q_total_bin*100, q_plus_peak/q_plus_bin, q_minus_bin/q_total_bin*100, q_minus_peak/q_minus_bin, q_plus_peak/q_plus_bin-q_minus_peak/q_minus_bin), file=open(output, 'a+'))  
