#!/bin/bash
#module load juicer
#module load java
#module load bwa

#juicer.sh -g mm10 -s MboI -a $1 -y /opt/apps/juicer/1.6.2/restriction_sites/mm10_MboI.txt -D /opt/apps/juicer/1.6.2/ -q general -l general
/wanglab/yliang/softwares/juicer/scripts/juicer_HPC.sh -g mm10 -s Arima -a $1 -y /wanglab/yliang/softwares/juicer/restriction_sites/mm10_Arima.txt -D /wanglab/yliang/softwares/juicer/ -q general -l general
