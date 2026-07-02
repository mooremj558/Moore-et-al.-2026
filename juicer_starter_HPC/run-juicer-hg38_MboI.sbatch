#!/bin/bash

#module load juicer
#module load java
#module load bwa
#module load cuda/8.0.61

#juicer.sh -g hg38 -s MboI -a $1 -y /opt/apps/juicer/1.6.2/restriction_sites/hg38_MboI.txt -D /opt/apps/juicer/1.6.2/ -q general -l general
/wanglab/yliang/softwares/juicer/scripts/juicer_HPC.sh -g hg38 -s MboI -a $1 -y /wanglab/yliang/softwares/juicer/restriction_sites/hg38_MboI.txt -D /wanglab/yliang/softwares/juicer -q general -l general

