#!/bin/bash
module load juicer
module load java
module load bwa

#juicer.sh -g mm9 -s MboI -a $1 -y /opt/apps/juicer/1.6.2/restriction_sites/mm9_MboI.txt -D /opt/apps/juicer/1.6.2/ -q general -l general
/home/yonghao/softwares/juicer/scripts/juicer.sh -g mm9 -s MboI -a $1 -y /home/yonghao/softwares/juicer/restriction_sites/mm9_MboI.txt -D /home/yonghao/softwares/juicer/ -q general -l general
