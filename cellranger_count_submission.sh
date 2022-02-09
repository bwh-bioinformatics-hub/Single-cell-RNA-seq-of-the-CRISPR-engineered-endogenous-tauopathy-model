#!/bin/bash
for i in $(cat /data/bioinformatics/projects/hassan2022/data/FC_07189/220118_A01061_0232_BH72NCDMXY/BPF_library_Feany_lab_mod.csv | awk -F "," '(NR>1){print $2}'); do
    bsub -q big -e ${i}_count.log "sh cellranger_count.sh $i"
done
