#!/bin/bash
for i in $(cat /data/bioinformatics/projects/hassan2022/scr/SampleSheet_FC007379_mod.csv | awk -F "," '(NR>1){print $2}'); do
    bsub -q big -e ${i}_count.log -u tzhao7@bwh.harvard.edu "sh cellranger_count.sh $i"
done
