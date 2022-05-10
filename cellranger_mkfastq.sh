
# cellranger mkfastq
conda activate scrnaseq
module load cellranger/6.0
bsub -Is -R 'rusage[mem=10000]' -n 4 /bin/bash

cd /data/bioinformatics/projects/hassan2022/output/FC_07379
nohup cellranger mkfastq --id=01_cellranger_mkfastq \
                     --run=/data/bioinformatics/projects/hassan2022/data/FC_07379/220419_A01061_0302_AHWKHHDRXY \
                     --csv=/data/bioinformatics/projects/hassan2022/data/FC_07379/220419_A01061_0302_AHWKHHDRXY/SampleSheet_FC007379_mod.csv

                                        
