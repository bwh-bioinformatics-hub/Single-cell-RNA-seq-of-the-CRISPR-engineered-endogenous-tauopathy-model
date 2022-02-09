# 02012022
# download data to ERIS
scp -r SCCcollabIGenomix@bpfngs.med.harvard.edu:./FC_07189 .
Password: 0ny2O9BMrjgF

# the above code will download a dir called "FC_07189"
cd FC_07189/
ls
220118_A01061_0232_BH72NCDMXY.tar.gz
220118_A01061_0232_BH72NCDMXY.tar.gz.md5

# unzip gz file
tar -xvf 220118_A01061_0232_BH72NCDMXY.tar.gz

# the above line of code will generate a dir "n"
/data/bioinformatics/projects/hassan2022/FC_07189/n/files/Genetics/BPF-NGS/novaseq/220118_A01061_0232_BH72NCDMXY
cd /data/bioinformatics/projects/hassan2022/
mkdir data
cd /data/bioinformatics/projects/hassan2022/
mv 220118_A01061_0232_BH72NCDMXY /data/bioinformatics/projects/hassan2022/data/FC_07189/

module load bcl2fastq/2.22.0
module load cellranger/6.1.1

# error while using bcl2fastq cellranger mkfastq, email hpcsupport
# install bcl2fastq on conda
conda create --name scrnaseq
conda activate scrnaseq
conda install -c dranew bcl2fastq

# cellranger mkfastq
# modify sample sheet and rerun
cellranger mkfastq --id=01_cellranger_mkfastq \
                     --run=/data/bioinformatics/projects/hassan2022/data/FC_07189/220118_A01061_0232_BH72NCDMXY \
                     --csv=/data/bioinformatics/projects/hassan2022/data/FC_07189/220118_A01061_0232_BH72NCDMXY/BPF_library_Feany_lab_mod.csv

# Location of Dm reference genome, ERISONE
/data/bioinformatics/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/Drosophila_melanogaster.genome

# Location of fastq files.
/data/bioinformatics/projects/hassan2022/output/01_cellranger_mkfastq/outs/fastq_path/

# Ready for cellranger count, ERISone
cd /data/bioinformatics/projects/hassan2022/output
mkdir 02_cellranger_count
cd 02_cellranger_count

# Run bash script for multiple samples, this script will call another bash file, csv file is the same as used for cellranger mkfastq index file
bash cellranger_count_submission.sh

############################ cellranger_count_submission.sh
#!/bin/bash
for i in $(cat /data/bioinformatics/projects/hassan2022/data/FC_07189/220118_A01061_0232_BH72NCDMXY/BPF_library_Feany_lab_mod.csv | awk -F "," '(NR>1){print $2}'); do
    bsub -q big -e ${i}_count.log "sh cellranger_count.sh $i"
done
###########################################################

############################ cellranger_count.sh
#!/bin/bash
annotation=$1
conda activate scrnaseq
module load cellranger/6.1.1
cellranger count  --id=$annotation \
                  --transcriptome=/data/bioinformatics/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/Drosophila_melanogaster.genome \
                  --fastqs=/data/bioinformatics/projects/hassan2022/output/01_cellranger_mkfastq/outs/fastq_path/ \
                  --sample=$annotation \
                  --include-introns \
                  --localcores=8 \
                  --localmem=64
###########################################################


# to find out mapping parameters

