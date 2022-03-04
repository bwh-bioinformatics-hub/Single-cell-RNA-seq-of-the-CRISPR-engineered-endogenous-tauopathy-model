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
