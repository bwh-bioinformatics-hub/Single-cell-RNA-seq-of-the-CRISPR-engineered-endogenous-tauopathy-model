# title: script for single-cell data analysis
# author: Tingting Zhao
# email: tzhao7@bwh.harvard.edu
# date: "02/05/2022"
# usage: bash main.sh /data/bioinformatics/projects/hassan2022/data/FC_07189/220118_A01061_0232_BH72NCDMXY/BPF_library_Feany_lab_mod.csv /data/bioinformatics/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/Drosophila_melanogaster.genome /data/bioinformatics/projects/hassan2022/ /data/bioinformatics/projects/hassan2022/data/ /data/bioinformatics/projects/hassan2022/data/FC_07189/220118_A01061_0232_BH72NCDMXY /data/bioinformatics/projects/hassan2022/scr/ /data/bioinformatics/projects/hassan2022/output/02_cellranger_count_expectedCells/ /data/bioinformatics/projects/hassan2022/output/04_seurat_individual/ /data/bioinformatics/projects/hassan2022/output/03_scrublet/ CRN00224913,CRN00224914,CRN00224915,CRN00224916,CRN00224917,CRN00224918,CRN00224919,CRN00224920 hassan2022 https://docs.google.com/spreadsheets/d/1gCzAeVe9Ekpyt8XdNOyBrOcw7Letr-WTpBh4-37ySpA/edit#gid=446579886

# folder structure: project folder (data, output, result, scr)

# setting
args <- commandArgs(trailingOnly = TRUE)

index_file = args[1]
reference_folder = args[2]
project_folder = args[3]
data_folder = args[4]
#baseSpace_ID =
bcl_files = args[5]
expect_cells = 10000
fastq_folder = "01_cellranger_mkfastq"
countTable_folder = "02_cellranger_count_expectedCells"
pwd = args[6]
indir = args[7]
outdir = args[8]
scrubletdir = args[9]
samples = unlist(strsplit(args[10], ','))
projectName = args[11]
marker_link = args[12]
marker_sheet = "MarkerGenesFiltered"

#Options for cell clustering algorithm, 1=louvain, 2=GLMPCA, 3= leiden, 0=louvain and GLMPCA and leiden
flag=1

conda activate scrnaseq
module load cellranger/6.1.1 #install under conda

# step1: download data from BaseSpace
#bs download run -i $baseSpace_ID -o $data_folder

# step2: bcl to fastq
cd $project_folder/output
bsub -q big cellranger mkfastq --id =$fastq_folder \
                   --run =$bcl_files \
                   --csv =$index_file

# step3: making count table
cd $project_folder/output/$countTable_folder
for i in $(cat $index_file | awk -F "," '(NR>1){print $2}'); do
bsub -q big cellranger count --id = $i \
--transcriptome = $reference_folder \
--fastqs =$project_folder/output/$fastq_folder/outs/fastq_path/ \
--sample =$i \
--expect-cells =$expect_cells \
--include-introns
done

# step4: running Seurat
bsub -q big -e seurat_individual.log Rscript seurat_individual.R $pwd $indir $outdir $scrubletdir $samples $projectName $markers

