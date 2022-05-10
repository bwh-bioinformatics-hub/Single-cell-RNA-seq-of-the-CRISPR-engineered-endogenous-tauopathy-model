#!/bin/bash
annotation=$1
conda activate scrnaseq
module load cellranger/6.0
cd /data/bioinformatics/projects/hassan2022/output/FC_07379/02_cellranger_count_expectedCells
cellranger count  --id=$annotation \
                  --transcriptome=/data/bioinformatics/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/Drosophila_melanogaster.genome \
                  --fastqs=/data/bioinformatics/projects/hassan2022/output/FC_07379/01_cellranger_mkfastq/outs/fastq_path/ \
                  --sample=$annotation \
                  --expect-cells=10000 

