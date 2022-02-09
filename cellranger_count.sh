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


