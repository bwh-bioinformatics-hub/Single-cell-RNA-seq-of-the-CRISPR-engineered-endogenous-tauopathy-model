#!/bin/bash
#conda activate scrnaseq
module load cellranger/6.0
cellranger count  --id=P251L1merged \
                  --transcriptome=/data/bioinformatics/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/Drosophila_melanogaster.genome \
                  --fastqs=/data/bioinformatics/projects/hassan2022/output/FC_07189/01_cellranger_mkfastq/outs/fastq_path/H72NCDMXY/,/data/bioinformatics/projects/rachel2022/output/FC_07296/01_cellranger_mkfastq/outs/fastq_path/H7KGTDMXY/ \
                  --sample=CRN00224913,P251L-1 \
                  --expect-cells=10000

cellranger count  --id=CTRL2merged \
                  --transcriptome=/data/bioinformatics/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/Drosophila_melanogaster.genome \
                  --fastqs=/data/bioinformatics/projects/hassan2022/output/FC_07189/01_cellranger_mkfastq/outs/fastq_path/H72NCDMXY/,/data/bioinformatics/projects/rachel2022/output/FC_07296/01_cellranger_mkfastq/outs/fastq_path/H7KGTDMXY/ \
                  --sample=CRN00224914,CTRL-2 \
                  --expect-cells=10000
                  
cellranger count  --id=TauKO1merged \
                  --transcriptome=/data/bioinformatics/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/Drosophila_melanogaster.genome \
                  --fastqs=/data/bioinformatics/projects/hassan2022/output/FC_07189/01_cellranger_mkfastq/outs/fastq_path/H72NCDMXY/,/data/bioinformatics/projects/rachel2022/output/FC_07296/01_cellranger_mkfastq/outs/fastq_path/H7KGTDMXY/ \
                  --sample=CRN00224915,Tau-KO-1 \
                  --expect-cells=10000

cellranger count  --id=P251L2merged \
                  --transcriptome=/data/bioinformatics/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/Drosophila_melanogaster.genome \
                  --fastqs=/data/bioinformatics/projects/hassan2022/output/FC_07189/01_cellranger_mkfastq/outs/fastq_path/H72NCDMXY/,/data/bioinformatics/projects/rachel2022/output/FC_07296/01_cellranger_mkfastq/outs/fastq_path/H7KGTDMXY/ \
                  --sample=CRN00224916,P251L-2 \
                  --expect-cells=10000
