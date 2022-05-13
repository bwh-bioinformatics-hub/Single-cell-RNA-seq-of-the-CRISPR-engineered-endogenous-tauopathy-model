scRNA-seq data analysis
Tingting Zhao
05092022

Step1:
bash cellranger_mkfastq.sh
This will generate fastq files.

Step2:
bash cellranger_count_submission.sh
This will call another bash file “cellranger_count.sh” and generate count tables

Step3:
bash run_scrublet_multi.sh
This will call “scrublet_multi.py” and generate scrublet results. 
If seeing error refer to “scrublet_multi_conditional.py”.

Step4:
bash scRNA_seq.sh
This will call “seurat_individual.R” and generate QC plots for each sample.

Step5:
hassan_merged.rmd in Rstudio
This will generate Seurat results.

