scRNA-seq data analysis
Anthony Cicalo
09222023

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
Run Rscript hassan_merged_seurat.r -l expectedCells/ -s scrublet/ -k outdir/ -j hassan2022 -r refdir path_to_ref_directory
This will generate Seurat results.

Step6:
Proceed with the trajectory and the fly phone DB analysis, using their respective code, on the integrated Seurat object. 

Step7: 
pySCENIC

Step8: 
OmicsIntegrator
