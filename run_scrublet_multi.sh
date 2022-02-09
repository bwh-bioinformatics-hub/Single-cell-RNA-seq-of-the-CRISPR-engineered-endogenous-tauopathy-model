PREFIX="/Volumes/BIOINFORMATICS/projects/hassan2022/output/02_cellranger_count/"
SUFFIX="/outs/filtered_feature_bc_matrix"
OUTS="/Volumes/BIOINFORMATICS/projects/hassan2022/output/03_scrublet/"
SAMPLES=("CRN00224913" "CRN00224914" "CRN00224915" "CRN00224916" "CRN00224917" "CRN00224918" "CRN00224919" "CRN00224920")

#~/anaconda3/bin/python scrublet_multi.py $PREFIX $SUFFIX $OUTS "${SAMPLES[@]}" # if run scripts on server
/Users/tingtingzhao/anaconda3/bin/python3 scrublet_multi.py $PREFIX $SUFFIX $OUTS "${SAMPLES[@]}"
