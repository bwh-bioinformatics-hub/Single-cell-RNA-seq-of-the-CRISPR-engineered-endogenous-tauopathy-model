PREFIX="/Volumes/BIOINFORMATICS/projects/hassan2022/output/merged/expectedCells/"
SUFFIX="/outs/filtered_feature_bc_matrix"
OUTS="/Volumes/BIOINFORMATICS/projects/hassan2022/output/merged/scrublet/"
SAMPLES=("control1" "CRN00224919" "CRN00224920" "P251L2merged" "CTRL2merged" "TauKO1merged" "LIB055588_CRN00233457" "LIB055588_CRN00233458")

~/anaconda3/bin/python scrublet_multi.py $PREFIX $SUFFIX $OUTS "${SAMPLES[@]}"
