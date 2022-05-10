import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

args = sys.argv

input_dir_prefix = args[1]
input_dir_suffix = args[2]
out_dir = args[3]
samples = args[4:]

for sample in samples:
  input_dir = input_dir_prefix + sample + input_dir_suffix
  print("reading from: ", input_dir)
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
  print("Processing ", sample)
  scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
  doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)
  np.savetxt(out_dir + sample + "_srublet.score", doublet_scores)
  print("Showing predicted doublets")
  print(str(predicted_doublets))
  if predicted_doublets is None:
      predicted_doublets = scrub.call_doublets(threshold=6)
  np.savetxt(out_dir + sample + "_srublet.logic", predicted_doublets)

