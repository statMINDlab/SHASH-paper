# ABIDE analysis

This folder contains code and selected results for the ABIDE application
in the SHASH manuscript.

## Contents

- `scripts/ABIDE_direct_ICA_SHASHi_MCD075_q9973_analysis.R`:
  Performs ICA, SHASH-i transformation, and MCD-based outlier detection
  for both participants. It also generates QQ plots and robust-distance
  plots.

- `scripts/ABIDE1_mediandiff.R`:
  Creates temporal median and absolute median-difference NIfTI images
  for Subject 1. These images were viewed in MANGO for Figure 13.

- `results/`: Selected saved volume-level results and flagged-volume lists.

- `figures/`: Selected QQ plots, the robust-distance plot, and the
  assembled brain-image figure.

- `sessionInfo.txt`: R session information from the original analysis.

- `figures/figure13_source_images/`: Individual MANGO images and legends
  used to assemble Figure 13 in Adobe Illustrator.

## Main analysis

The script uses `fMRIscrub::Dat1` and `fMRIscrub::Dat2` as inputs.
It applies SHASH-i to ICA time courses and compares MCD robust distances
before and after transformation. MCD uses a subset fraction of 0.75 and
a GM14 scaled-F cutoff with probability 0.9973.

The main analysis requires R packages including `fMRIscrub`, `rrobot`,
`robustbase`, `CerioliOutlierDetection`, `dplyr`, `tibble`, `ggplot2`,
`ggrepel`, `scales`, and `tidyr`. See `sessionInfo.txt` for the original
R and package versions.

Run the script from the `ABIDE_analysis` directory:

```r
source("scripts/ABIDE_direct_ICA_SHASHi_MCD075_q9973_analysis.R")
```

The script creates its output folder under the current working directory.
Running it again may overwrite files in that output folder.

## Brain-image visualization

`scripts/ABIDE1_mediandiff.R` requires `RNifti`, the Subject 1 fMRI NIfTI
file, and its brain mask. Set `data_dir` in the script to the directory
containing `fwarest.nii.gz` and `brain_mask.nii.gz` before running it.

The generated images were viewed in MANGO and the figure was assembled
in Adobe Illustrator. `figures/ABIDE1_fig.pdf` is an assembled figure;
`figures/ABIDE_fig.ai` is an editable Illustrator file.