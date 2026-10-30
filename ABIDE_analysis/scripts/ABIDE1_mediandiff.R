# =============================================================================
# ABIDE1 absolute median-difference 4D NIfTI for Mango visualization
# =============================================================================
# Goal:
#   Compute, for every fMRI volume t:
#
#       |Y_t(v) - median_s Y_s(v)|
#
#   and save the result as one 4D NIfTI file.
#
#   This allows navigation across volumes/time points in Mango.
#
# Outputs:
#   1. Full-brain masked 4D absolute median-difference image
#   2. Dat1 sagittal-slice-only 4D absolute median-difference image
#
# Important:
#   The Dat1 analysis is slice-level. Therefore, the sagittal-slice-only
#   output is the most directly aligned with the analyzed data object.
# =============================================================================

library(RNifti)

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

# Set this to the directory containing the Subject 1 fMRI and mask files.
data_dir <- "path/to/Pitt_0050048"

fmri_file <- file.path(data_dir, "fwarest.nii.gz")
mask_file <- file.path(data_dir, "brain_mask.nii.gz")

out_dir <- file.path(data_dir, "ABIDE1_meddiff_NifTi")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(fmri_file))
stopifnot(file.exists(mask_file))

# -----------------------------------------------------------------------------
# 2. Read fMRI and mask
# -----------------------------------------------------------------------------

fmri <- RNifti::readNifti(fmri_file)
mask <- RNifti::readNifti(mask_file)

fmri_arr <- as.array(fmri)          # 91 x 109 x 91 x 193
mask_arr <- as.array(mask) > 0      # use > 0, not == 1

stopifnot(length(dim(fmri_arr)) == 4)
stopifnot(all(dim(fmri_arr)[1:3] == dim(mask_arr)))

n_time <- dim(fmri_arr)[4]

cat("Image dimensions:\n")
print(dim(fmri_arr))

cat("Number of time volumes:", n_time, "\n")

# -----------------------------------------------------------------------------
# 3. Volumes to include in the 4D output
# -----------------------------------------------------------------------------
# Option A: include all fMRI volumes
# This creates a 4D file with the same number of time points as the original fMRI.

volumes_to_include <- seq_len(n_time)

# Option B: include only detected volumes
# Uncomment this if you want a shorter 4D file containing only detected volumes.
#
# volumes_to_include <- c(60, 61, 68, 69, 70, 71)

volumes_to_include <- sort(unique(as.integer(volumes_to_include)))

if (any(volumes_to_include < 1 | volumes_to_include > n_time)) {
  stop("At least one selected volume is outside the valid time range.")
}

cat("Number of volumes included in output:", length(volumes_to_include), "\n")

# -----------------------------------------------------------------------------
# 4. Identify the Dat1 sagittal slice
# -----------------------------------------------------------------------------
# Mango/world coordinate discussed for Dat1:
#   x = 2, y = 0, z = 18
#
# These are world coordinates in mm, not R array indices.
# Since Dat1 corresponds to a sagittal slice, we use x = 2 mm.

world_coord <- c(2, 0, 18)

A <- RNifti::xform(fmri)

voxel_coord_0based <- solve(A) %*% c(world_coord, 1)

voxel_index_R <- round(voxel_coord_0based[1:3]) + 1
voxel_index_R <- as.integer(voxel_index_R)

x_index <- voxel_index_R[1]

if (x_index < 1 || x_index > dim(mask_arr)[1]) {
  stop("Converted x-index is outside the image dimension.")
}

cat("Selected sagittal slice x-index in R:", x_index, "\n")
cat("Number of in-mask voxels in selected slice:", sum(mask_arr[x_index, , ]), "\n")

# -----------------------------------------------------------------------------
# 5. Create sagittal slice mask
# -----------------------------------------------------------------------------

slice_mask <- array(FALSE, dim = dim(mask_arr))
slice_mask[x_index, , ] <- mask_arr[x_index, , ]

if (sum(slice_mask) == 0) {
  stop("The selected sagittal slice mask contains zero voxels.")
}

# -----------------------------------------------------------------------------
# 6. Compute voxelwise temporal median image
# -----------------------------------------------------------------------------

cat("Computing voxelwise temporal median image...\n")

median_img <- apply(
  fmri_arr,
  c(1, 2, 3),
  median,
  na.rm = TRUE
)

# Save temporal median images for reference

median_img_full_masked <- median_img
median_img_full_masked[!mask_arr] <- 0

RNifti::writeNifti(
  RNifti::updateNifti(
    median_img_full_masked,
    template = fmri
  ),
  file.path(
    out_dir,
    "ABIDE1_temporal_median_full_masked.nii.gz"
  )
)

median_img_sagittal <- array(0, dim = dim(mask_arr))
median_img_sagittal[slice_mask] <- median_img[slice_mask]


RNifti::writeNifti(
  RNifti::updateNifti(
    median_img_sagittal,
    template = fmri
  ),
  file.path(
    out_dir,
    paste0(
      "ABIDE1_temporal_median_sagittal_xindex_",
      x_index,
      ".nii.gz"
    )
  )
)

# -----------------------------------------------------------------------------
# 7. Create 4D absolute median-difference arrays
# -----------------------------------------------------------------------------

cat("Computing 4D absolute median-difference images...\n")

spatial_dim <- dim(fmri_arr)[1:3]
n_out <- length(volumes_to_include)

abs_diff_full_4D <- array(
  0,
  dim = c(spatial_dim, n_out)
)

abs_diff_sagittal_4D <- array(
  0,
  dim = c(spatial_dim, n_out)
)

for (i in seq_along(volumes_to_include)) {
  
  vol <- volumes_to_include[i]
  
  cat("Processing original volume", vol, "as output frame", i, "\n")
  
  vol_img <- fmri_arr[, , , vol]
  
  abs_diff_img <- abs(
    vol_img - median_img
  )
  
  # Full-brain masked absolute difference
  abs_diff_full <- abs_diff_img
  abs_diff_full[!mask_arr] <- 0
  
  # Sagittal-slice-only absolute difference
  abs_diff_sagittal <- array(0, dim = spatial_dim)
  abs_diff_sagittal[slice_mask] <- abs_diff_img[slice_mask]
  
  abs_diff_full_4D[, , , i] <- abs_diff_full
  abs_diff_sagittal_4D[, , , i] <- abs_diff_sagittal
}

# -----------------------------------------------------------------------------
# 8. Save the 4D NIfTI files
# -----------------------------------------------------------------------------

full_4D_file <- file.path(
  out_dir,
  "ABIDE1_abs_minus_temporal_median_all_volumes_full_masked_4D.nii.gz"
)

sagittal_4D_file <- file.path(
  out_dir,
  paste0(
    "ABIDE1_abs_minus_temporal_median_all_volumes_sagittal_xindex_",
    x_index,
    "_4D.nii.gz"
  )
)

RNifti::writeNifti(
  RNifti::updateNifti(
    abs_diff_full_4D,
    template = fmri
  ),
  full_4D_file
)

RNifti::writeNifti(
  RNifti::updateNifti(
    abs_diff_sagittal_4D,
    template = fmri
  ),
  sagittal_4D_file
)

# -----------------------------------------------------------------------------
# 9. Save volume mapping
# -----------------------------------------------------------------------------
# This is important if you include only detected volumes.
# Mango frame 1 may correspond to original volume 60, etc.

volume_mapping <- data.frame(
  Mango_frame = seq_along(volumes_to_include),
  Original_volume = volumes_to_include
)

write.csv(
  volume_mapping,
  file = file.path(
    out_dir,
    "ABIDE1_4D_volume_mapping.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 10. Print outputs
# -----------------------------------------------------------------------------

cat("\nSaved 4D Mango visualization files to:\n")
cat(out_dir, "\n\n")

cat("Full-brain 4D absolute median-difference file:\n")
cat(full_4D_file, "\n\n")

cat("Sagittal-slice-only 4D absolute median-difference file:\n")
cat(sagittal_4D_file, "\n\n")

cat("Volume mapping saved as:\n")
cat(file.path(out_dir, "ABIDE1_4D_volume_mapping.csv"), "\n")

