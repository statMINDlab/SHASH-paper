# =============================================================================
# ABIDE DIRECT ICA ANALYSIS: SHASH-i QQ PLOTS AND MCD DISTANCES
# =============================================================================
#
# Workflow
#   1. ICA projection with fMRIscrub::pscrub(projection = "ICA")
#   2. Use the ICA time-course matrix M directly:
#        rows = volumes, columns = ICs
#
#   3. Univariate IC-wise analysis:
#        - compute robust-z scores for each original IC
#        - apply SHASH-i directly to each IC time course
#        - create QQ plots before and after SHASH-i
#
#   4. Multivariate detection using Minimum Covariance Determinant (MCD):
#        - robustbase::covMcd(alpha = 0.75)
#        - before SHASH-i: MCD applied directly to the original ICA matrix M
#        - after SHASH-i: MCD applied to the SHASH-i transformed ICA matrix
#
#      No additional robust-z standardization is required before MCD,
#      because the MCD covariance estimate accounts for variable scaling
#      and correlation when computing robust Mahalanobis distances.
#
#   5. MCD outlier rule:
#        - keep the MCD robust distances exactly as originally computed
#        - replace only the chi-square cutoff with the Green-Martin ("GM14")
#          scaled-F cutoff for MCD alpha = 0.75
#        - classify a volume as an outlier when RD^2 exceeds the
#          0.9973 GM14 scaled-F cutoff
#
# Important
#   - No PCA or ROBPCA is used in this script.
#   - SHASH-i is applied directly to each ICA time course.
#   - Robust-z standardization is used only for the univariate analysis,
#     not as preprocessing for MCD.
#   - QQ values outside [-plot_limit, plot_limit] are capped only for display.
#   - MCD robust distances are upper-capped only for visualization.
#   - MCD uses robustbase::covMcd() with alpha = 0.75.
#   - Raw, uncapped values are retained for all outlier decisions
#     and numerical summaries.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(fMRIscrub)
library(robustbase)
library(rrobot)
library(dplyr)
library(tibble)
library(ggplot2)

library(ggrepel)
library(scales)

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Package 'ggrepel' is required. Install it with install.packages('ggrepel').")
}


if (!requireNamespace("CerioliOutlierDetection", quietly = TRUE)) {
  stop(
    "Package 'CerioliOutlierDetection' is required for the GM14 scaled-F cutoff. ",
    "Install it with install.packages('CerioliOutlierDetection')."
  )
}


# -----------------------------------------------------------------------------
# 2. Reproducibility and settings
# -----------------------------------------------------------------------------

set.seed(20260728)

# Display range only. Analytical values are never capped.
plot_limit <- 10

# Number of facet columns in the all-IC QQ figures.
all_IC_columns <- 6

# Univariate score threshold for IC-level robust-z and SHASH-i flags.
univariate_threshold <- 3

# MCD settings.
# alpha = 0.75 is used for the MCD fit.
# No robust-z standardization is applied before the multivariate MCD step.
mcd_alpha <- 0.75

# Cutoff probability for squared robust distances.
# The MCD fit and robust distances remain unchanged from the original analysis.
# Only the reference cutoff is changed from chi-square to the Green-Martin
# ("GM14") finite-sample scaled-F approximation.
mcd_cutoff_probability <- 0.9973
mcd_significance_level <- 1 - mcd_cutoff_probability
mcd_cutoff_method <- "GM14"
mcd_cutoff_variant <- "cutoff.pred"

# The covariance used for the MCD robust distances is the same corrected
# final/reweighted covariance returned by robustbase::covMcd() in the original
# analysis. Therefore the cutoff function does not apply the MCD consistency
# correction a second time.
mcd_use_consistency_correction <- FALSE

# Plotting-only upper cap for MCD robust-distance figures.
# MCD distances are nonnegative, so only an upper cap is needed.
# The actual cap used in each dataset is:
#   max(mcd_plot_minimum_cap, mcd_plot_cutoff_multiplier * MCD_cutoff)
# Raw distances are never modified for detection or numerical summaries.
mcd_plot_minimum_cap <- 50
mcd_plot_cutoff_multiplier <- 5

root_out <- file.path(
  getwd(),
  "ABIDE_direct_ICA_SHASHi_MCD075_GM14F_q9973_analysis"
)

dir.create(
  root_out,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 3. Helper: SHASH-i fitting and transformation
# -----------------------------------------------------------------------------

fit_and_transform <- function(x) {
  
  x <- as.numeric(x)
  
  if (anyNA(x)) {
    stop("The input vector contains missing values.")
  }
  
  if (!all(is.finite(x))) {
    stop("The input vector contains non-finite values.")
  }
  
  fit <- rrobot::SHASH_out(
    x,
    thr1 = univariate_threshold,
    tail = "both",
    method_init = "isoplus",
    thr_iso = 0.6,
    thr = univariate_threshold
  )
  
  coefficients <- fit$SHASH_coef
  
  transformed <- rrobot::SHASH_to_normal(
    x = x,
    mu = coefficients$mu,
    sigma = coefficients$sigma,
    nu = coefficients$nu,
    tau = coefficients$tau
  )
  
  transformed <- as.numeric(transformed)
  
  out_idx <- fit$out_idx
  
  if (is.null(out_idx) || length(out_idx) == 0) {
    out_idx <- integer(0)
  } else {
    out_idx <- sort(unique(as.integer(out_idx)))
    out_idx <- out_idx[
      out_idx >= 1 &
        out_idx <= length(x)
    ]
  }
  
  list(
    original = x,
    transformed = transformed,
    out_idx = out_idx,
    fit = fit,
    coefficients = coefficients
  )
}


# -----------------------------------------------------------------------------
# 4. Helper: robust z-score
# -----------------------------------------------------------------------------

robust_z_score <- function(x) {
  
  x <- as.numeric(x)
  
  centre <- median(
    x,
    na.rm = TRUE
  )
  
  scale_value <- mad(
    x,
    center = centre,
    constant = 1.4826,
    na.rm = TRUE
  )
  
  if (!is.finite(scale_value) || scale_value == 0) {
    stop("The robust MAD is zero or non-finite.")
  }
  
  (x - centre) / scale_value
}


# -----------------------------------------------------------------------------
# 5. Helper: collapse IC numbers into a readable string
# -----------------------------------------------------------------------------

collapse_IC_numbers <- function(x) {
  
  x <- sort(unique(as.integer(x)))
  
  if (length(x) == 0) {
    return("")
  }
  
  paste(x, collapse = ", ")
}


# -----------------------------------------------------------------------------
# 6. Helper: add QQ coordinates and plotting-only capped values
# -----------------------------------------------------------------------------

prepare_QQ_coordinates <- function(qq_data, limit = 10) {
  
  label_limit <- limit - 0.4
  
  qq_data %>%
    group_by(
      Dataset,
      IC_number,
      Group
    ) %>%
    arrange(
      Value,
      .by_group = TRUE
    ) %>%
    mutate(
      Theoretical_quantile = qnorm(ppoints(n())),
      
      # Used only as the vertical plotting position.
      Plot_value = pmax(pmin(Value, limit), -limit),
      
      # Keeps labels just inside the panel border.
      Label_value = pmax(pmin(Value, label_limit), -label_limit),
      
      # Arrow indicates that the true score lies beyond the display range.
      Plot_label = case_when(
        Value > limit ~ paste0(Volume, "\u2191"),
        Value < -limit ~ paste0(Volume, "\u2193"),
        TRUE ~ as.character(Volume)
      ),
      
      Was_capped_for_plotting = abs(Value) > limit
    ) %>%
    ungroup()
}


# -----------------------------------------------------------------------------
# 7. Helper: MCD robust distances
# -----------------------------------------------------------------------------

compute_MCD_distances <- function(X, dataset_name, subject_label, stage_label) {
  
  X <- as.matrix(X)
  
  if (anyNA(X)) {
    stop("MCD input contains missing values for ", dataset_name, ", ", stage_label, ".")
  }
  
  if (!all(is.finite(X))) {
    stop("MCD input contains non-finite values for ", dataset_name, ", ", stage_label, ".")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (p >= n) {
    stop(
      "MCD requires fewer variables than observations. For ",
      dataset_name,
      ", p = ", p,
      " and n = ", n,
      "."
    )
  }
  
  # Minimum Covariance Determinant using robustbase.
  # By default covMcd() returns the final reweighted MCD estimates.
  mcd_fit <- robustbase::covMcd(
    X,
    alpha = mcd_alpha
  )
  
  mcd_center <- mcd_fit$center
  mcd_cov <- mcd_fit$cov
  mcd_h <- mcd_fit$quan
  
  # Robust squared Mahalanobis distances.
  robust_distance_squared <- mahalanobis(
    X,
    center = mcd_center,
    cov = mcd_cov
  )
  
  robust_distance <- sqrt(robust_distance_squared)
  
  # Green-Martin ("GM14") scaled-F cutoff.
  # IMPORTANT: the MCD fit, center, covariance, and robust distances above
  # are unchanged from the original chi-square analysis.
  cutoff_info <- CerioliOutlierDetection::hr05CutoffMvnormal(
    n.obs = n,
    p.dim = p,
    mcd.alpha = mcd_alpha,
    signif.alpha = mcd_significance_level,
    method = mcd_cutoff_method,
    use.consistency.correction = mcd_use_consistency_correction
  )
  
  distance_cutoff_squared <- cutoff_info[[mcd_cutoff_variant]]
  distance_cutoff <- sqrt(distance_cutoff_squared)
  
  tibble(
    Dataset = dataset_name,
    Subject = subject_label,
    Stage = stage_label,
    Volume = seq_len(n),
    Number_of_ICs = p,
    MCD_package = "robustbase::covMcd",
    MCD_alpha = mcd_alpha,
    MCD_h = mcd_h,
    MCD_distance = robust_distance,
    MCD_distance_squared = robust_distance_squared,
    MCD_cutoff_probability = mcd_cutoff_probability,
    MCD_cutoff_distribution = "scaled F",
    MCD_cutoff_method = "GM14 (Green-Martin extension)",
    MCD_cutoff_variant = mcd_cutoff_variant,
    MCD_cutoff_m_pred = cutoff_info$m.pred,
    MCD_cutoff_m_asy = cutoff_info$m.asy,
    MCD_cutoff_consistency_factor = cutoff_info$c.alpha,
    MCD_cutoff = distance_cutoff,
    MCD_cutoff_squared = distance_cutoff_squared,
    MCD_outlier = robust_distance_squared > distance_cutoff_squared
  )
}


# -----------------------------------------------------------------------------
# 8. Helper: add plotting-only capped MCD distances
# -----------------------------------------------------------------------------

prepare_MCD_plot_data <- function(mcd_data) {
  
  if (!all(c("MCD_distance", "MCD_cutoff") %in% names(mcd_data))) {
    stop("MCD data must contain MCD_distance and MCD_cutoff columns.")
  }
  
  if (!all(is.finite(mcd_data$MCD_distance))) {
    stop("MCD_distance contains non-finite values.")
  }
  
  if (!all(is.finite(mcd_data$MCD_cutoff))) {
    stop("MCD_cutoff contains non-finite values.")
  }
  
  mcd_plot_cap <- max(
    mcd_plot_minimum_cap,
    mcd_plot_cutoff_multiplier * max(mcd_data$MCD_cutoff, na.rm = TRUE)
  )
  
  mcd_data %>%
    mutate(
      MCD_distance_plot = pmin(MCD_distance, mcd_plot_cap),
      MCD_distance_was_capped_for_plotting = MCD_distance > mcd_plot_cap,
      MCD_plot_cap = mcd_plot_cap,
      MCD_plot_label = if_else(
        MCD_distance_was_capped_for_plotting,
        paste0(Volume, "\u2191"),
        as.character(Volume)
      )
    )
}


# -----------------------------------------------------------------------------
# 9. Plot function: all IC QQ plots for one dataset
# -----------------------------------------------------------------------------

make_all_IC_QQ_plot <- function(
    qq_data,
    subject_label,
    number_ICs,
    limit = 10,
    ncol_facets = 6
) {
  
  IC_levels <- paste0(
    "IC",
    seq_len(number_ICs)
  )
  
  plot_data <- qq_data %>%
    mutate(
      IC = factor(
        IC,
        levels = IC_levels
      )
    )
  
  ggplot(
    plot_data,
    aes(
      x = Theoretical_quantile,
      y = Plot_value,
      color = Group
    )
  ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -univariate_threshold,
      ymax = univariate_threshold,
      fill = "gray85",
      alpha = 0.6
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "black",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = c(-univariate_threshold, univariate_threshold),
      color = "black",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    geom_point(
      size = 0.65,
      alpha = 0.60
    ) +
    ggrepel::geom_text_repel(
      data = plot_data %>%
        filter(
          Is_outlier,
          Group == "robust-z ICA (before)"
        ),
      aes(
        y = Label_value,
        label = Plot_label
      ),
      size = 1.9,
      fontface = "bold",
      direction = "y",
      nudge_x = -0.08,
      hjust = 1,
      box.padding = 0.15,
      point.padding = 0.08,
      min.segment.length = 0,
      segment.color = "grey55",
      segment.alpha = 0.6,
      segment.size = 0.2,
      max.overlaps = Inf,
      max.time = 2,
      max.iter = 20000,
      seed = 20260728,
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = plot_data %>%
        filter(
          Is_outlier,
          Group == "SHASH-i ICA (after)"
        ),
      aes(
        y = Label_value,
        label = Plot_label
      ),
      size = 1.9,
      fontface = "bold",
      direction = "y",
      nudge_x = 0.08,
      hjust = 0,
      box.padding = 0.15,
      point.padding = 0.08,
      min.segment.length = 0,
      segment.color = "grey55",
      segment.alpha = 0.6,
      segment.size = 0.2,
      max.overlaps = Inf,
      max.time = 2,
      max.iter = 20000,
      seed = 20260728,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "robust-z ICA (before)" = "steelblue",
        "SHASH-i ICA (after)" = "tomato"
      )
    ) +
    facet_wrap(
      ~ IC,
      ncol = ncol_facets,
      scales = "fixed"
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.10, 0.10))
    ) +
    scale_y_continuous(
      limits = c(-limit, limit),
      breaks = seq(-limit, limit, by = 2),
      expand = expansion(mult = 0.02)
    ) +
    labs(
      title = paste0(subject_label, ": QQ plots for all ICA time courses"),
      subtitle = paste0(
        number_ICs,
        " ICs; SHASH-i applied directly to ICA time courses; values outside +/-",
        limit,
        " are capped only for display"
      ),
      x = "Theoretical quantiles",
      y = "Sample quantiles",
      color = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.7
      ),
      panel.spacing = grid::unit(0.9, "lines"),
      strip.text = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 10),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9)
    )
}


# -----------------------------------------------------------------------------
# 10. Plot function: MCD robust distances
# -----------------------------------------------------------------------------

make_MCD_distance_plot <- function(mcd_data, subject_label) {
  
  plot_data <- mcd_data
  
  if (!"MCD_distance_plot" %in% names(plot_data)) {
    plot_data <- prepare_MCD_plot_data(plot_data)
  }
  
  cutoff_data <- plot_data %>%
    distinct(
      Dataset,
      Subject,
      Stage,
      MCD_cutoff,
      MCD_plot_cap
    )
  
  plot_cap <- max(
    plot_data$MCD_plot_cap,
    na.rm = TRUE
  )
  
  y_upper_limit <- plot_cap + max(
    1,
    0.06 * plot_cap
  )
  
  ggplot(
    plot_data,
    aes(
      x = Volume,
      y = MCD_distance_plot,
      color = Stage
    )
  ) +
    geom_line(
      linewidth = 0.55,
      alpha = 0.85
    ) +
    geom_point(
      data = plot_data %>% filter(MCD_outlier),
      size = 1.7,
      alpha = 0.90,
      show.legend = FALSE
    ) +
    geom_hline(
      data = cutoff_data,
      aes(yintercept = MCD_cutoff),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.45,
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = plot_data %>% filter(MCD_outlier),
      aes(
        y = MCD_distance_plot,
        label = MCD_plot_label
      ),
      size = 2.5,
      fontface = "bold",
      max.overlaps = Inf,
      seed = 20260728,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "MCD on original ICA before SHASH-i" = "steelblue",
        "MCD on SHASH-i transformed ICA" = "tomato"
      )
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 10),
      minor_breaks = scales::pretty_breaks(n = 20),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_continuous(
      limits = c(0, y_upper_limit),
      breaks = pretty(c(0, plot_cap), n = 8),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste0(subject_label, ": MCD robust distances in ICA space"),
      subtitle = paste0(
        "robustbase::covMcd(alpha = ",
        format(mcd_alpha, nsmall = 2),
        "); dashed line = GM14 scaled-F q",
        mcd_cutoff_probability,
        " cutoff; distances above ",
        round(plot_cap, 1),
        " are capped only for display"
      ),
      x = "Volume number",
      y = "MCD robust distance",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.6
      ),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    )
}



# -----------------------------------------------------------------------------
# 11. ABIDE datasets
# -----------------------------------------------------------------------------

datasets <- list(
  ABIDE1 = fMRIscrub::Dat1,
  ABIDE2 = fMRIscrub::Dat2
)

dataset_labels <- c(
  ABIDE1 = "Subject 1 (high noise)",
  ABIDE2 = "Subject 2 (low noise)"
)


# -----------------------------------------------------------------------------
# 12. Storage objects
# -----------------------------------------------------------------------------

analysis_results <- list()

dimension_summary_rows <- list()
timing_rows <- list()

IC_outlier_summary_rows <- list()
before_outlier_detail_rows <- list()
after_outlier_detail_rows <- list()
before_after_comparison_rows <- list()
visual_inspection_summary_rows <- list()
final_univariate_union_summary_rows <- list()

QQ_rows_all <- list()
outlier_consistency_rows <- list()
MCD_distance_rows <- list()
MCD_volume_summary_rows <- list()

all_IC_plots <- list()
MCD_distance_plots <- list()


# -----------------------------------------------------------------------------
# 13. Complete direct-ICA analysis for each dataset
# -----------------------------------------------------------------------------

for (dataset_name in names(datasets)) {
  
  cat(
    "\n========================================\n",
    "Processing ",
    dataset_name,
    "\n========================================\n",
    sep = ""
  )
  
  dataset_output_directory <- file.path(
    root_out,
    dataset_name
  )
  
  dir.create(
    dataset_output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dat <- datasets[[dataset_name]]
  subject_label <- dataset_labels[[dataset_name]]
  
  
  # ---------------------------------------------------------------------------
  # 12.1 ICA projection
  # ---------------------------------------------------------------------------
  
  ICA_time <- system.time({
    ICA_result <- fMRIscrub::pscrub(
      dat,
      projection = "ICA"
    )
  })
  
  M <- as.matrix(ICA_result$ICA$M)
  
  number_ICs <- ncol(M)
  number_volumes <- nrow(M)
  
  colnames(M) <- paste0("IC", seq_len(number_ICs))
  rownames(M) <- seq_len(number_volumes)
  
  cat("Number of volumes:", number_volumes, "\n")
  cat("Number of ICA components:", number_ICs, "\n")
  
  dimension_summary_rows[[dataset_name]] <- tibble(
    Dataset = dataset_name,
    Subject = subject_label,
    Number_of_volumes = number_volumes,
    Number_of_ICA_components = number_ICs,
    PCA_used = FALSE,
    Analysis_space = "ICA time courses directly",
    MCD_package = "robustbase::covMcd",
    MCD_alpha = mcd_alpha,
    MCD_cutoff_probability = mcd_cutoff_probability,
    MCD_cutoff_distribution = "scaled F",
    MCD_cutoff_method = "GM14 (Green-Martin extension)"
  )
  
  
  # ---------------------------------------------------------------------------
  # 12.2 Initialize IC-level objects
  # ---------------------------------------------------------------------------
  
  before_outlier_matrix <- matrix(
    FALSE,
    nrow = number_volumes,
    ncol = number_ICs
  )
  
  after_outlier_matrix <- matrix(
    FALSE,
    nrow = number_volumes,
    ncol = number_ICs
  )
  
  colnames(before_outlier_matrix) <- paste0("IC", seq_len(number_ICs))
  colnames(after_outlier_matrix) <- paste0("IC", seq_len(number_ICs))
  
  rownames(before_outlier_matrix) <- seq_len(number_volumes)
  rownames(after_outlier_matrix) <- seq_len(number_volumes)
  
  # Robust-z scores are retained for the univariate IC-wise analysis and QQ plots.
  # They are NOT used as input to the multivariate MCD analysis.
  ICA_robz <- matrix(
    NA_real_,
    nrow = number_volumes,
    ncol = number_ICs
  )
  
  ICA_shash <- matrix(
    NA_real_,
    nrow = number_volumes,
    ncol = number_ICs
  )
  
  colnames(ICA_robz) <- paste0("IC", seq_len(number_ICs))
  colnames(ICA_shash) <- paste0("IC", seq_len(number_ICs))
  rownames(ICA_robz) <- seq_len(number_volumes)
  rownames(ICA_shash) <- seq_len(number_volumes)
  
  IC_SHASH_results <- vector(
    mode = "list",
    length = number_ICs
  )
  
  names(IC_SHASH_results) <- paste0("IC", seq_len(number_ICs))
  
  IC_summary_dataset <- list()
  IC_volume_dataset <- list()
  QQ_dataset <- list()
  consistency_dataset <- list()
  
  SHASH_elapsed_total <- 0
  
  
  # ---------------------------------------------------------------------------
  # 12.3 Analyze every ICA time course directly
  # ---------------------------------------------------------------------------
  
  for (IC_number in seq_len(number_ICs)) {
    
    IC_name <- paste0("IC", IC_number)
    IC_values <- as.numeric(M[, IC_number])
    volume_number <- seq_along(IC_values)
    
    # Before transformation: robust-z scores of raw ICA time course
    before_scores <- robust_z_score(IC_values)
    ICA_robz[, IC_number] <- before_scores
    
    before_outlier_indices <- which(abs(before_scores) > univariate_threshold)
    
    if (length(before_outlier_indices) > 0) {
      before_outlier_matrix[before_outlier_indices, IC_number] <- TRUE
    }
    
    # After transformation: SHASH-i applied directly to the ICA time course
    SHASH_time <- system.time({
      SHASH_result <- fit_and_transform(IC_values)
    })
    
    SHASH_elapsed_total <-
      SHASH_elapsed_total +
      unname(SHASH_time["elapsed"])
    
    IC_SHASH_results[[IC_number]] <- SHASH_result
    
    after_scores <- SHASH_result$transformed
    ICA_shash[, IC_number] <- after_scores
    
    after_outlier_indices <- SHASH_result$out_idx
    
    if (length(after_outlier_indices) > 0) {
      after_outlier_matrix[after_outlier_indices, IC_number] <- TRUE
    }
    
    # Check whether SHASH_out()$out_idx agrees with final +/- threshold scores.
    transformed_threshold_indices <- which(abs(after_scores) > univariate_threshold)
    
    consistency_dataset[[IC_number]] <- tibble(
      Dataset = dataset_name,
      IC = IC_number,
      SHASH_out_indices = paste(after_outlier_indices, collapse = ", "),
      Transformed_values_outside_threshold = paste(transformed_threshold_indices, collapse = ", "),
      Same_indices = setequal(after_outlier_indices, transformed_threshold_indices)
    )
    
    # Per-IC counts and volume lists
    IC_summary_dataset[[IC_number]] <- tibble(
      Dataset = dataset_name,
      Subject = subject_label,
      IC = IC_number,
      Number_before_outliers = length(before_outlier_indices),
      Before_outlier_volumes = if (length(before_outlier_indices) == 0) {
        ""
      } else {
        paste(before_outlier_indices, collapse = ", ")
      },
      Number_after_outliers = length(after_outlier_indices),
      After_outlier_volumes = if (length(after_outlier_indices) == 0) {
        ""
      } else {
        paste(after_outlier_indices, collapse = ", ")
      }
    )
    
    # Complete volume-level record for this IC
    IC_volume_dataset[[IC_number]] <- tibble(
      Dataset = dataset_name,
      Subject = subject_label,
      IC = IC_number,
      Volume = volume_number,
      Original_ICA_score = IC_values,
      Before_score = before_scores,
      Absolute_before_score = abs(before_scores),
      Before_direction = if_else(before_scores >= 0, "positive", "negative"),
      Flagged_before = volume_number %in% before_outlier_indices,
      After_score = after_scores,
      Absolute_after_score = abs(after_scores),
      After_direction = if_else(after_scores >= 0, "positive", "negative"),
      Flagged_after = volume_number %in% after_outlier_indices
    ) %>%
      mutate(
        Status = case_when(
          Flagged_before & Flagged_after ~ "flagged before and after",
          Flagged_before & !Flagged_after ~ "flagged before only",
          !Flagged_before & Flagged_after ~ "flagged after only",
          TRUE ~ "not flagged"
        )
      )
    
    # QQ data before transformation
    before_QQ <- tibble(
      Dataset = dataset_name,
      Subject = subject_label,
      IC_number = IC_number,
      IC = IC_name,
      Group = "robust-z ICA (before)",
      Volume = volume_number,
      Value = before_scores,
      Is_outlier = volume_number %in% before_outlier_indices
    )
    
    # QQ data after transformation
    after_QQ <- tibble(
      Dataset = dataset_name,
      Subject = subject_label,
      IC_number = IC_number,
      IC = IC_name,
      Group = "SHASH-i ICA (after)",
      Volume = volume_number,
      Value = after_scores,
      Is_outlier = volume_number %in% after_outlier_indices
    )
    
    QQ_dataset[[IC_number]] <- bind_rows(before_QQ, after_QQ)
  }
  
  
  # ---------------------------------------------------------------------------
  # 12.4 MCD robust distances in direct ICA space
  # ---------------------------------------------------------------------------
  
  MCD_time <- system.time({
    mcd_before_current <- compute_MCD_distances(
      X = M,
      dataset_name = dataset_name,
      subject_label = subject_label,
      stage_label = "MCD on original ICA before SHASH-i"
    )
    
    mcd_after_current <- compute_MCD_distances(
      X = ICA_shash,
      dataset_name = dataset_name,
      subject_label = subject_label,
      stage_label = "MCD on SHASH-i transformed ICA"
    )
  })
  
  MCD_distances_current <- bind_rows(
    mcd_before_current,
    mcd_after_current
  ) %>%
    prepare_MCD_plot_data()
  
  
  # ---------------------------------------------------------------------------
  # 12.5 Construct dataset-level tables
  # ---------------------------------------------------------------------------
  
  IC_outlier_summary_current <- bind_rows(IC_summary_dataset)
  IC_volume_all_current <- bind_rows(IC_volume_dataset)
  
  before_outlier_details_current <-
    IC_volume_all_current %>%
    filter(Flagged_before) %>%
    select(
      Dataset,
      Subject,
      IC,
      Volume,
      Original_ICA_score,
      Before_score,
      Absolute_before_score,
      Before_direction
    ) %>%
    arrange(IC, Volume)
  
  after_outlier_details_current <-
    IC_volume_all_current %>%
    filter(Flagged_after) %>%
    select(
      Dataset,
      Subject,
      IC,
      Volume,
      Original_ICA_score,
      After_score,
      Absolute_after_score,
      After_direction
    ) %>%
    arrange(IC, Volume)
  
  before_after_comparison_current <-
    IC_volume_all_current %>%
    filter(Flagged_before | Flagged_after) %>%
    select(
      Dataset,
      Subject,
      IC,
      Volume,
      Original_ICA_score,
      Before_score,
      Absolute_before_score,
      Before_direction,
      Flagged_before,
      After_score,
      Absolute_after_score,
      After_direction,
      Flagged_after,
      Status
    ) %>%
    arrange(IC, Volume)
  
  # Unique-volume table for visual inspection of brain images.
  visual_inspection_summary_current <-
    IC_volume_all_current %>%
    filter(Flagged_before | Flagged_after) %>%
    group_by(
      Dataset,
      Subject,
      Volume
    ) %>%
    summarise(
      Flagged_before_any = any(Flagged_before),
      Number_of_ICs_before = sum(Flagged_before),
      ICs_before = collapse_IC_numbers(IC[Flagged_before]),
      Maximum_absolute_before_score = if (any(Flagged_before)) {
        max(Absolute_before_score[Flagged_before], na.rm = TRUE)
      } else {
        NA_real_
      },
      Flagged_after_any = any(Flagged_after),
      Number_of_ICs_after = sum(Flagged_after),
      ICs_after = collapse_IC_numbers(IC[Flagged_after]),
      Maximum_absolute_after_score = if (any(Flagged_after)) {
        max(Absolute_after_score[Flagged_after], na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      Volume_status = case_when(
        Flagged_before_any & Flagged_after_any ~ "flagged before and after",
        Flagged_before_any & !Flagged_after_any ~ "flagged before only",
        !Flagged_before_any & Flagged_after_any ~ "flagged after only",
        TRUE ~ "not flagged"
      )
    ) %>%
    arrange(Volume)
  
  # Final direct-ICA univariate detection set: union across ICs after SHASH-i.
  final_univariate_union_summary_current <-
    after_outlier_details_current %>%
    group_by(
      Dataset,
      Subject,
      Volume
    ) %>%
    summarise(
      Number_of_ICs_flagging_volume = n_distinct(IC),
      ICs = collapse_IC_numbers(IC),
      Maximum_absolute_score = max(Absolute_after_score, na.rm = TRUE),
      Most_extreme_score = After_score[which.max(Absolute_after_score)][1],
      Direction_of_most_extreme_score = After_direction[which.max(Absolute_after_score)][1],
      .groups = "drop"
    ) %>%
    arrange(Volume)
  
  # MCD volume-level summary.
  mcd_before_volumes <- mcd_before_current %>%
    filter(MCD_outlier) %>%
    pull(Volume)
  
  mcd_after_volumes <- mcd_after_current %>%
    filter(MCD_outlier) %>%
    pull(Volume)
  
  before_union_volumes <- which(rowSums(before_outlier_matrix) > 0)
  after_union_volumes <- which(rowSums(after_outlier_matrix) > 0)
  
  MCD_volume_summary_current <- tibble(
    Dataset = dataset_name,
    Subject = subject_label,
    Volume = seq_len(number_volumes),
    Univariate_ICA_before_any = seq_len(number_volumes) %in% before_union_volumes,
    Univariate_ICA_after_any = seq_len(number_volumes) %in% after_union_volumes,
    MCD_before_outlier = seq_len(number_volumes) %in% mcd_before_volumes,
    MCD_after_outlier = seq_len(number_volumes) %in% mcd_after_volumes
  ) %>%
    left_join(
      mcd_before_current %>%
        select(Volume, MCD_distance_before = MCD_distance),
      by = "Volume"
    ) %>%
    left_join(
      mcd_after_current %>%
        select(Volume, MCD_distance_after = MCD_distance),
      by = "Volume"
    ) %>%
    arrange(Volume)
  
  # QQ coordinates and display-capped values.
  QQ_current <- bind_rows(QQ_dataset) %>%
    mutate(
      Group = factor(
        Group,
        levels = c(
          "robust-z ICA (before)",
          "SHASH-i ICA (after)"
        )
      )
    ) %>%
    prepare_QQ_coordinates(limit = plot_limit)
  
  
  # ---------------------------------------------------------------------------
  # 12.6 Store dataset results
  # ---------------------------------------------------------------------------
  
  IC_outlier_summary_rows[[dataset_name]] <- IC_outlier_summary_current
  before_outlier_detail_rows[[dataset_name]] <- before_outlier_details_current
  after_outlier_detail_rows[[dataset_name]] <- after_outlier_details_current
  before_after_comparison_rows[[dataset_name]] <- before_after_comparison_current
  visual_inspection_summary_rows[[dataset_name]] <- visual_inspection_summary_current
  final_univariate_union_summary_rows[[dataset_name]] <- final_univariate_union_summary_current
  QQ_rows_all[[dataset_name]] <- QQ_current
  outlier_consistency_rows[[dataset_name]] <- bind_rows(consistency_dataset)
  MCD_distance_rows[[dataset_name]] <- MCD_distances_current
  MCD_volume_summary_rows[[dataset_name]] <- MCD_volume_summary_current
  
  timing_rows[[dataset_name]] <- tibble(
    Dataset = dataset_name,
    ICA_elapsed_seconds = unname(ICA_time["elapsed"]),
    SHASH_all_ICs_elapsed_seconds = SHASH_elapsed_total,
    MCD_elapsed_seconds = unname(MCD_time["elapsed"])
  )
  
  
  # ---------------------------------------------------------------------------
  # 12.7 Make and save all-IC QQ figure
  # ---------------------------------------------------------------------------
  
  plot_all_ICs <- make_all_IC_QQ_plot(
    qq_data = QQ_current,
    subject_label = subject_label,
    number_ICs = number_ICs,
    limit = plot_limit,
    ncol_facets = all_IC_columns
  )
  
  all_IC_plots[[dataset_name]] <- plot_all_ICs
  
  print(plot_all_ICs)
  
  number_of_facet_rows <- ceiling(number_ICs / all_IC_columns)
  figure_height <- max(6, 2.7 * number_of_facet_rows + 2)
  
  ggsave(
    filename = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_all_ICs_QQ_direct_ICA_SHASHi.png")
    ),
    plot = plot_all_ICs,
    width = 14,
    height = figure_height,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_all_ICs_QQ_direct_ICA_SHASHi.pdf")
    ),
    plot = plot_all_ICs,
    width = 14,
    height = figure_height,
    bg = "white"
  )
  
  
  # ---------------------------------------------------------------------------
  # 12.8 Make and save MCD robust-distance figure
  # ---------------------------------------------------------------------------
  
  plot_MCD_distances <- make_MCD_distance_plot(
    mcd_data = MCD_distances_current,
    subject_label = subject_label
  )
  
  MCD_distance_plots[[dataset_name]] <- plot_MCD_distances
  
  print(plot_MCD_distances)
  
  ggsave(
    filename = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_MCD_robust_distances_direct_ICA.png")
    ),
    plot = plot_MCD_distances,
    width = 10,
    height = 5.5,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_MCD_robust_distances_direct_ICA.pdf")
    ),
    plot = plot_MCD_distances,
    width = 10,
    height = 5.5,
    bg = "white"
  )
  
  
  # ---------------------------------------------------------------------------
  # 12.9 Save dataset-specific data objects as RDS
  # ---------------------------------------------------------------------------
  
  saveRDS(
    IC_outlier_summary_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_IC_outlier_summary_before_after.rds")
    )
  )
  
  saveRDS(
    before_outlier_details_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_before_transformation_outliers_all_ICs.rds")
    )
  )
  
  saveRDS(
    after_outlier_details_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_after_transformation_outliers_all_ICs.rds")
    )
  )
  
  saveRDS(
    before_after_comparison_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_before_after_IC_volume_comparison.rds")
    )
  )
  
  saveRDS(
    visual_inspection_summary_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_visual_inspection_volume_summary_direct_ICA.rds")
    )
  )
  
  saveRDS(
    final_univariate_union_summary_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_final_SHASHi_union_outlier_summary_direct_ICA.rds")
    )
  )
  
  saveRDS(
    MCD_distances_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_MCD_distances_direct_ICA.rds")
    )
  )
  
  saveRDS(
    MCD_volume_summary_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_MCD_volume_summary_direct_ICA.rds")
    )
  )
  
  before_outlier_matrix_df <- as.data.frame(before_outlier_matrix) %>%
    rownames_to_column(var = "Volume")
  
  after_outlier_matrix_df <- as.data.frame(after_outlier_matrix) %>%
    rownames_to_column(var = "Volume")
  
  saveRDS(
    before_outlier_matrix_df,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_before_outlier_matrix_ICs.rds")
    )
  )
  
  saveRDS(
    after_outlier_matrix_df,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_after_outlier_matrix_ICs.rds")
    )
  )
  
  saveRDS(
    QQ_current,
    file = file.path(
      dataset_output_directory,
      paste0(dataset_name, "_QQ_plot_data_all_ICs.rds")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # 12.10 Save complete fitted objects
  # ---------------------------------------------------------------------------
  
  analysis_results[[dataset_name]] <- list(
    ICA_result = ICA_result,
    M = M,
    PCA_used = FALSE,
    number_ICs = number_ICs,
    ICA_robz = ICA_robz,
    ICA_shash = ICA_shash,
    IC_SHASH_results = IC_SHASH_results,
    before_outlier_matrix = before_outlier_matrix,
    after_outlier_matrix = after_outlier_matrix,
    before_union_volumes = before_union_volumes,
    after_union_volumes = after_union_volumes,
    MCD_distances = MCD_distances_current,
    MCD_volume_summary = MCD_volume_summary_current,
    IC_outlier_summary = IC_outlier_summary_current,
    before_outlier_details = before_outlier_details_current,
    after_outlier_details = after_outlier_details_current,
    before_after_comparison = before_after_comparison_current,
    visual_inspection_summary = visual_inspection_summary_current,
    final_univariate_union_summary = final_univariate_union_summary_current,
    QQ_data = QQ_current,
    all_IC_QQ_plot = plot_all_ICs,
    MCD_distance_plot = plot_MCD_distances
  )
}


# -----------------------------------------------------------------------------
# 13. Combine results across ABIDE1 and ABIDE2
# -----------------------------------------------------------------------------

dimension_summary <- bind_rows(dimension_summary_rows)
timing_summary <- bind_rows(timing_rows)
IC_outlier_summary <- bind_rows(IC_outlier_summary_rows)
before_outlier_details <- bind_rows(before_outlier_detail_rows)
after_outlier_details <- bind_rows(after_outlier_detail_rows)
before_after_comparison <- bind_rows(before_after_comparison_rows)
visual_inspection_volume_summary <- bind_rows(visual_inspection_summary_rows)
final_univariate_union_summary <- bind_rows(final_univariate_union_summary_rows)
outlier_consistency_check <- bind_rows(outlier_consistency_rows)
QQ_data_all <- bind_rows(QQ_rows_all)
MCD_distances_all <- bind_rows(MCD_distance_rows)
MCD_volume_summary_all <- bind_rows(MCD_volume_summary_rows)


# Exactly the labels used in the QQ plots, before and after, for all ICs.
QQ_outlier_labels_all_ICs <- QQ_data_all %>%
  filter(Is_outlier) %>%
  select(
    Dataset,
    Subject,
    IC_number,
    IC,
    Group,
    Volume,
    Value,
    Theoretical_quantile,
    Plot_value,
    Label_value,
    Plot_label,
    Was_capped_for_plotting
  ) %>%
  arrange(Dataset, IC_number, Group, Volume)


# Record every score that was capped only for display.
values_capped_for_plotting <- QQ_data_all %>%
  filter(Was_capped_for_plotting) %>%
  select(
    Dataset,
    Subject,
    IC_number,
    IC,
    Group,
    Volume,
    Value,
    Plot_value,
    Plot_label
  ) %>%
  arrange(Dataset, IC_number, Group, Volume)


# Combined MCD outlier list.
MCD_outlier_list <- MCD_distances_all %>%
  filter(MCD_outlier) %>%
  arrange(Dataset, Stage, Volume)


# Record every MCD distance that was capped only for display.
MCD_values_capped_for_plotting <- MCD_distances_all %>%
  filter(MCD_distance_was_capped_for_plotting) %>%
  select(
    Dataset,
    Subject,
    Stage,
    Volume,
    Number_of_ICs,
    MCD_package,
    MCD_alpha,
    MCD_h,
    MCD_distance,
    MCD_distance_plot,
    MCD_plot_cap,
    MCD_cutoff,
    MCD_outlier,
    MCD_plot_label
  ) %>%
  arrange(Dataset, Stage, Volume)


# Compact comparison of direct-ICA univariate SHASH-i and MCD-after-SHASH-i.
combined_detection_summary <- MCD_volume_summary_all %>%
  mutate(
    Final_status = case_when(
      Univariate_ICA_after_any & MCD_after_outlier ~ "univariate SHASH-i and MCD",
      Univariate_ICA_after_any & !MCD_after_outlier ~ "univariate SHASH-i only",
      !Univariate_ICA_after_any & MCD_after_outlier ~ "MCD only",
      TRUE ~ "not flagged"
    )
  )


# =============================================================================
# Manuscript QQ plots: first six ICs only
# =============================================================================
# These figures use the already-computed QQ_data_all object.
# No ICA or SHASH-i analysis is rerun.
#
# Display choices:
#   - IC1--IC6 only
#   - 3 columns x 2 rows
#   - no volume-number labels
#   - no ggrepel connecting segments
#   - same robust-z before / SHASH-i after values used in the all-IC QQ figures
# =============================================================================

QQ_first6_all <- QQ_data_all %>%
  filter(IC_number %in% 1:6) %>%
  mutate(
    IC = factor(
      IC,
      levels = paste0("IC", 1:6)
    )
  )

saveRDS(
  QQ_first6_all,
  file.path(root_out, "ABIDE_direct_ICA_first6_ICs_QQ_plot_data.rds")
)

first6_IC_plots <- list()

for (dataset_name in c("ABIDE1", "ABIDE2")) {
  
  subject_label <- dataset_labels[[dataset_name]]
  
  QQ_first6_current <- QQ_first6_all %>%
    filter(Dataset == dataset_name)
  
  plot_first6 <- ggplot(
    QQ_first6_current,
    aes(
      x = Theoretical_quantile,
      y = Plot_value,
      color = Group
    )
  ) +
    
    # Gray region corresponding to the non-outlier range
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -univariate_threshold,
      ymax = univariate_threshold,
      fill = "gray85",
      alpha = 0.6
    ) +
    
    # Standard Normal QQ reference line
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "black",
      linewidth = 0.4
    ) +
    
    # Outlier thresholds
    geom_hline(
      yintercept = c(
        -univariate_threshold,
        univariate_threshold
      ),
      color = "black",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    
    # QQ points
    geom_point(
      size = 1.4,
      alpha = 0.7
    ) +
    
    # Consistent before/after colors and legend labels
    scale_color_manual(
      values = c(
        "robust-z ICA (before)" = "steelblue",
        "SHASH-i ICA (after)" = "tomato"
      ),
      breaks = c(
        "robust-z ICA (before)",
        "SHASH-i ICA (after)"
      ),
      labels = c(
        "without normality transformation",
        "with normality transformation"
      )
    ) +
    
    # First six independent components
    facet_wrap(
      ~ IC,
      ncol = 6,
      scales = "fixed"
    ) +
    
    # X-axis
    scale_x_continuous(
      breaks = seq(-3, 3, by = 1),
      expand = expansion(mult = 0)
    ) +
    
    # Y-axis
    scale_y_continuous(
      breaks = seq(
        -plot_limit,
        plot_limit,
        by = 2
      ),
      expand = expansion(mult = 0.02)
    ) +
    
    # Display limits
    coord_cartesian(
      xlim = c(-3, 3),
      ylim = c(-plot_limit, plot_limit)
    ) +
    
    # Labels
    labs(
      title = paste0(
        subject_label,
        ": QQ plots for the first six ICA time courses"
      ),
      subtitle = paste0(
        "IC1--IC6; SHASH-i applied directly to ICA time courses; ",
        "values outside +/-",
        plot_limit,
        " are capped only for display"
      ),
      x = "Theoretical quantiles",
      y = "Sample quantiles",
      color = NULL
    ) +
    
    theme_minimal(base_size = 14) +
    
    guides(
      color = guide_legend(
        override.aes = list(
          size = 4,
          alpha = 1
        )
      )
    ) +
    
    theme(
      # Remove internal grid
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Border around each IC panel
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.7
      ),
      
      # Space between IC panels
      panel.spacing = grid::unit(
        0.9,
        "lines"
      ),
      
      # IC labels
      strip.text = element_text(
        face = "bold",
        size = 14
      ),
      
      # Axis numbers
      axis.text.x = element_text(
        size = 14
      ),
      axis.text.y = element_text(
        size = 14
      ),
      
      # Axis titles
      axis.title.x = element_text(
        size = 14
      ),
      axis.title.y = element_text(
        size = 14
      ),
      
      # Legend
      legend.position = "bottom",
      legend.text = element_text(
        size = 12
      ),
      
      # Main title
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      
      # Subtitle
      plot.subtitle = element_text(
        size = 14
      )
    )
  
  
  # Store plot
  first6_IC_plots[[dataset_name]] <- plot_first6
  
  # Display plot
  print(plot_first6)
  
  
  # -----------------------------------------------------------------------
  # Save PNG
  # -----------------------------------------------------------------------
  
  ggsave(
    filename = file.path(
      root_out,
      dataset_name,
      paste0(
        dataset_name,
        "_6ICs_QQ.png"
      )
    ),
    plot = plot_first6,
    width = 14,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  
  # -----------------------------------------------------------------------
  # Save PDF
  # -----------------------------------------------------------------------
  
  ggsave(
    filename = file.path(
      root_out,
      dataset_name,
      paste0(
        dataset_name,
        "_first6_ICs_QQ_direct_ICA_SHASHi_no_volume_labels.pdf"
      )
    ),
    plot = plot_first6,
    width = 14,
    height = 5,
    bg = "white"
  )
}


# Analysis settings saved explicitly for reproducibility and later comparison.
analysis_settings <- tibble(
  MCD_package = "robustbase::covMcd",
  MCD_alpha = mcd_alpha,
  MCD_cutoff_probability = mcd_cutoff_probability,
  MCD_cutoff_distribution = "scaled F",
  MCD_cutoff_method = "GM14 (Green-Martin extension)",
  MCD_cutoff_variant = mcd_cutoff_variant,
  MCD_before_input = "Original ICA matrix M; no robust-z standardization",
  MCD_after_input = "SHASH-i transformed ICA matrix",
  Univariate_threshold = univariate_threshold,
  MCD_plot_minimum_cap = mcd_plot_minimum_cap,
  MCD_plot_cutoff_multiplier = mcd_plot_cutoff_multiplier
)


# -----------------------------------------------------------------------------
# 14. Save combined data objects as RDS
# -----------------------------------------------------------------------------

saveRDS(
  dimension_summary,
  file.path(root_out, "ABIDE_direct_ICA_dimension_summary.rds")
)

saveRDS(
  timing_summary,
  file.path(root_out, "ABIDE_direct_ICA_timing_summary.rds")
)

saveRDS(
  IC_outlier_summary,
  file.path(root_out, "ABIDE_direct_ICA_outlier_summary_before_after.rds")
)

saveRDS(
  before_outlier_details,
  file.path(root_out, "ABIDE_direct_ICA_before_transformation_outliers_all_ICs.rds")
)

saveRDS(
  after_outlier_details,
  file.path(root_out, "ABIDE_direct_ICA_after_transformation_outliers_all_ICs.rds")
)

saveRDS(
  before_after_comparison,
  file.path(root_out, "ABIDE_direct_ICA_before_after_IC_volume_comparison.rds")
)

saveRDS(
  visual_inspection_volume_summary,
  file.path(root_out, "ABIDE_direct_ICA_visual_inspection_volume_summary.rds")
)

saveRDS(
  final_univariate_union_summary,
  file.path(root_out, "ABIDE_direct_ICA_final_SHASHi_union_outlier_summary.rds")
)

saveRDS(
  outlier_consistency_check,
  file.path(root_out, "ABIDE_direct_ICA_SHASH_outlier_consistency_check.rds")
)

saveRDS(
  QQ_outlier_labels_all_ICs,
  file.path(root_out, "ABIDE_direct_ICA_QQ_outlier_labels_before_after_all_ICs.rds")
)

saveRDS(
  QQ_data_all,
  file.path(root_out, "ABIDE_direct_ICA_all_QQ_plot_data.rds")
)

saveRDS(
  values_capped_for_plotting,
  file.path(root_out, "ABIDE_direct_ICA_QQ_values_capped_only_for_plotting.rds")
)

saveRDS(
  MCD_distances_all,
  file.path(root_out, "ABIDE_direct_ICA_MCD_distances_all.rds")
)

saveRDS(
  MCD_outlier_list,
  file.path(root_out, "ABIDE_direct_ICA_MCD_outlier_list.rds")
)

saveRDS(
  MCD_values_capped_for_plotting,
  file.path(root_out, "ABIDE_direct_ICA_MCD_values_capped_only_for_plotting.rds")
)

saveRDS(
  analysis_settings,
  file.path(root_out, "ABIDE_direct_ICA_analysis_settings.rds")
)

saveRDS(
  MCD_volume_summary_all,
  file.path(root_out, "ABIDE_direct_ICA_MCD_volume_summary_all.rds")
)

saveRDS(
  combined_detection_summary,
  file.path(root_out, "ABIDE_direct_ICA_combined_detection_summary.rds")
)


# -----------------------------------------------------------------------------
# 15. Save complete R objects
# -----------------------------------------------------------------------------

saveRDS(
  analysis_results,
  file = file.path(
    root_out,
    "ABIDE_direct_ICA_SHASHi_MCD075_GM14F_q9973_analysis_results.rds"
  )
)

complete_analysis_bundle <- list(
  analysis_results = analysis_results,
  dimension_summary = dimension_summary,
  timing_summary = timing_summary,
  IC_outlier_summary = IC_outlier_summary,
  before_outlier_details = before_outlier_details,
  after_outlier_details = after_outlier_details,
  before_after_comparison = before_after_comparison,
  visual_inspection_volume_summary = visual_inspection_volume_summary,
  final_univariate_union_summary = final_univariate_union_summary,
  outlier_consistency_check = outlier_consistency_check,
  QQ_outlier_labels_all_ICs = QQ_outlier_labels_all_ICs,
  QQ_data_all = QQ_data_all,
  values_capped_for_plotting = values_capped_for_plotting,
  MCD_distances_all = MCD_distances_all,
  MCD_outlier_list = MCD_outlier_list,
  MCD_values_capped_for_plotting = MCD_values_capped_for_plotting,
  analysis_settings = analysis_settings,
  MCD_volume_summary_all = MCD_volume_summary_all,
  combined_detection_summary = combined_detection_summary,
  all_IC_plots = all_IC_plots,
  first6_IC_plots = first6_IC_plots,
  QQ_first6_all = QQ_first6_all,
  MCD_distance_plots = MCD_distance_plots
)

saveRDS(
  complete_analysis_bundle,
  file = file.path(
    root_out,
    "ABIDE_direct_ICA_SHASHi_MCD075_GM14F_q9973_complete_bundle.rds"
  )
)


# -----------------------------------------------------------------------------
# 16. Record package and session information
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    root_out,
    "sessionInfo.txt"
  )
)


# -----------------------------------------------------------------------------
# 17. Print principal outputs
# -----------------------------------------------------------------------------

cat(
  "\n========================================\n",
  "Dimension summary\n",
  "========================================\n",
  sep = ""
)

print(dimension_summary, n = Inf, width = Inf)

cat(
  "\n========================================\n",
  "Direct ICA SHASH-i union of detected volumes\n",
  "========================================\n",
  sep = ""
)

print(final_univariate_union_summary, n = Inf, width = Inf)

cat(
  "\n========================================\n",
  "MCD outlier list\n",
  "========================================\n",
  sep = ""
)

print(MCD_outlier_list, n = Inf, width = Inf)

cat(
  "\n========================================\n",
  "MCD distances capped only for plotting\n",
  "========================================\n",
  sep = ""
)

print(MCD_values_capped_for_plotting, n = Inf, width = Inf)

cat(
  "\n========================================\n",
  "Combined direct-ICA detection summary\n",
  "========================================\n",
  sep = ""
)

print(
  combined_detection_summary %>% filter(Final_status != "not flagged"),
  n = Inf,
  width = Inf
)

cat(
  "\n========================================\n",
  "QQ-plot outlier labels: before and after\n",
  "========================================\n",
  sep = ""
)

print(QQ_outlier_labels_all_ICs, n = Inf, width = Inf)

cat(
  "\n========================================\n",
  "Values capped only for plotting\n",
  "========================================\n",
  sep = ""
)

print(values_capped_for_plotting, n = Inf, width = Inf)

cat(
  "\n========================================\n",
  "Analysis complete.\n",
  "Outputs saved in:\n",
  root_out,
  "\n========================================\n",
  sep = ""
)

# =============================================================================
# 18. Comparison-ready summaries
# =============================================================================

# Unique volumes flagged before/after SHASH-i by at least one IC.
univariate_volume_count_summary <- visual_inspection_volume_summary %>%
  group_by(Dataset, Subject) %>%
  summarise(
    before_SHASHi = sum(Flagged_before_any, na.rm = TRUE),
    after_SHASHi = sum(Flagged_after_any, na.rm = TRUE),
    flagged_both = sum(Flagged_before_any & Flagged_after_any, na.rm = TRUE),
    before_only = sum(Flagged_before_any & !Flagged_after_any, na.rm = TRUE),
    after_only = sum(!Flagged_before_any & Flagged_after_any, na.rm = TRUE),
    .groups = "drop"
  )

print(univariate_volume_count_summary, n = Inf, width = Inf)

volume_number_list <- visual_inspection_volume_summary %>%
  group_by(Dataset, Subject) %>%
  summarise(
    Volumes_before_SHASHi = paste(Volume[Flagged_before_any], collapse = ", "),
    Volumes_after_SHASHi = paste(Volume[Flagged_after_any], collapse = ", "),
    Volumes_flagged_both = paste(Volume[Flagged_before_any & Flagged_after_any], collapse = ", "),
    Volumes_before_only = paste(Volume[Flagged_before_any & !Flagged_after_any], collapse = ", "),
    Volumes_after_only = paste(Volume[!Flagged_before_any & Flagged_after_any], collapse = ", "),
    .groups = "drop"
  )

print(volume_number_list, n = Inf, width = Inf)

saveRDS(
  univariate_volume_count_summary,
  file.path(root_out, "ABIDE_univariate_volume_count_summary.rds")
)

saveRDS(
  volume_number_list,
  file.path(root_out, "ABIDE_volume_numbers_before_after_SHASHi.rds")
)

# ICA version of number-of-components flagging each volume.
ica_plot_data <- visual_inspection_volume_summary %>%
  select(
    Dataset, Subject, Volume,
    Number_of_ICs_before, Number_of_ICs_after
  ) %>%
  tidyr::pivot_longer(
    cols = c(Number_of_ICs_before, Number_of_ICs_after),
    names_to = "Stage",
    values_to = "Number_of_ICs"
  ) %>%
  mutate(
    Stage = recode(
      Stage,
      "Number_of_ICs_before" = "Before SHASH-i",
      "Number_of_ICs_after" = "After SHASH-i"
    ),
    Stage = factor(Stage, levels = c("Before SHASH-i", "After SHASH-i")),
    Label = ifelse(Number_of_ICs >= 2, as.character(Volume), "")
  )

saveRDS(
  ica_plot_data,
  file.path(root_out, "ABIDE_direct_ICA_number_of_ICs_flagging_each_volume_data.rds")
)

p_ica_counts <- ggplot(
  ica_plot_data %>% filter(Number_of_ICs > 0),
  aes(x = Volume, y = Number_of_ICs, color = Stage)
) +
  geom_segment(aes(xend = Volume, yend = 0), linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = Label),
    vjust = -0.45,
    size = 3.3,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Before SHASH-i" = "steelblue",
      "After SHASH-i" = "tomato"
    )
  ) +
  facet_grid(Stage ~ Subject, scales = "free_x") +
  labs(
    title = "Number of ICs flagging each volume",
    subtitle = "Volumes detected in at least two ICs are labelled",
    x = "Volume number",
    y = "Number of ICs",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11)
  )

print(p_ica_counts)

ggsave(
  filename = file.path(root_out, "ABIDE_direct_ICA_number_of_ICs_flagging_each_volume.png"),
  plot = p_ica_counts,
  width = 12, height = 7, dpi = 300
)

ggsave(
  filename = file.path(root_out, "ABIDE_direct_ICA_number_of_ICs_flagging_each_volume.pdf"),
  plot = p_ica_counts,
  width = 12, height = 7
)

# MCD-based detected-volume lists before and after SHASH-i.
mcd_volume_number_list <- MCD_distances_all %>%
  group_by(Dataset, Subject) %>%
  summarise(
    MCD_before_SHASHi = paste(
      Volume[Stage == "MCD on original ICA before SHASH-i" & MCD_outlier],
      collapse = ", "
    ),
    MCD_after_SHASHi = paste(
      Volume[Stage == "MCD on SHASH-i transformed ICA" & MCD_outlier],
      collapse = ", "
    ),
    MCD_flagged_both = paste(
      intersect(
        Volume[Stage == "MCD on original ICA before SHASH-i" & MCD_outlier],
        Volume[Stage == "MCD on SHASH-i transformed ICA" & MCD_outlier]
      ),
      collapse = ", "
    ),
    MCD_before_only = paste(
      setdiff(
        Volume[Stage == "MCD on original ICA before SHASH-i" & MCD_outlier],
        Volume[Stage == "MCD on SHASH-i transformed ICA" & MCD_outlier]
      ),
      collapse = ", "
    ),
    MCD_after_only = paste(
      setdiff(
        Volume[Stage == "MCD on SHASH-i transformed ICA" & MCD_outlier],
        Volume[Stage == "MCD on original ICA before SHASH-i" & MCD_outlier]
      ),
      collapse = ", "
    ),
    .groups = "drop"
  )

print(mcd_volume_number_list, n = Inf, width = Inf)

saveRDS(
  mcd_volume_number_list,
  file = file.path(
    root_out,
    "ABIDE_direct_ICA_MCD_volume_numbers_before_after_SHASHi.rds"
  )
)

## to plot the zoomed in line plots from the selected time segments
# for ABIDE 1
p_zoom_ABIDE1 <- make_MCD_distance_plot(
  mcd_data = MCD_distances_all %>%
    filter( Dataset == "ABIDE1", Volume >= 40, Volume <= 160),
  subject_label = "Subject 1 (high noise)"
) + scale_x_continuous(breaks = 40:160, minor_breaks = NULL
) 

print(p_zoom_ABIDE1)
ggsave(file.path(root_out, "ABIDE1", "ABIDE1_MCD_zoom_40_160.png"),
       p_zoom_ABIDE1, width = 10, height = 5.5, dpi = 300, bg = "white")



# for ABIDE 2
p_zoom_ABIDE2 <- make_MCD_distance_plot(
  mcd_data = MCD_distances_all %>%
    filter(Dataset == "ABIDE2", Volume >= 99, Volume <= 110),
  subject_label = "Subject 2 (low noise)"
) +
  scale_x_continuous(
    breaks = 99:110,
    minor_breaks = NULL
  ) +
  geom_vline(
    xintercept = 99:110,
    color = "gray85",
    linewidth = 0.35
  )

print(p_zoom_ABIDE2)

ggsave(
  file.path(root_out, "ABIDE2", "ABIDE2_MCD_zoom_99_110.png"),
  p_zoom_ABIDE2,
  width = 10,
  height = 5.5,
  dpi = 300,
  bg = "white"
)


### loading the saved data to see the list of identified volumes:
library(dplyr)
library(tidyr)

# Optional: reload this RDS object in a fresh R session if needed.
# mcd_volume_number_list <- readRDS(
#   file.path(
#     root_out,
#     "ABIDE_direct_ICA_MCD_volume_numbers_before_after_SHASHi.rds"
#   )
# )

mcd_long <- bind_rows(
  analysis_results$ABIDE1$MCD_distances,
  analysis_results$ABIDE2$MCD_distances
)

View(mcd_long)

mcd_outliers_long <- mcd_long %>%
  filter(MCD_outlier) %>%
  select(
    Dataset,
    Subject,
    Stage,
    Volume,
    MCD_distance,
    MCD_cutoff
  ) %>%
  arrange(
    Dataset,
    Stage,
    Volume
  )

View(mcd_outliers_long)

mcd_outliers_long <- mcd_outliers_long %>%
  mutate(
    Stage = case_when(
      Stage == "MCD on original ICA before SHASH-i" ~ "Before SHASH-i",
      Stage == "MCD on SHASH-i transformed ICA"    ~ "After SHASH-i",
      TRUE ~ Stage
    )
  )

# Before SHASH-i
abide1_before <- mcd_outliers_long %>%
  filter(
    Dataset == "ABIDE1",
    Stage == "Before SHASH-i"
  ) %>%
  pull(Volume) %>%
  sort() %>%
  unique()

abide1_before
# After SHASH-i
abide1_after <- mcd_outliers_long %>%
  filter(
    Dataset == "ABIDE1",
    Stage == "After SHASH-i"
  ) %>%
  pull(Volume) %>%
  sort() %>%
  unique()

abide1_after


#identified volumes to compare
both <- intersect(abide1_before, abide1_after)

before_only <- setdiff(abide1_before, abide1_after)

after_only <- setdiff(abide1_after, abide1_before)


# I wonder their MCD distance magnitudes:
library(dplyr)
library(tidyr)

abide1_mcd <- analysis_results$ABIDE1$MCD_distances

selected_volumes <- c(58, 59, 60, 61,
                      89, 90, 91, 92,
                      141, 142, 143, 144, 145)

abide1_mcd_selected <- abide1_mcd %>%
  filter(Volume %in% selected_volumes) %>%
  mutate(
    Method = case_when(
      grepl("before SHASH-i", Stage) ~ "Before",
      grepl("transformed ICA", Stage) ~ "After",
      TRUE ~ Stage
    )
  ) %>%
  select(Volume, Method, MCD_distance, MCD_cutoff, MCD_outlier) %>%
  pivot_wider(
    names_from = Method,
    values_from = c(MCD_distance, MCD_outlier),
    names_sep = "_"
  ) %>%
  mutate(
    Status = case_when(
      MCD_outlier_Before & MCD_outlier_After ~ "both",
      MCD_outlier_Before & !MCD_outlier_After ~ "before only",
      !MCD_outlier_Before & MCD_outlier_After ~ "after only",
      TRUE ~ "none"
    )
  ) %>%
  arrange(Volume)

abide1_mcd_selected

saveRDS(
  abide1_mcd_selected,
  file.path(root_out, "ABIDE1_selected_volumes_MCD_distances.rds")
)


library(dplyr)
library(ggplot2)

# -------------------------------------------------------------------------
# ABIDE1: MCD robust distances, volumes 40--160
# -------------------------------------------------------------------------

abide1_zoom_data <- MCD_distances_all %>%
  filter(
    Dataset == "ABIDE1",
    Volume >= 40,
    Volume <= 160
  )

# Same idea as the full plot: cap very large distances only for display
plot_cap <- max(
  50,
  5 * max(abide1_zoom_data$MCD_cutoff, na.rm = TRUE)
)

abide1_zoom_data <- abide1_zoom_data %>%
  mutate(
    Was_capped = MCD_distance > plot_cap,
    MCD_distance_display = pmin(MCD_distance, plot_cap)
  )

capped_label_data <- abide1_zoom_data %>%
  filter(Was_capped) %>%
  mutate(
    label_x = case_when(
      Stage == "MCD on original ICA before SHASH-i" ~ Volume - 0.8,
      Stage == "MCD on SHASH-i transformed ICA"     ~ Volume + 0.8
    ),
    label_y = MCD_distance_display - 0.3
  )

# Cutoff information
cutoff_data <- abide1_zoom_data %>%
  distinct(
    Stage,
    MCD_cutoff
  )


# -------------------------------------------------------------------------
# Plot directly with ggplot
# -------------------------------------------------------------------------
selected_volumes <- c(
  58:61,
  89:92,
  141:145
)

p_ABIDE1_40_160 <- ggplot(
  abide1_zoom_data,
  aes(
    x = Volume,
    y = MCD_distance_display,
    color = Stage
  )
) +
  geom_vline(
    xintercept = selected_volumes,
    color = "gray85",
    linewidth = 0.5
  ) +
  
  # Connect robust distances across volumes
  geom_line(
    linewidth = 0.6,
    alpha = 0.85
  ) +
  
  # Show detected outliers as points
  geom_point(
    data = abide1_zoom_data %>%
      filter(MCD_outlier),
    size = 2.2,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_text(
    data = capped_label_data,
    aes(
      x = label_x,
      y = label_y,
      label = "\u2191",
      color = Stage
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  # MCD cutoff
  geom_hline(
    data = cutoff_data,
    aes(yintercept = MCD_cutoff),
    inherit.aes = FALSE,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "MCD on original ICA before SHASH-i" = "steelblue",
      "MCD on SHASH-i transformed ICA" = "tomato"
    ),
    labels = c(
      "without normality transformation",
      "with normality transformation"
    )
  ) +
  scale_x_continuous(
    breaks = seq(40, 160, by = 10),
    minor_breaks = NULL
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  labs(
    title = "Subject 1 (high noise): MCD robust distances in ICA space",
    subtitle = paste0(
      "Volumes 40--160; ",
      "robustbase::covMcd(alpha = 0.75); ",
      "dashed line = GM14 scaled-F q0.9973 cutoff"
    ),
    x = "Volume number",
    y = "MCD robust distance",
    color = NULL
  ) +
  
  theme_minimal(base_size = 15) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(size = 14)
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1.2
      )
    )
  )

print(p_ABIDE1_40_160)
ggsave(
  filename = "ABIDE1_MCD_vol40_160.png",
  plot = p_ABIDE1_40_160,
  path = file.path(root_out, "ABIDE1"),
  width = 10,
  height = 5,
  dpi = 300
)

