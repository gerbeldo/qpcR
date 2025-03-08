#' @importFrom dplyr mutate group_by ungroup row_number filter summarise left_join
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar geom_point position_jitter
#' @importFrom stringr str_remove str_starts str_detect fixed
#' @importFrom vroom vroom
#' @importFrom janitor clean_names
#' @importFrom lme4 lmer
#' @import lmerTest
NULL

#' Add Sample Metadata
#'
#' Adds metadata to a wells data frame by extracting the sample type, biological replicate,
#' and calculating technical replicate numbers.
#'
#' @param wells A data frame containing wells information.
#' @param bio_rep_re A regular expression pattern used to remove the biological replicate number from the sample name. Defaults to `"\\s\\d+$"`.
#'
#' @return A modified data frame with additional columns: `sample_type`, `bio_rep`, and `tech_rep`.
add_sample_metadata <- function(wells, bio_rep_re = "\\s\\d+$") {
  wells |>
    mutate(
      sample_type = stringr::str_remove(sample, bio_rep_re),
      bio_rep = as.integer(stringr::str_remove(sample, stringr::fixed(sample_type)))
    ) |>
    group_by(sample, target) |>
    mutate(tech_rep = row_number()) |>
    ungroup()
}

#' Calculate Average Housekeeping Cq
#'
#' Filters the wells for a specific housekeeping target and computes the average Cq value for each sample.
#'
#' @param wells A data frame containing wells information.
#' @param housekeeping_target A character string specifying the housekeeping gene target.
#'
#' @return A data frame with columns `sample` and `mean_cq`, representing the average Cq value per sample.
average_housekeeping <- function(wells, housekeeping_target) {
  wells |>
    filter(target == housekeeping_target) |>
    group_by(sample) |>
    summarise(mean_cq = mean(cq))
}

#' Calculate Delta Cq (DCq)
#'
#' Computes the Delta Cq (DCq) values by subtracting the mean housekeeping Cq from each sample's Cq.
#'
#' @param wells A data frame containing wells information.
#' @param housekeeping_target A character string specifying the housekeeping gene target.
#'
#' @return A data frame with an added column `DCq` representing the difference between the Cq and the mean housekeeping Cq.
calculate_DCq <- function(wells, housekeeping_target) {
  mean_hk_per_sample <- average_housekeeping(wells, housekeeping_target)

  wells |>
    left_join(mean_hk_per_sample, by = "sample") |>
    mutate(DCq = cq - mean_cq)
}

#' Calculate Relative Expression
#'
#' Computes the relative expression values using the formula 2^(-DCq).
#'
#' @param wells A data frame containing wells information with a `DCq` column.
#'
#' @return A data frame with an additional column `rel_exp` representing the relative expression.
calculate_relative_expression <- function(wells) {
  wells |> mutate(rel_exp = 2^(-DCq))
}

#' Normalize Relative Expression to Mean Wild-Type
#'
#' Normalizes the relative expression values by dividing by the mean relative expression of wild-type samples.
#'
#' @param wells A data frame containing wells information with a `rel_exp` column.
#' @param target A character string specifying the target gene to use for normalization.
#' @param WT_string A character string used to identify wild-type samples. Defaults to `"WT"`.
#'
#' @return A data frame with an additional column `rel_exp_norm` representing the normalized relative expression.
normalize_to_mean_WT <- function(wells, target, WT_string = "WT") {
  WT_norm_factor <- wells |>
    filter(stringr::str_starts(sample, WT_string), target == target) |>
    summarise(mean_rel_exp = mean(rel_exp)) |>
    pull(mean_rel_exp)

  wells |> mutate(rel_exp_norm = rel_exp / WT_norm_factor)
}

#' Create a Superplot
#'
#' Generates a superplot that displays summary statistics (mean and standard error) as well as individual data points
#' for a specified target.
#'
#' @param wells A data frame containing wells information.
#' @param target A character string specifying the target gene to plot.
#' @param plot_order An optional character vector specifying the order of sample types. If provided, sample types are
#' converted to a factor with the specified levels.
#'
#' @return A ggplot object representing the superplot.
make_superplot <- function(wells, target, plot_order) {
  if (!is.null(plot_order)) {
    wells$sample_type <- factor(wells$sample_type, levels = plot_order)
  }
  stats <- wells |>
    filter(target == target) |>
    group_by(sample_type) |>
    summarise(
      mean_rel_exp_norm = mean(rel_exp_norm),
      sd_rel_exp_norm = sd(rel_exp_norm),
      se_rel_exp_norm = sd_rel_exp_norm / sqrt(n())
    )

  ggplot(stats, aes(
    x = sample_type, y = mean_rel_exp_norm,
    ymin = mean_rel_exp_norm - se_rel_exp_norm,
    ymax = mean_rel_exp_norm + se_rel_exp_norm,
    fill = sample_type
  )) +
    geom_col() +
    geom_errorbar(width = 0.2) +
    geom_point(
      data = wells,
      aes(x = sample_type, y = rel_exp_norm, shape = as.factor(bio_rep)),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.2)
    )
}

#' Perform Statistical Analysis
#'
#' Applies a linear mixed effects model (using lmer) to the normalized relative expression data,
#' analyzing the effect of sample type with biological replicates as random effects.
#'
#' @param wells A data frame containing wells information with `rel_exp_norm` computed.
#' @param sample_type_levels An optional character vector specifying the factor levels for `sample_type`. If `NULL`,
#' the wild-type level (matched by `WT_string`) is placed first.
#' @param WT_string A character string used to identify wild-type samples. Defaults to `"WT"`.
#'
#' @return A summary of the fitted linear mixed effects model.
#'
#' @export
perform_statistics <- function(wells, sample_type_levels = NULL, WT_string = "WT") {
  if (is.null(sample_type_levels)) {
    unique_samples <- unique(wells$sample_type)
    WT_level <- unique_samples[stringr::str_detect(unique_samples, WT_string)][1]
    sample_type_levels <- c(WT_level, setdiff(unique_samples, WT_level))
  }

  wells$sample_type <- factor(wells$sample_type, levels = sample_type_levels)
  wells$bio_rep <- factor(wells$bio_rep)

  model <- lmer(rel_exp_norm ~ sample_type + (1 | bio_rep), data = wells)

  return(summary(model))
}

#' Preprocess qPCR Data
#'
#' Reads and preprocesses qPCR data from a file. The function cleans column names,
#' filters out omitted data, calculates Delta Cq (DCq), computes relative expression,
#' and normalizes the expression to the wild-type (WT) samples.
#'
#' @param file_path A character string specifying the path to the input data file.
#' @param housekeeping_target A character string specifying the housekeeping gene target. Defaults to `"GAPDH"`.
#' @param WT_string A character string used to identify wild-type samples. Defaults to `"WT"`.
#'
#' @return A processed data frame with computed columns including `DCq`, `rel_exp`, and `rel_exp_norm`.
#'
#' @export
preprocess_data <- function(file_path, housekeeping_target = "GAPDH", WT_string = "WT") {
  wells <- vroom::vroom(file_path, na = c("N/A", "Undetermined")) |>
    janitor::clean_names() |>
    add_sample_metadata()

  wells <- wells |> filter(omit == FALSE)

  wells <- calculate_DCq(wells, housekeeping_target)

  clean_wells <- wells |> filter(target != housekeeping_target)

  clean_wells <- calculate_relative_expression(clean_wells)

  clean_wells <- normalize_to_mean_WT(clean_wells, housekeeping_target, WT_string)

  return(clean_wells)
}

#' Plot qPCR Data
#'
#' Creates a superplot for a specified target using preprocessed data.
#'
#' @param wells A processed data frame containing wells information.
#' @param target A character string specifying the target gene to plot.
#' @param plot_order An optional character vector specifying the order of sample types.
#'
#' @return A ggplot object representing the superplot.
#'
#' @export
plot_data <- function(wells, target, plot_order = NULL) {
  make_superplot(wells, target, plot_order)
}
