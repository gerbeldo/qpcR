library(dplyr)
library(ggplot2)
library(vroom)
library(lme4)
library(lmerTest)

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

average_housekeeping <- function(wells, housekeeping_target) {
  wells |>
    filter(target == housekeeping_target) |>
    group_by(sample) |>
    summarise(mean_cq = mean(cq))
}

calculate_DCq <- function(wells, housekeeping_target) {
  mean_hk_per_sample <- average_housekeeping(wells, housekeeping_target)

  wells |>
    left_join(mean_hk_per_sample, by = "sample") |>
    mutate(DCq = cq - mean_cq)
}

calculate_relative_expression <- function(wells) {
  wells |> mutate(rel_exp = 2^(-DCq))
}

normalize_to_mean_WT <- function(wells, target, WT_string = "WT") {
  WT_norm_factor <- wells |>
    filter(stringr::str_starts(sample, WT_string), target == target) |>
    summarise(mean_rel_exp = mean(rel_exp)) |>
    pull(mean_rel_exp)

  wells |> mutate(rel_exp_norm = rel_exp / WT_norm_factor)
}

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

plot_data <- function(wells, target, plot_order = NULL) {
  make_superplot(wells, target, plot_order)
}
