#### Preamble ####
# Purpose: Conducts Kitagawa decomposition of Black-White infant mortality gap
# Author: Monica Alexander
# Date: 2026-01-31
# Contact: monica.alexander@utoronto.ca
# License: MIT
# Pre-requisites:
#   - Run 01-clean_data.R first to create data/analysis_data/analysis_data.csv
#   - Required R packages: tidyverse
# Outputs:
#   - data/analysis_data/decomposition_results.csv (annual decomposition)
#   - data/analysis_data/detailed_decomposition.csv (by gestational age)
#   - data/analysis_data/decomposition_by_mother_age.csv (by mother's age)

#### Workspace setup ####
library(tidyverse)

set.seed(853)  # For reproducibility

# Create output directory
dir.create("data/analysis_data", recursive = TRUE, showWarnings = FALSE)

#### Load data ####

analysis_data <- read_csv(
  "data/analysis_data/analysis_data.csv",
  show_col_types = FALSE
) |>
  filter(year >= 2003)

cat("Data loaded:", nrow(analysis_data), "observations\n")
cat("Years:", min(analysis_data$year), "-", max(analysis_data$year), "\n\n")

#### Kitagawa Decomposition ####
#
# Method: Kitagawa (1955) decomposition of rate differences
#
# The difference between two rates can be decomposed into:
#
# R_B - R_W = Composition Effect + Rate Effect
#
# Where:
# Composition Effect = Sum_i[(p_B_i - p_W_i) * (M_B_i + M_W_i)/2]
# Rate Effect = Sum_i[(M_B_i - M_W_i) * (p_B_i + p_W_i)/2]
#
# - R_B, R_W = overall infant mortality rates for Black and White
# - p_B_i, p_W_i = proportion of births in gestational age group i
# - M_B_i, M_W_i = infant mortality rate in gestational age group i

#### Monte Carlo CI helpers ####
#
# 95% confidence intervals around rates and Kitagawa decomposition results are
# obtained by treating observed death counts as Poisson with mean equal to the
# observed count, simulating death counts for each (race, gestational age,
# [cause]) cell, recomputing rates and the decomposition, and taking the 2.5%
# and 97.5% quantiles of the resulting distributions.
#
# Births are treated as fixed (large counts; sampling variability dominated by
# the death side).

n_sim_default <- 1000

# Simulate n_sim Kitagawa decompositions for one year. Returns a tibble of
# n_sim rows with imr_black, imr_white, total_difference, composition_effect,
# rate_effect, composition_pct, rate_pct.
simulate_kitagawa_year <- function(year_data, n_sim = n_sim_default) {
  black <- year_data |>
    filter(race == "Non-Hispanic Black") |>
    arrange(gestational_age_group)
  white <- year_data |>
    filter(race == "Non-Hispanic White") |>
    arrange(gestational_age_group)
  stopifnot(identical(as.character(black$gestational_age_group),
                      as.character(white$gestational_age_group)))

  p_b <- black$birth_proportion
  p_w <- white$birth_proportion
  d_b <- black$deaths
  d_w <- white$deaths
  b_b <- black$births
  b_w <- white$births

  # Simulated death counts: n_sim x n_ga matrices
  d_b_sim <- sapply(d_b, function(d) rpois(n_sim, d))
  d_w_sim <- sapply(d_w, function(d) rpois(n_sim, d))
  if (is.null(dim(d_b_sim))) d_b_sim <- matrix(d_b_sim, nrow = n_sim)
  if (is.null(dim(d_w_sim))) d_w_sim <- matrix(d_w_sim, nrow = n_sim)

  m_b_sim <- sweep(d_b_sim, 2, b_b, "/") * 1000
  m_w_sim <- sweep(d_w_sim, 2, b_w, "/") * 1000

  R_b_sim <- as.numeric(m_b_sim %*% p_b)
  R_w_sim <- as.numeric(m_w_sim %*% p_w)
  total_sim <- R_b_sim - R_w_sim

  dp <- p_b - p_w
  sp <- p_b + p_w
  comp_sim <- as.numeric(((m_b_sim + m_w_sim) / 2) %*% dp)
  rate_sim <- as.numeric(((m_b_sim - m_w_sim) / 2) %*% sp)

  tibble(
    sim_id = seq_len(n_sim),
    imr_black = R_b_sim,
    imr_white = R_w_sim,
    total_difference = total_sim,
    composition_effect = comp_sim,
    rate_effect = rate_sim,
    composition_pct = if_else(total_sim == 0, NA_real_, comp_sim / total_sim * 100),
    rate_pct = if_else(total_sim == 0, NA_real_, rate_sim / total_sim * 100)
  )
}

# Simulate n_sim cause-specific Kitagawa decompositions for one (year, cause).
# year_cause_data must have: race, gestational_age_group, p_b, p_w,
# deaths, births (cause-specific deaths, total births in cell).
simulate_kitagawa_cause_year <- function(black, white, n_sim = n_sim_default) {
  # black/white are tibbles with rows in matching GA order containing
  # birth_proportion (p), deaths (cause-specific), births (cell totals)
  stopifnot(identical(as.character(black$gestational_age_group),
                      as.character(white$gestational_age_group)))
  p_b <- black$birth_proportion
  p_w <- white$birth_proportion
  d_b <- black$deaths
  d_w <- white$deaths
  b_b <- black$births
  b_w <- white$births

  d_b_sim <- sapply(d_b, function(d) rpois(n_sim, d))
  d_w_sim <- sapply(d_w, function(d) rpois(n_sim, d))
  if (is.null(dim(d_b_sim))) d_b_sim <- matrix(d_b_sim, nrow = n_sim)
  if (is.null(dim(d_w_sim))) d_w_sim <- matrix(d_w_sim, nrow = n_sim)

  m_b_sim <- sweep(d_b_sim, 2, b_b, "/") * 1000
  m_w_sim <- sweep(d_w_sim, 2, b_w, "/") * 1000

  dp <- p_b - p_w
  sp <- p_b + p_w
  comp_sim <- as.numeric(((m_b_sim + m_w_sim) / 2) %*% dp)
  rate_sim <- as.numeric(((m_b_sim - m_w_sim) / 2) %*% sp)
  gap_sim <- comp_sim + rate_sim

  tibble(
    sim_id = seq_len(n_sim),
    composition_effect = comp_sim,
    rate_effect = rate_sim,
    gap = gap_sim
  )
}

# 95% CI bound helpers (2.5% / 97.5% quantiles), used by every CI summary
# block so the na.rm behaviour and quantile type stay consistent.
q_lo <- function(x) quantile(x, 0.025, names = FALSE, na.rm = TRUE)
q_hi <- function(x) quantile(x, 0.975, names = FALSE, na.rm = TRUE)

# Simulate n_sim per-gestational-age composition and rate contributions for
# each year by Poisson-resampling cell death counts. Returns a long tibble
# (year x gestational_age_group x sim_id) used to derive per-GA CIs under any
# gestational-age scheme. Shared by the main analysis and the 7-category
# sensitivity analysis.
simulate_detailed_contributions <- function(data, n_sim = n_sim_default) {
  years_d <- sort(unique(data$year))
  map_dfr(years_d, function(yr) {
    year_data <- data |> filter(year == yr)
    black <- year_data |> filter(race == "Non-Hispanic Black") |>
      arrange(gestational_age_group)
    white <- year_data |> filter(race == "Non-Hispanic White") |>
      arrange(gestational_age_group)
    # Guard against Black/White GA-group mismatch, which would silently
    # misalign the sweep() columns below.
    stopifnot(identical(as.character(black$gestational_age_group),
                        as.character(white$gestational_age_group)))
    p_b <- black$birth_proportion; p_w <- white$birth_proportion
    d_b <- black$deaths; d_w <- white$deaths
    b_b <- black$births; b_w <- white$births
    d_b_sim <- sapply(d_b, function(d) rpois(n_sim, d))
    d_w_sim <- sapply(d_w, function(d) rpois(n_sim, d))
    if (is.null(dim(d_b_sim))) d_b_sim <- matrix(d_b_sim, nrow = n_sim)
    if (is.null(dim(d_w_sim))) d_w_sim <- matrix(d_w_sim, nrow = n_sim)
    m_b_sim <- sweep(d_b_sim, 2, b_b, "/") * 1000
    m_w_sim <- sweep(d_w_sim, 2, b_w, "/") * 1000
    # Per-GA composition and rate contributions (n_sim x n_ga matrices)
    comp_mat <- sweep((m_b_sim + m_w_sim) / 2, 2, p_b - p_w, "*")
    rate_mat <- sweep((m_b_sim - m_w_sim) / 2, 2, p_b + p_w, "*")
    tibble(
      year = yr,
      gestational_age_group = rep(as.character(black$gestational_age_group),
                                  each = n_sim),
      sim_id = rep(seq_len(n_sim), times = ncol(comp_mat)),
      composition_contribution = as.numeric(comp_mat),
      rate_contribution = as.numeric(rate_mat)
    )
  })
}

# Reduce a long per-GA simulation tibble to per-(year, GA) 95% CI bounds for
# the composition, rate, and total contributions.
summarize_detailed_ci <- function(detailed_sims) {
  detailed_sims |>
    summarise(
      comp_lo = q_lo(composition_contribution),
      comp_hi = q_hi(composition_contribution),
      rate_lo = q_lo(rate_contribution),
      rate_hi = q_hi(rate_contribution),
      total_lo = q_lo(composition_contribution + rate_contribution),
      total_hi = q_hi(composition_contribution + rate_contribution),
      .by = c(year, gestational_age_group)
    )
}

kitagawa_decomposition <- function(data, year_val,
                                   prop_col = "birth_proportion",
                                   rate_col = "mortality_rate") {
  # Filter to specific year and join Black/White by GA group
  year_data <- data |>
    filter(year == year_val)

  joined <- year_data |>
    filter(race == "Non-Hispanic Black") |>
    select(gestational_age_group,
           p_b = all_of(prop_col), m_b = all_of(rate_col)) |>
    inner_join(
      year_data |>
        filter(race == "Non-Hispanic White") |>
        select(gestational_age_group,
               p_w = all_of(prop_col), m_w = all_of(rate_col)),
      by = "gestational_age_group"
    )

  # Overall infant mortality rates (weighted by birth distribution)
  R_b <- sum(joined$p_b * joined$m_b)
  R_w <- sum(joined$p_w * joined$m_w)
  total_diff <- R_b - R_w

  # Kitagawa decomposition
  composition_effect <- sum((joined$p_b - joined$p_w) * (joined$m_b + joined$m_w) / 2)
  rate_effect <- sum((joined$m_b - joined$m_w) * (joined$p_b + joined$p_w) / 2)

  tibble(
    year = year_val,
    imr_black = R_b,
    imr_white = R_w,
    total_difference = total_diff,
    composition_effect = composition_effect,
    rate_effect = rate_effect,
    composition_pct = composition_effect / total_diff * 100,
    rate_pct = rate_effect / total_diff * 100
  )
}

#### Run decomposition for all years ####

years <- sort(unique(analysis_data$year))
decomposition_results <- map_dfr(years, ~kitagawa_decomposition(analysis_data, .x))

#### Detailed decomposition by gestational age group ####

detailed_decomposition <- function(data, year_val) {
  year_data <- data |> filter(year == year_val)

  joined <- year_data |>
    filter(race == "Non-Hispanic Black") |>
    select(gestational_age_group,
           prop_black = birth_proportion, mortality_black = mortality_rate) |>
    inner_join(
      year_data |>
        filter(race == "Non-Hispanic White") |>
        select(gestational_age_group,
               prop_white = birth_proportion, mortality_white = mortality_rate),
      by = "gestational_age_group"
    )

  joined |>
    mutate(
      year = year_val,
      composition_contribution = (prop_black - prop_white) * (mortality_black + mortality_white) / 2,
      rate_contribution = (mortality_black - mortality_white) * (prop_black + prop_white) / 2,
      total_contribution = composition_contribution + rate_contribution
    ) |>
    select(year, gestational_age_group, prop_black, prop_white,
           mortality_black, mortality_white,
           composition_contribution, rate_contribution, total_contribution)
}

detailed_results <- map_dfr(years, ~detailed_decomposition(analysis_data, .x))

#### Display results ####

cat("=== Kitagawa Decomposition Results ===\n\n")

cat("First and last years:\n")
decomposition_results |>
  filter(year %in% c(min(year), max(year))) |>
  select(year, imr_black, imr_white, total_difference, composition_pct, rate_pct) |>
  mutate(across(where(is.numeric), ~round(.x, 2))) |>
  print()

cat("\nAverage contribution over study period:\n")
decomposition_results |>
  summarise(
    mean_composition_pct = mean(composition_pct),
    mean_rate_pct = mean(rate_pct)
  ) |>
  mutate(across(everything(), ~round(.x, 1))) |>
  print()

#### Create figures ####

# Color palette - coral and teal for better visual distinction
colors <- c("Gestational age distribution" = "#E07A5F", "Mortality risk" = "#3D7068")

# Figure 1: Stacked area chart of absolute contributions
# Area charts work better for longer time series (23 years)
# Rate effect on bottom, composition on top
fig1 <- decomposition_results |>
  select(year, `Gestational age distribution` = composition_effect,
         `Mortality risk` = rate_effect) |>
  pivot_longer(-year, names_to = "Component", values_to = "Contribution") |>
  mutate(Component = factor(Component,
                            levels = c("Mortality risk",
                                       "Gestational age distribution"))) |>
  ggplot(aes(x = year, y = Contribution, fill = Component)) +
  geom_area() +
  geom_line(
    data = decomposition_results |>
      mutate(total = composition_effect + rate_effect),
    aes(x = year, y = total, fill = NULL),
    linewidth = 0.5, color = "black"
  ) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = seq(2005, 2025, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Year",
    y = "Contribution to gap\n(deaths per 1,000 live births)",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Figure 2: Horizontal bar chart showing contribution by gestational age group
# Uses most recent year to show the breakdown of where the gap comes from
# Order gestational age from earliest to latest
ga_order <- c("<28 weeks", "28-31 weeks", "32-33 weeks", "34-36 weeks",
              "37-38 weeks", "39-41 weeks", "42+ weeks")

fig2 <- detailed_results |>
  filter(year == max(years)) |>
  select(gestational_age_group,
         `Gestational age distribution` = composition_contribution,
         `Mortality risk` = rate_contribution) |>
  pivot_longer(-gestational_age_group, names_to = "Component",
               values_to = "Contribution") |>
  mutate(
    gestational_age_group = factor(gestational_age_group, levels = rev(ga_order)),
    Component = factor(Component, levels = c("Mortality risk",
                                              "Gestational age distribution"))
  ) |>
  ggplot(aes(x = Contribution, y = gestational_age_group, fill = Component)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Contribution to gap (deaths per 1,000 live births)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Figure 3: Bar chart by gestational age for first and last years
fig3 <- detailed_results |>
  filter(year %in% c(min(years), max(years))) |>
  select(year, gestational_age_group, Composition = composition_contribution,
         Rate = rate_contribution) |>
  pivot_longer(c(Composition, Rate), names_to = "Component", values_to = "Contribution") |>
  mutate(
    Component = recode(Component,
                       "Composition" = "Gestational age distribution",
                       "Rate" = "Mortality risk"),
    year = factor(year)
  ) |>
  ggplot(aes(x = gestational_age_group, y = Contribution, fill = Component)) +
  geom_col(position = "dodge") +
  facet_wrap(~year, ncol = 1) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Gestational age group",
    y = "Contribution to gap\n(deaths per 1,000 live births)",
    fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

# Figure 4: Percent of disparity due to mortality risk over time
# Shows the declining contribution of the rate effect
fig4 <- decomposition_results |>
  ggplot(aes(x = year, y = rate_pct)) +
  geom_line(color = "#3D7068", linewidth = 1) +
  geom_point(color = "#3D7068", size = 2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = seq(2005, 2025, by = 5)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 25)) +
  labs(
    x = "Year",
    y = "Percent of disparity due to\nmortality risk differences"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

#### Stratification by Mother's Age ####

# Load mother's age data
analysis_data_by_age <- read_csv(
  "data/analysis_data/analysis_data_by_mother_age.csv",
  show_col_types = FALSE
) |>
  filter(year >= 2003)

# Order mother's age groups
mage_levels <- c("Under 20", "20-24", "25-29", "30-34", "35-39", "40 and over")

# Calculate decomposition by mother's age group for most recent year
decomp_by_mother_age <- function(data, year_val, age_group) {
  year_data <- data |>
    filter(year == year_val, mother_age_group == age_group)

  # Calculate birth proportions within each race
  year_data <- year_data |>
    mutate(
      birth_proportion = births / sum(births),
      mortality_rate = deaths / births * 1000,
      .by = race
    )

  # Join Black and White by gestational age group
  joined <- year_data |>
    filter(race == "Non-Hispanic Black") |>
    select(gestational_age_group,
           p_b = birth_proportion, m_b = mortality_rate) |>
    inner_join(
      year_data |>
        filter(race == "Non-Hispanic White") |>
        select(gestational_age_group,
               p_w = birth_proportion, m_w = mortality_rate),
      by = "gestational_age_group"
    )

  if (nrow(joined) == 0) return(NULL)

  R_b <- sum(joined$p_b * joined$m_b)
  R_w <- sum(joined$p_w * joined$m_w)
  total_diff <- R_b - R_w

  composition_effect <- sum((joined$p_b - joined$p_w) * (joined$m_b + joined$m_w) / 2)
  rate_effect <- sum((joined$m_b - joined$m_w) * (joined$p_b + joined$p_w) / 2)

  tibble(
    year = year_val,
    mother_age_group = age_group,
    imr_black = R_b,
    imr_white = R_w,
    total_difference = total_diff,
    composition_effect = composition_effect,
    rate_effect = rate_effect,
    composition_pct = if_else(total_diff == 0, NA_real_, composition_effect / total_diff * 100),
    rate_pct = if_else(total_diff == 0, NA_real_, rate_effect / total_diff * 100)
  )
}

# Calculate decomposition by mother's age for most recent year
last_year <- max(analysis_data_by_age$year)
decomp_by_age_results <- map_dfr(
  mage_levels,
  ~decomp_by_mother_age(analysis_data_by_age, last_year, .x)
) |>
  filter(!is.na(total_difference))

cat("\n=== Decomposition by Mother's Age (", last_year, ") ===\n")
decomp_by_age_results |>
  select(mother_age_group, imr_black, imr_white, total_difference, composition_pct, rate_pct) |>
  mutate(across(where(is.numeric), ~round(.x, 1))) |>
  print()

# Figure 5: Decomposition by mother's age group
fig5 <- decomp_by_age_results |>
  select(mother_age_group,
         `Gestational age distribution` = composition_effect,
         `Mortality risk` = rate_effect) |>
  pivot_longer(-mother_age_group, names_to = "Component", values_to = "Contribution") |>
  mutate(
    mother_age_group = factor(mother_age_group, levels = mage_levels),
    Component = factor(Component, levels = c("Mortality risk", "Gestational age distribution"))
  ) |>
  ggplot(aes(x = mother_age_group, y = Contribution, fill = Component)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Mother's age group",
    y = "Contribution to gap\n(deaths per 1,000 live births)",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# Figure 6: Percent composition effect by mother's age
fig6 <- decomp_by_age_results |>
  mutate(mother_age_group = factor(mother_age_group, levels = mage_levels)) |>
  ggplot(aes(x = mother_age_group, y = composition_pct)) +
  geom_col(fill = "#E07A5F") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  labs(
    x = "Mother's age group",
    y = "Percent of disparity due to\ngestational age distribution"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

#### Cause-of-Death Decomposition ####
# Kitagawa decomposition applied separately to each cause of death.
# For each cause j, the cause-specific gap decomposes into:
#   Composition_j = sum_i (p_B,i - p_W,i) * (M_B,i,j + M_W,i,j) / 2
#   Rate_j       = sum_i (M_B,i,j - M_W,i,j) * (p_B,i + p_W,i) / 2
# where p_r,i = birth proportion in GA group i, M_r,i,j = cause-specific
# mortality rate in GA group i.
#
# Raw ICD-10 codes are saved in cause_of_death_by_ga.csv by 01-clean_data.R.
# Categorization is done here so groupings can be adjusted without re-running
# the expensive data processing step.

cause_ga_raw <- read_csv("data/analysis_data/cause_of_death_by_ga.csv",
                         show_col_types = FALSE) |>
  filter(year >= 2003)

# Load shared cause-of-death categorization function
source("scripts/categorize_cause.R")

# Apply categorization and aggregate deaths by year/race/GA/cause
# Verify births is constant within each year-race-GA group before using first()
births_check <- cause_ga_raw |>
  distinct(year, race, gestational_age_group, births) |>
  count(year, race, gestational_age_group) |>
  filter(n > 1)
if (nrow(births_check) > 0) {
  stop("births is not constant within year-race-GA groups; cannot use first()")
}

cause_ga_data <- cause_ga_raw |>
  mutate(cause_of_death = categorize_cause(ucod)) |>
  summarise(deaths = sum(deaths), births = first(births),
            .by = c(year, race, gestational_age_group, cause_of_death))

# Get birth proportions from analysis_data (same source as main decomposition)
birth_props <- analysis_data |>
  select(year, race, gestational_age_group, birth_proportion)

# Compute cause-specific mortality rate per GA group: M_r,i,j = deaths/births * 1000
cause_ga_rates <- cause_ga_data |>
  mutate(rate = deaths / births * 1000)

# For each (year, cause), apply Kitagawa decomposition
# Need all combinations of year × GA group × cause × race (fill missing with 0)
all_combos <- expand_grid(
  year = unique(cause_ga_rates$year),
  race = c("Non-Hispanic Black", "Non-Hispanic White"),
  gestational_age_group = unique(birth_props$gestational_age_group),
  cause_of_death = unique(cause_ga_rates$cause_of_death)
)

cause_ga_full <- all_combos |>
  left_join(cause_ga_rates |> select(year, race, gestational_age_group, cause_of_death, rate),
            by = c("year", "race", "gestational_age_group", "cause_of_death")) |>
  mutate(rate = replace_na(rate, 0)) |>
  left_join(birth_props, by = c("year", "race", "gestational_age_group"))

# Pivot to wide by race for Kitagawa formula
cause_kitagawa <- cause_ga_full |>
  pivot_wider(
    names_from = race,
    values_from = c(rate, birth_proportion),
    names_sep = "_"
  ) |>
  rename(
    M_B = `rate_Non-Hispanic Black`,
    M_W = `rate_Non-Hispanic White`,
    p_B = `birth_proportion_Non-Hispanic Black`,
    p_W = `birth_proportion_Non-Hispanic White`
  ) |>
  mutate(
    comp_i = (p_B - p_W) * (M_B + M_W) / 2,
    rate_i = (M_B - M_W) * (p_B + p_W) / 2
  )

# Sum across GA groups for each (year, cause) to get cause-level effects
cause_decomp <- cause_kitagawa |>
  summarise(
    composition_effect = sum(comp_i),
    rate_effect = sum(rate_i),
    .by = c(year, cause_of_death)
  ) |>
  mutate(
    gap = composition_effect + rate_effect,
    composition_pct = if_else(gap == 0, NA_real_, composition_effect / gap * 100),
    rate_pct = if_else(gap == 0, NA_real_, rate_effect / gap * 100)
  )

# Also compute CSIMR for backwards compatibility
# Total births per (year, race) from distinct GA groups (births repeats per cause)
total_births_by_race <- cause_ga_data |>
  distinct(year, race, gestational_age_group, births) |>
  summarise(total_births = sum(births), .by = c(year, race))

cause_csimr <- cause_ga_data |>
  summarise(deaths = sum(deaths), .by = c(year, race, cause_of_death)) |>
  left_join(total_births_by_race, by = c("year", "race")) |>
  mutate(csimr = deaths / total_births * 1000) |>
  select(year, race, cause_of_death, csimr) |>
  pivot_wider(names_from = race, values_from = csimr, values_fill = 0) |>
  rename(csimr_black = `Non-Hispanic Black`, csimr_white = `Non-Hispanic White`)

cause_decomp <- cause_decomp |>
  left_join(cause_csimr, by = c("year", "cause_of_death"))

cat("\n=== Cause-of-Death Decomposition (Kitagawa per cause) ===\n")
cause_decomp |>
  filter(year == max(year)) |>
  select(cause_of_death, csimr_black, csimr_white, gap,
         composition_effect, rate_effect) |>
  mutate(across(where(is.numeric), ~round(.x, 3))) |>
  arrange(desc(gap)) |>
  print()

# Cause order for figures (by gap contribution in most recent year, descending)
cause_order <- cause_decomp |>
  filter(year == max(year)) |>
  arrange(desc(gap)) |>
  pull(cause_of_death)

# Top 8 causes for trend figure
top_causes <- head(cause_order, 8)

# Figure 7: Stacked bar chart showing composition vs rate per cause (most recent year)
cause_plot_data <- cause_decomp |>
  filter(year == max(year)) |>
  mutate(cause_of_death = factor(cause_of_death, levels = rev(cause_order)))

cause_plot_long <- cause_plot_data |>
  select(cause_of_death, gap, Composition = composition_effect,
         Rate = rate_effect) |>
  pivot_longer(cols = c(Rate, Composition), names_to = "Component",
               values_to = "value") |>
  mutate(Component = factor(Component, levels = c("Composition", "Rate")))

fig7 <- ggplot(cause_plot_long, aes(x = value, y = cause_of_death, fill = Component)) +
  geom_col(width = 0.7) +
  geom_text(data = cause_plot_data,
            aes(x = gap + 0.02, y = cause_of_death, label = round(gap, 2), fill = NULL),
            hjust = 0, size = 3, color = "gray30") +
  scale_fill_manual(values = c("Composition" = "#d4a373", "Rate" = "#588157")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Contribution to gap (deaths per 1,000 live births)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Figure 8: Single-panel line chart for top 8 causes (shared y-axis)
fig8 <- cause_decomp |>
  filter(cause_of_death %in% top_causes) |>
  mutate(cause_of_death = factor(cause_of_death, levels = top_causes)) |>
  ggplot(aes(x = year, y = gap, color = cause_of_death)) +
  geom_vline(xintercept = 2020, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  scale_x_continuous(breaks = seq(2005, 2025, by = 5)) +
  labs(x = NULL, y = "Deaths per 1,000 live births", color = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 3))

#### Monte Carlo 95% CIs ####

cat("\n=== Running Monte Carlo 95% CIs (n_sim =", n_sim_default, ") ===\n")

# Overall decomposition CIs ----
set.seed(853)
overall_sims <- map_dfr(years, function(yr) {
  simulate_kitagawa_year(analysis_data |> filter(year == yr)) |>
    mutate(year = yr, .before = 1)
})

overall_ci_long <- overall_sims |>
  summarise(
    across(c(imr_black, imr_white, total_difference,
             composition_effect, rate_effect, composition_pct, rate_pct),
           list(lo = q_lo, hi = q_hi)),
    .by = year
  )

decomposition_results_with_ci <- decomposition_results |>
  left_join(overall_ci_long, by = "year")

write_csv(decomposition_results_with_ci,
          "data/analysis_data/decomposition_results_with_ci.csv")
cat("Overall CIs saved: decomposition_results_with_ci.csv\n")

# IMR-only CIs (long format, for Fig 1A ribbons) ----
imr_with_ci <- bind_rows(
  overall_sims |>
    summarise(
      lo = q_lo(imr_black),
      hi = q_hi(imr_black),
      .by = year
    ) |>
    left_join(decomposition_results |> select(year, imr = imr_black), by = "year") |>
    mutate(race = "Non-Hispanic Black"),
  overall_sims |>
    summarise(
      lo = q_lo(imr_white),
      hi = q_hi(imr_white),
      .by = year
    ) |>
    left_join(decomposition_results |> select(year, imr = imr_white), by = "year") |>
    mutate(race = "Non-Hispanic White")
) |>
  select(year, race, imr, lo, hi) |>
  arrange(year, race)

write_csv(imr_with_ci, "data/analysis_data/imr_with_ci.csv")
cat("IMR CIs saved: imr_with_ci.csv\n")

# Per-GA decomposition CIs (for Fig 1B error bars on each GA group) ----
set.seed(853)
detailed_sims <- simulate_detailed_contributions(analysis_data)
detailed_ci <- summarize_detailed_ci(detailed_sims)

detailed_decomposition_with_ci <- detailed_results |>
  mutate(gestational_age_group = as.character(gestational_age_group)) |>
  left_join(detailed_ci, by = c("year", "gestational_age_group"))

write_csv(detailed_decomposition_with_ci,
          "data/analysis_data/detailed_decomposition_with_ci.csv")
cat("Per-GA decomposition CIs saved: detailed_decomposition_with_ci.csv\n")

# Cause-specific CIs ----
# Pre-arrange wide tables: per (year, cause), Black and White GA tibbles with
# deaths (cause-specific) and births (GA cell totals from analysis_data).
ga_cells <- analysis_data |>
  select(year, race, gestational_age_group, births, birth_proportion)

cause_with_cells <- cause_ga_data |>
  rename(cause_deaths = deaths, cause_births = births) |>
  left_join(ga_cells |> rename(cell_births = births),
            by = c("year", "race", "gestational_age_group")) |>
  mutate(births = coalesce(cell_births, cause_births)) |>
  select(year, race, gestational_age_group, cause_of_death,
         deaths = cause_deaths, births, birth_proportion)

# Ensure every (year, race, GA, cause) cell exists (fill missing with 0 deaths)
full_cells <- expand_grid(
  year = unique(cause_with_cells$year),
  race = c("Non-Hispanic Black", "Non-Hispanic White"),
  gestational_age_group = unique(cause_with_cells$gestational_age_group),
  cause_of_death = unique(cause_with_cells$cause_of_death)
) |>
  left_join(cause_with_cells,
            by = c("year", "race", "gestational_age_group", "cause_of_death")) |>
  left_join(ga_cells, by = c("year", "race", "gestational_age_group"),
            suffix = c("", ".cell")) |>
  mutate(
    deaths = replace_na(deaths, 0),
    births = coalesce(births, births.cell),
    birth_proportion = coalesce(birth_proportion, birth_proportion.cell)
  ) |>
  select(year, race, gestational_age_group, cause_of_death, deaths, births, birth_proportion)

set.seed(853)
cause_sims <- map_dfr(years, function(yr) {
  causes_yr <- full_cells |>
    filter(year == yr) |>
    distinct(cause_of_death) |>
    pull(cause_of_death)
  map_dfr(causes_yr, function(cse) {
    sub <- full_cells |> filter(year == yr, cause_of_death == cse)
    black <- sub |> filter(race == "Non-Hispanic Black") |>
      arrange(gestational_age_group)
    white <- sub |> filter(race == "Non-Hispanic White") |>
      arrange(gestational_age_group)
    simulate_kitagawa_cause_year(black, white) |>
      mutate(year = yr, cause_of_death = cse, .before = 1)
  })
})

cause_ci <- cause_sims |>
  summarise(
    composition_lo = q_lo(composition_effect),
    composition_hi = q_hi(composition_effect),
    rate_lo = q_lo(rate_effect),
    rate_hi = q_hi(rate_effect),
    gap_lo = q_lo(gap),
    gap_hi = q_hi(gap),
    .by = c(year, cause_of_death)
  )

cause_decomposition_with_ci <- cause_decomp |>
  left_join(cause_ci, by = c("year", "cause_of_death"))

write_csv(cause_decomposition_with_ci,
          "data/analysis_data/cause_decomposition_with_ci.csv")
cat("Cause-specific CIs saved: cause_decomposition_with_ci.csv\n")

#### Save outputs ####

# Data
write_csv(decomposition_results, "data/analysis_data/decomposition_results.csv")
write_csv(detailed_results, "data/analysis_data/detailed_decomposition.csv")
write_csv(decomp_by_age_results, "data/analysis_data/decomposition_by_mother_age.csv")
write_csv(cause_decomp, "data/analysis_data/cause_decomposition.csv")

cat("\n=== Output files ===\n")
cat("Data: data/analysis_data/decomposition_results.csv\n")
cat("Data: data/analysis_data/detailed_decomposition.csv\n")
cat("Data: data/analysis_data/decomposition_by_mother_age.csv\n")
cat("Data: data/analysis_data/cause_decomposition.csv\n")
cat("Data: data/analysis_data/decomposition_results_with_ci.csv\n")
cat("Data: data/analysis_data/detailed_decomposition_with_ci.csv\n")
cat("Data: data/analysis_data/cause_decomposition_with_ci.csv\n")
cat("Data: data/analysis_data/imr_with_ci.csv\n")
