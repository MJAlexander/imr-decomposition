#### Preamble ####
# Purpose: Sensitivity analyses for the brief report supplement:
#   (a) Alternative gestational-age decomposition schemes compared to the
#       7-category scheme used in the main report:
#         - 5-category: <32, 32-33, 34-36, 37-41, 42+ weeks (coarser)
#         - 16-category: <28, each week 28-41 individually, 42+ weeks (finer)
#   (b) Cause-of-death decomposition with the combined Sudden Unexpected
#       Infant Death (SUID) category broken out into its three ICD-10
#       components: SIDS (R95), accidental suffocation in bed (W75), and
#       ill-defined/unspecified (R99).
# Each output ships with Monte Carlo 95% CIs computed with the same Poisson
# resampling procedure as the main analysis.
# Author: Monica Alexander
# License: MIT
# Pre-requisites:
#   - scripts/01-clean_data.R has produced analysis_data_by_week.csv (for the
#     5- and 16-category GA re-aggregation) and cause_of_death_by_ga.csv (for
#     the SUID-split cause decomposition)
#   - scripts/02-analyze_data.R has produced the main-analysis outputs
#     (used here for side-by-side comparison) and defines the Monte Carlo
#     helper functions

#### Workspace setup ####
library(tidyverse)

dir.create("data/analysis_data", recursive = TRUE, showWarnings = FALSE)

set.seed(853)

# Reuse the Kitagawa + Monte Carlo machinery from the main analysis script.
# scripts/02-analyze_data.R is sourced for its function definitions; the
# data writes are repeated but harmless.
source("scripts/02-analyze_data.R")
source("scripts/categorize_cause.R")

#### GA-scheme decomposition helpers ####
# Re-aggregate the per-week birth/death counts into an arbitrary gestational-age
# scheme and run the overall + per-GA Kitagawa decompositions with Monte Carlo
# 95% CIs. Used for both the 5- and 16-category sensitivity schemes; structurally
# identical to the 7-category main analysis (02-analyze_data.R).

aggregate_ga_scheme <- function(data_by_week, categorize_fn, ga_levels) {
  data_by_week |>
    mutate(gestational_age_group = categorize_fn(gest_weeks)) |>
    filter(!is.na(gestational_age_group)) |>
    summarise(births = sum(births), deaths = sum(deaths),
              .by = c(year, race, gestational_age_group)) |>
    mutate(gestational_age_group = factor(gestational_age_group, levels = ga_levels)) |>
    mutate(
      mortality_rate = deaths / births * 1000,
      birth_proportion = births / sum(births),
      .by = c(year, race)
    ) |>
    arrange(year, race, gestational_age_group)
}

# Returns list(overall = decomposition_results + CI columns,
#              detailed = per-GA detailed decomposition + CI columns).
# Reseeds before each Monte Carlo block (matching the main-analysis blocks) so
# results are reproducible regardless of which schemes run before it.
decompose_ga_scheme <- function(analysis_data, seed = 853) {
  years_s <- sort(unique(analysis_data$year))

  overall <- map_dfr(years_s, ~kitagawa_decomposition(analysis_data, .x))
  set.seed(seed)
  overall_sims <- map_dfr(years_s, function(yr) {
    simulate_kitagawa_year(analysis_data |> filter(year == yr)) |>
      mutate(year = yr, .before = 1)
  })
  overall_ci <- overall_sims |>
    summarise(across(c(imr_black, imr_white, total_difference,
                       composition_effect, rate_effect, composition_pct, rate_pct),
                     list(lo = q_lo, hi = q_hi)),
              .by = year)
  overall_with_ci <- left_join(overall, overall_ci, by = "year")

  detailed <- map_dfr(years_s, ~detailed_decomposition(analysis_data, .x)) |>
    arrange(year, gestational_age_group) |>
    mutate(gestational_age_group = as.character(gestational_age_group))
  set.seed(seed)
  detailed_ci <- summarize_detailed_ci(simulate_detailed_contributions(analysis_data))
  detailed_with_ci <- left_join(detailed, detailed_ci,
                                by = c("year", "gestational_age_group"))

  list(overall = overall_with_ci, detailed = detailed_with_ci)
}

#### Define 5-category GA scheme (sensitivity) ####
# The coarser clinical scheme used in earlier drafts: very preterm (<32),
# early preterm (32-33), late preterm (34-36), term (37-41), post-term (42+).
ga_levels_5 <- c("<32 weeks", "32-33 weeks", "34-36 weeks", "37-41 weeks",
                 "42+ weeks")

categorize_gest_age_5 <- function(weeks) {
  case_when(
    is.na(weeks) | weeks == 99 ~ NA_character_,
    weeks < 32 ~ "<32 weeks",
    weeks >= 32 & weeks <= 33 ~ "32-33 weeks",
    weeks >= 34 & weeks <= 36 ~ "34-36 weeks",
    weeks >= 37 & weeks <= 41 ~ "37-41 weeks",
    weeks >= 42 ~ "42+ weeks",
    TRUE ~ NA_character_
  )
}

#### Define 16-category GA scheme (sensitivity) ####
# <28 weeks, each week 28-41 individually, and 42+ weeks (16 groups).
ga_levels_16 <- c("<28 weeks", paste0(28:41, " weeks"), "42+ weeks")

categorize_gest_age_16 <- function(weeks) {
  case_when(
    is.na(weeks) | weeks == 99 ~ NA_character_,
    weeks < 28 ~ "<28 weeks",
    weeks >= 28 & weeks <= 41 ~ paste0(weeks, " weeks"),
    weeks >= 42 ~ "42+ weeks",
    TRUE ~ NA_character_
  )
}

#### 5- and 16-category sensitivity decompositions ####
# The main analysis (02-analyze_data.R) uses the 7-category scheme. Here we
# re-aggregate the same per-week counts into the coarser 5-category and finer
# 16-category schemes to show the composition/rate split is scheme-dependent.

data_by_week <- read_csv("data/analysis_data/analysis_data_by_week.csv",
                          show_col_types = FALSE) |>
  filter(year >= 2003)

analysis_data_5 <- aggregate_ga_scheme(data_by_week, categorize_gest_age_5, ga_levels_5)
res_5 <- decompose_ga_scheme(analysis_data_5)
write_csv(res_5$overall, "data/analysis_data/decomposition_results_5cat.csv")
write_csv(res_5$detailed, "data/analysis_data/detailed_decomposition_5cat.csv")

analysis_data_16 <- aggregate_ga_scheme(data_by_week, categorize_gest_age_16, ga_levels_16)
res_16 <- decompose_ga_scheme(analysis_data_16)
write_csv(res_16$overall, "data/analysis_data/decomposition_results_16cat.csv")
write_csv(res_16$detailed, "data/analysis_data/detailed_decomposition_16cat.csv")

#### Cause-specific decomposition (SUID split, main 7-category scheme) ####

#### GA-scheme comparison (7-category main vs 5- and 16-category sensitivities) ####

scheme_cols <- function(d, label) {
  d |> transmute(year, scheme = label, composition_effect, rate_effect,
                 total_difference, composition_pct, rate_pct)
}

ga_scheme_comparison <- bind_rows(
  scheme_cols(decomposition_results, "7-category"),  # main, from sourcing 02
  scheme_cols(res_5$overall, "5-category"),
  scheme_cols(res_16$overall, "16-category")
) |>
  arrange(year, scheme)

write_csv(ga_scheme_comparison,
          "data/analysis_data/ga_scheme_comparison.csv")

#### SUID-split cause decomposition (main 7-category GA scheme) ####
# Uses the main 7-category GA scheme and the SUID-split categorizer.

cause_ga_raw_suid <- read_csv("data/analysis_data/cause_of_death_by_ga.csv",
                              show_col_types = FALSE) |>
  filter(year >= 2003)

cause_ga_suid <- cause_ga_raw_suid |>
  mutate(cause_of_death = categorize_cause_suid_split(ucod)) |>
  summarise(deaths = sum(deaths), births = first(births),
            .by = c(year, race, gestational_age_group, cause_of_death))

birth_props_main <- analysis_data |>
  select(year, race, gestational_age_group, birth_proportion, births_total = births)

cause_ga_suid_rates <- cause_ga_suid |>
  left_join(birth_props_main, by = c("year", "race", "gestational_age_group")) |>
  mutate(births = coalesce(births_total, births),
         rate = deaths / births * 1000) |>
  select(year, race, gestational_age_group, cause_of_death,
         deaths, births, birth_proportion, rate)

all_combos_suid <- expand_grid(
  year = unique(cause_ga_suid_rates$year),
  race = c("Non-Hispanic Black", "Non-Hispanic White"),
  gestational_age_group = unique(cause_ga_suid_rates$gestational_age_group),
  cause_of_death = unique(cause_ga_suid_rates$cause_of_death)
)

cause_ga_suid_full <- all_combos_suid |>
  left_join(cause_ga_suid_rates |>
              select(year, race, gestational_age_group, cause_of_death, rate),
            by = c("year", "race", "gestational_age_group", "cause_of_death")) |>
  mutate(rate = replace_na(rate, 0)) |>
  left_join(birth_props_main |> select(year, race, gestational_age_group, birth_proportion, births_total),
            by = c("year", "race", "gestational_age_group"))

cause_decomp_suid_split <- cause_ga_suid_full |>
  pivot_wider(names_from = race,
              values_from = c(rate, birth_proportion, births_total),
              names_sep = "_") |>
  rename(
    M_B = `rate_Non-Hispanic Black`,
    M_W = `rate_Non-Hispanic White`,
    p_B = `birth_proportion_Non-Hispanic Black`,
    p_W = `birth_proportion_Non-Hispanic White`
  ) |>
  mutate(
    comp_i = (p_B - p_W) * (M_B + M_W) / 2,
    rate_i = (M_B - M_W) * (p_B + p_W) / 2
  ) |>
  summarise(composition_effect = sum(comp_i),
            rate_effect = sum(rate_i),
            .by = c(year, cause_of_death)) |>
  mutate(gap = composition_effect + rate_effect)

# Monte Carlo CIs for SUID-split decomposition
ga_cells_main <- analysis_data |>
  select(year, race, gestational_age_group, births, birth_proportion)

cause_suid_cells_full <- expand_grid(
  year = unique(cause_ga_suid$year),
  race = c("Non-Hispanic Black", "Non-Hispanic White"),
  gestational_age_group = unique(ga_cells_main$gestational_age_group),
  cause_of_death = unique(cause_ga_suid$cause_of_death)
) |>
  left_join(cause_ga_suid |>
              select(year, race, gestational_age_group, cause_of_death, deaths),
            by = c("year", "race", "gestational_age_group", "cause_of_death")) |>
  left_join(ga_cells_main,
            by = c("year", "race", "gestational_age_group")) |>
  mutate(deaths = replace_na(deaths, 0))

set.seed(853)
suid_cause_sims <- map_dfr(sort(unique(cause_suid_cells_full$year)), function(yr) {
  causes_yr <- cause_suid_cells_full |>
    filter(year == yr) |>
    distinct(cause_of_death) |>
    pull(cause_of_death)
  map_dfr(causes_yr, function(cse) {
    sub <- cause_suid_cells_full |> filter(year == yr, cause_of_death == cse)
    black <- sub |> filter(race == "Non-Hispanic Black") |>
      arrange(gestational_age_group)
    white <- sub |> filter(race == "Non-Hispanic White") |>
      arrange(gestational_age_group)
    simulate_kitagawa_cause_year(black, white) |>
      mutate(year = yr, cause_of_death = cse, .before = 1)
  })
})

suid_cause_ci <- suid_cause_sims |>
  summarise(
    composition_lo = q_lo(composition_effect),
    composition_hi = q_hi(composition_effect),
    rate_lo = q_lo(rate_effect),
    rate_hi = q_hi(rate_effect),
    gap_lo = q_lo(gap),
    gap_hi = q_hi(gap),
    .by = c(year, cause_of_death)
  )

cause_decomp_suid_split_with_ci <- cause_decomp_suid_split |>
  left_join(suid_cause_ci, by = c("year", "cause_of_death"))

write_csv(cause_decomp_suid_split_with_ci,
          "data/analysis_data/cause_decomposition_suid_split.csv")

cat("\n=== Sensitivity Analysis Outputs ===\n")
cat("Data: data/analysis_data/decomposition_results_5cat.csv\n")
cat("Data: data/analysis_data/detailed_decomposition_5cat.csv\n")
cat("Data: data/analysis_data/decomposition_results_16cat.csv\n")
cat("Data: data/analysis_data/detailed_decomposition_16cat.csv\n")
cat("Data: data/analysis_data/ga_scheme_comparison.csv\n")
cat("Data: data/analysis_data/cause_decomposition_suid_split.csv\n")
