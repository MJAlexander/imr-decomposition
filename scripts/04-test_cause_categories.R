#### Preamble ####
# Purpose: Tests that cause-of-death categories match Wolf et al. (2025,
#          JAMA Pediatrics) eTable 1 and that no old category names remain
# Author: Monica Alexander
# Date: 2026-02-13
# Contact: monica.alexander@utoronto.ca
# License: MIT
# Pre-requisites:
#   - Run 02-analyze_data.R first to create cause_decomposition.csv
#   - Required R packages: tidyverse

#### Workspace setup ####
library(tidyverse)

pass <- 0
fail <- 0

check <- function(description, condition) {
  if (isTRUE(condition)) {
    cat("  PASS:", description, "\n")
    pass <<- pass + 1
  } else {
    cat("  FAIL:", description, "\n")
    fail <<- fail + 1
  }
}

cat("=== Test: Cause-of-death categories (Wolf et al. 2025) ===\n\n")

#### Test 1: categorize_cause ICD-10 mappings ####
cat("Test 1: ICD-10 code mappings\n")

# Source the shared function
source("scripts/categorize_cause.R")

# 1. Congenital malformations: Q00-Q99
check("Q00 -> Congenital malformations", categorize_cause("Q000") == "Congenital malformations")
check("Q99 -> Congenital malformations", categorize_cause("Q999") == "Congenital malformations")
check("Q24 -> Congenital malformations", categorize_cause("Q248") == "Congenital malformations")

# 2. Short gestation/LBW: P07
check("P07 -> Short gestation/LBW", categorize_cause("P070") == "Short gestation/LBW")
check("P07 -> Short gestation/LBW", categorize_cause("P073") == "Short gestation/LBW")

# 3. SUID: R95, R99, W75
check("R95 (SIDS) -> SUID", categorize_cause("R950") == "Sudden unexpected infant death")
check("R99 (ill-defined) -> SUID", categorize_cause("R990") == "Sudden unexpected infant death")
check("W75 (suffocation in bed) -> SUID", categorize_cause("W750") == "Sudden unexpected infant death")
check("R95 without 4th char -> SUID", categorize_cause("R95") == "Sudden unexpected infant death")

# 4. Unintentional injuries: V01-V99, W00-W74, W76-W99, X00-X59
check("V01 -> Unintentional injuries", categorize_cause("V010") == "Unintentional injuries")
check("V99 -> Unintentional injuries", categorize_cause("V990") == "Unintentional injuries")
check("W00 -> Unintentional injuries", categorize_cause("W000") == "Unintentional injuries")
check("W74 -> Unintentional injuries", categorize_cause("W740") == "Unintentional injuries")
check("W76 -> Unintentional injuries", categorize_cause("W760") == "Unintentional injuries")
check("W99 -> Unintentional injuries", categorize_cause("W990") == "Unintentional injuries")
check("X00 -> Unintentional injuries", categorize_cause("X000") == "Unintentional injuries")
check("X59 -> Unintentional injuries", categorize_cause("X590") == "Unintentional injuries")
# W75 should NOT be unintentional injuries (it's SUID)
check("W75 is NOT Unintentional injuries", categorize_cause("W750") != "Unintentional injuries")

# 5. Maternal complications: P01
check("P01 -> Maternal complications", categorize_cause("P010") == "Maternal complications")

# 6. Placenta/cord/membranes: P02
check("P02 -> Placenta/cord/membranes", categorize_cause("P020") == "Placenta/cord/membranes")

# 7. Infection: A00-B99, P35-P39 (includes bacterial sepsis P36)
check("A00 -> Infection", categorize_cause("A000") == "Infection")
check("A41 -> Infection", categorize_cause("A410") == "Infection")
check("B99 -> Infection", categorize_cause("B990") == "Infection")
check("P35 -> Infection", categorize_cause("P350") == "Infection")
check("P36 (bacterial sepsis) -> Infection", categorize_cause("P360") == "Infection")
check("P37 -> Infection", categorize_cause("P370") == "Infection")
check("P38 -> Infection", categorize_cause("P380") == "Infection")
check("P39 -> Infection", categorize_cause("P390") == "Infection")
# P36 should NOT be its own category
check("P36 is NOT 'Bacterial sepsis of newborn'",
      categorize_cause("P360") != "Bacterial sepsis of newborn")

# 8. Respiratory distress: P22
check("P22 -> Respiratory distress", categorize_cause("P220") == "Respiratory distress")

# 9. Circulatory diseases: I00-I99
check("I00 -> Circulatory diseases", categorize_cause("I000") == "Circulatory diseases")
check("I99 -> Circulatory diseases", categorize_cause("I990") == "Circulatory diseases")

# 10. Hemorrhage: P50-P52, P54 (NOT P53)
check("P50 -> Hemorrhage", categorize_cause("P500") == "Hemorrhage")
check("P51 -> Hemorrhage", categorize_cause("P510") == "Hemorrhage")
check("P52 -> Hemorrhage", categorize_cause("P520") == "Hemorrhage")
check("P54 -> Hemorrhage", categorize_cause("P540") == "Hemorrhage")
check("P53 is NOT Hemorrhage", categorize_cause("P530") != "Hemorrhage")
check("P53 -> Other", categorize_cause("P530") == "Other")
# Should NOT use old name "Neonatal hemorrhage"
check("P50 is NOT 'Neonatal hemorrhage'", categorize_cause("P500") != "Neonatal hemorrhage")

# 11. Hypoxia and birth asphyxia: P20-P21
check("P20 -> Hypoxia and birth asphyxia", categorize_cause("P200") == "Hypoxia and birth asphyxia")
check("P21 -> Hypoxia and birth asphyxia", categorize_cause("P210") == "Hypoxia and birth asphyxia")
# Should NOT use old name "Birth asphyxia/hypoxia"
check("P20 is NOT 'Birth asphyxia/hypoxia'",
      categorize_cause("P200") != "Birth asphyxia/hypoxia")

# 12. Necrotizing enterocolitis: P77
check("P77 -> Necrotizing enterocolitis", categorize_cause("P770") == "Necrotizing enterocolitis")

# 13. Assault (homicide): U01, X85-Y09
check("U01 -> Assault (homicide)", categorize_cause("U010") == "Assault (homicide)")
check("X85 -> Assault (homicide)", categorize_cause("X850") == "Assault (homicide)")
check("Y09 -> Assault (homicide)", categorize_cause("Y090") == "Assault (homicide)")
check("Y01 -> Assault (homicide)", categorize_cause("Y010") == "Assault (homicide)")
# Should NOT be just "Assault" without "(homicide)"
check("X85 is NOT plain 'Assault'", categorize_cause("X850") != "Assault")

# 14. Old categories that should no longer exist
check("R99 is NOT 'Ill-defined/unspecified'",
      categorize_cause("R990") != "Ill-defined/unspecified")
check("R95 is NOT 'SIDS'", categorize_cause("R950") != "SIDS")
check("P29 -> Other (not 'Perinatal cardiovascular')",
      categorize_cause("P290") == "Other")
check("P25 -> Other (not 'Perinatal respiratory')",
      categorize_cause("P250") == "Other")
check("P23 -> Other (not separate 'Infections' category)",
      categorize_cause("P230") == "Other")
check("J18 -> Other (not 'Infections')",
      categorize_cause("J180") == "Other")

cat("\n")

#### Test 2: Analysis scripts source the shared categorize_cause.R ####
cat("Test 2: categorize_cause() sourced from shared file\n")

# The canonical function lives in scripts/categorize_cause.R and is sourced by
# the analysis scripts (no inline copy). The brief report reads pre-computed
# CSVs and does not define categorize_cause(), so there is no manuscript copy
# to keep in sync.
for (script in c("scripts/02-analyze_data.R", "scripts/03-sensitivity_analyses.R")) {
  script_lines <- readLines(script)
  check(paste0(basename(script), " sources categorize_cause.R"),
        any(str_detect(script_lines, 'source\\("scripts/categorize_cause\\.R"\\)')))
}

cat("\n")

#### Test 3: Output CSV has only valid category names ####
cat("Test 3: cause_decomposition.csv category names\n")

cause_csv <- read_csv("data/analysis_data/cause_decomposition.csv",
                       show_col_types = FALSE)

valid_categories <- c(
  "Congenital malformations",
  "Short gestation/LBW",
  "Sudden unexpected infant death",
  "Unintentional injuries",
  "Maternal complications",
  "Placenta/cord/membranes",
  "Infection",
  "Respiratory distress",
  "Circulatory diseases",
  "Hemorrhage",
  "Hypoxia and birth asphyxia",
  "Necrotizing enterocolitis",
  "Assault (homicide)",
  "Other"
)

old_categories <- c(
  "SIDS",
  "Ill-defined/unspecified",
  "Bacterial sepsis of newborn",
  "Perinatal cardiovascular",
  "Perinatal respiratory (other)",
  "Neonatal hemorrhage",
  "Birth asphyxia/hypoxia",
  "Assault",
  "Infections"
)

actual_categories <- sort(unique(cause_csv$cause_of_death))

check("All categories in CSV are valid",
      all(actual_categories %in% valid_categories))

check("No old category names in CSV",
      !any(actual_categories %in% old_categories))

unexpected <- setdiff(actual_categories, valid_categories)
if (length(unexpected) > 0) {
  cat("    Unexpected categories found:", paste(unexpected, collapse = ", "), "\n")
}

missing <- setdiff(valid_categories, actual_categories)
if (length(missing) > 0) {
  cat("    Categories not in data:", paste(missing, collapse = ", "), "\n")
}

check("CSV has exactly 14 categories", length(actual_categories) == 14)

cat("\n")

#### Test 4: No old category names in README prose ####
cat("Test 4: No old category names in README\n")

readme_text <- paste(readLines("README.md"), collapse = "\n")
check("No 'SIDS' in README", !str_detect(readme_text, "\\bSIDS\\b"))
check("No 'Ill-defined' in README", !str_detect(readme_text, "Ill-defined"))

cat("\n")

#### Test 5: Kitagawa columns exist in CSV ####
cat("Test 5: Kitagawa decomposition columns exist\n")

check("composition_effect column exists",
      "composition_effect" %in% names(cause_csv))
check("rate_effect column exists",
      "rate_effect" %in% names(cause_csv))
check("composition_pct column exists",
      "composition_pct" %in% names(cause_csv))
check("rate_pct column exists",
      "rate_pct" %in% names(cause_csv))
check("gap column exists",
      "gap" %in% names(cause_csv))
check("csimr_black column exists",
      "csimr_black" %in% names(cause_csv))
check("csimr_white column exists",
      "csimr_white" %in% names(cause_csv))

cat("\n")

#### Test 6: composition_effect + rate_effect ≈ gap ####
cat("Test 6: composition_effect + rate_effect ≈ gap for each row\n")

cause_csv <- cause_csv |>
  mutate(check_sum = composition_effect + rate_effect)

max_diff <- max(abs(cause_csv$check_sum - cause_csv$gap))
check(paste0("Max |comp + rate - gap| = ", round(max_diff, 8), " (should be ~0)"),
      max_diff < 1e-6)

for (yr in c(2003, 2010, 2020, 2024)) {
  yr_diff <- cause_csv |>
    filter(year == yr) |>
    mutate(d = abs(check_sum - gap)) |>
    pull(d) |>
    max()
  check(paste0("comp + rate = gap (", yr, "), max diff = ", round(yr_diff, 8)),
        yr_diff < 1e-6)
}

cat("\n")

#### Test 7: Cause-specific gaps sum to total gap ####
cat("Test 7: Cause-specific gaps sum to total gap\n")
cat("  (Exact additivity expected within floating-point precision)\n")

decomp_results <- read_csv("data/analysis_data/decomposition_results.csv",
                            show_col_types = FALSE)

for (yr in c(2003, 2010, 2020, 2024)) {
  cause_gap_sum <- cause_csv |>
    filter(year == yr) |>
    pull(gap) |>
    sum()
  total_gap <- decomp_results |>
    filter(year == yr) |>
    pull(total_difference)
  check(paste0("Cause gaps = total gap (", yr, "), diff = ",
               format(abs(cause_gap_sum - total_gap), scientific = TRUE)),
        abs(cause_gap_sum - total_gap) < 1e-10)
}

cat("\n")

#### Test 8: Sum of cause effects = overall effects ####
cat("Test 8: Sum of cause-specific effects = overall Kitagawa effects\n")
cat("  (Exact additivity expected within floating-point precision)\n")

for (yr in c(2003, 2010, 2020, 2024)) {
  cause_comp_sum <- cause_csv |>
    filter(year == yr) |>
    pull(composition_effect) |>
    sum()
  cause_rate_sum <- cause_csv |>
    filter(year == yr) |>
    pull(rate_effect) |>
    sum()
  overall_comp <- decomp_results |>
    filter(year == yr) |>
    pull(composition_effect)
  overall_rate <- decomp_results |>
    filter(year == yr) |>
    pull(rate_effect)
  check(paste0("Cause comp sum = overall comp (", yr, "), diff = ",
               format(abs(cause_comp_sum - overall_comp), scientific = TRUE)),
        abs(cause_comp_sum - overall_comp) < 1e-10)
  check(paste0("Cause rate sum = overall rate (", yr, "), diff = ",
               format(abs(cause_rate_sum - overall_rate), scientific = TRUE)),
        abs(cause_rate_sum - overall_rate) < 1e-10)
}

cat("\n")

#### Test 9: categorize_cause_suid_split() splits SUID correctly ####
cat("Test 9: SUID-split categorizer\n")

check("R95 -> 'SIDS (R95)'",
      categorize_cause_suid_split("R950") == "SIDS (R95)")
check("W75 -> 'Accidental suffocation in bed (W75)'",
      categorize_cause_suid_split("W750") == "Accidental suffocation in bed (W75)")
check("R99 -> 'Ill-defined/unspecified (R99)'",
      categorize_cause_suid_split("R990") == "Ill-defined/unspecified (R99)")
check("Q24 unchanged in suid_split",
      categorize_cause_suid_split("Q248") == "Congenital malformations")
check("P07 unchanged in suid_split",
      categorize_cause_suid_split("P073") == "Short gestation/LBW")
check("suid_split never produces combined SUID label",
      !any(c(categorize_cause_suid_split("R950"),
             categorize_cause_suid_split("W750"),
             categorize_cause_suid_split("R990")) == "Sudden unexpected infant death"))

if (file.exists("data/analysis_data/cause_decomposition_suid_split.csv")) {
  suid_split_csv <- read_csv(
    "data/analysis_data/cause_decomposition_suid_split.csv",
    show_col_types = FALSE
  )
  cause_csv_main <- read_csv(
    "data/analysis_data/cause_decomposition.csv",
    show_col_types = FALSE
  )

  for (yr in c(2003, 2020, 2024)) {
    suid_components <- suid_split_csv |>
      filter(year == yr,
             cause_of_death %in% c("SIDS (R95)",
                                    "Accidental suffocation in bed (W75)",
                                    "Ill-defined/unspecified (R99)"))
    suid_sum_gap <- sum(suid_components$gap)
    suid_sum_comp <- sum(suid_components$composition_effect)
    suid_sum_rate <- sum(suid_components$rate_effect)

    combined <- cause_csv_main |>
      filter(year == yr, cause_of_death == "Sudden unexpected infant death")
    if (nrow(combined) == 1) {
      check(paste0("SUID-split gap sum = combined SUID gap (", yr, ")"),
            abs(suid_sum_gap - combined$gap) < 1e-8)
      check(paste0("SUID-split composition sum = combined SUID composition (", yr, ")"),
            abs(suid_sum_comp - combined$composition_effect) < 1e-8)
      check(paste0("SUID-split rate sum = combined SUID rate (", yr, ")"),
            abs(suid_sum_rate - combined$rate_effect) < 1e-8)
    }
  }
} else {
  cat("  (skipped: cause_decomposition_suid_split.csv not yet generated)\n")
}

cat("\n")

#### Test 10: sensitivity GA-scheme decompositions match the main total ####
cat("Test 10: sensitivity GA-scheme decompositions (vs main 7-category)\n")

# Each sensitivity scheme is a pure re-aggregation of the same births/deaths,
# so the race-specific IMRs and the total gap must be identical to the main
# 7-category analysis; only the composition/rate split changes.
if (file.exists("data/analysis_data/decomposition_results.csv")) {
  decomp_main <- read_csv("data/analysis_data/decomposition_results.csv",
                          show_col_types = FALSE)
  for (scheme in c("5cat", "16cat")) {
    f <- paste0("data/analysis_data/decomposition_results_", scheme, ".csv")
    # Assert existence as a check (not a silent `next`): a missing output is a
    # failure, not a skip, so a half-run pipeline cannot report all-passed.
    f_exists <- file.exists(f)
    check(paste0(scheme, " output file exists"), f_exists)
    if (!f_exists) next
    d <- read_csv(f, show_col_types = FALSE)
    joined <- decomp_main |>
      select(year, imr_b_m = imr_black, imr_w_m = imr_white,
             total_m = total_difference) |>
      inner_join(d |> select(year, imr_b = imr_black, imr_w = imr_white,
                             total = total_difference), by = "year")
    # nrow guard: an empty or year-disjoint CSV would otherwise pass vacuously
    # because max(numeric(0)) == -Inf < tol.
    check(paste0(scheme, " has matching years"), nrow(joined) == nrow(decomp_main))
    check(paste0(scheme, " IMR Black matches main (tol 1e-6)"),
          nrow(joined) > 0 && max(abs(joined$imr_b_m - joined$imr_b)) < 1e-6)
    check(paste0(scheme, " IMR White matches main (tol 1e-6)"),
          nrow(joined) > 0 && max(abs(joined$imr_w_m - joined$imr_w)) < 1e-6)
    check(paste0(scheme, " total gap matches main (tol 1e-6)"),
          nrow(joined) > 0 && max(abs(joined$total_m - joined$total)) < 1e-6)
    check(paste0(scheme, ": composition + rate = total (tol 1e-8)"),
          nrow(d) > 0 &&
            max(abs((d$composition_effect + d$rate_effect) - d$total_difference)) < 1e-8)
  }
} else {
  cat("  (skipped: main decomposition not yet generated)\n")
}

cat("\n")

#### Test 11: Monte Carlo CIs contain the point estimates ####
cat("Test 11: Monte Carlo CIs contain the point estimates\n")

if (file.exists("data/analysis_data/decomposition_results_with_ci.csv")) {
  ci_csv <- read_csv("data/analysis_data/decomposition_results_with_ci.csv",
                      show_col_types = FALSE)
  # The point estimate is the deterministic Kitagawa value; the CI is the
  # 2.5/97.5% quantiles of a Poisson resample centred on the observed counts.
  # For these aggregate (sum/difference) metrics the point estimate must lie
  # inside its own CI in every year -- a loose threshold would hide a real
  # mis-pairing of point estimates and CI bounds.
  for (metric in c("imr_black", "imr_white", "total_difference",
                   "composition_effect", "rate_effect")) {
    lo <- ci_csv[[paste0(metric, "_lo")]]
    hi <- ci_csv[[paste0(metric, "_hi")]]
    pt <- ci_csv[[metric]]
    n_in <- sum(pt >= lo & pt <= hi, na.rm = TRUE)
    check(paste0("Point inside 95% CI for ", metric,
                 " (", n_in, "/", length(pt), " years)"),
          n_in == length(pt))
  }
} else {
  cat("  (skipped: decomposition_results_with_ci.csv not yet generated)\n")
}

cat("\n")

#### Test 12: Alexander 1996 birthweight exclusion column exists ####
cat("Test 12: exclusion_counts.csv has Alexander column\n")

if (file.exists("data/analysis_data/exclusion_counts.csv")) {
  excl_csv <- read_csv("data/analysis_data/exclusion_counts.csv",
                        show_col_types = FALSE)
  check("excluded_implausible_bwt_ga column exists",
        "excluded_implausible_bwt_ga" %in% names(excl_csv))
  if ("excluded_implausible_bwt_ga" %in% names(excl_csv)) {
    n_excl_2003plus <- excl_csv |>
      filter(year >= 2003) |>
      pull(excluded_implausible_bwt_ga) |>
      sum()
    check("Alexander exclusion > 0 post-2003",
          n_excl_2003plus > 0)
  }
  # Unknown/missing Hispanic origin is an exclusion reason (mhisp != 0); the
  # dropped White/Black-race records are a subset of the race/ethnicity exclusion.
  check("excluded_hispanic_unknown column exists",
        "excluded_hispanic_unknown" %in% names(excl_csv))
  if ("excluded_hispanic_unknown" %in% names(excl_csv)) {
    check("unknown-Hispanic exclusion is a subset of race exclusion",
          all(excl_csv$excluded_hispanic_unknown <= excl_csv$excluded_race))
  }
  # Partition identity: missing + low + high == total Alexander exclusion.
  # (GA is no longer an exclusion reason, so there is no ga_out_of_range term;
  # this identity would fail if any ga_known record were mis-bucketed.)
  partition_cols <- c("excluded_bwt_missing", "excluded_bwt_low",
                      "excluded_bwt_high")
  for (col in partition_cols) {
    check(paste0(col, " column exists"), col %in% names(excl_csv))
  }
  check("excluded_bwt_ga_out_of_range column was removed",
        !("excluded_bwt_ga_out_of_range" %in% names(excl_csv)))
  # The Talge z-score comparison was removed; none of its columns should remain.
  zscore_cols <- c("excluded_zscore", "excluded_zscore_ineligible",
                   "excluded_alexander_only_eligible", "excluded_zscore_only_eligible",
                   "excluded_both_methods_eligible", "n_zscore_eligible")
  check("z-score comparison columns were removed",
        !any(zscore_cols %in% names(excl_csv)))
  if (all(c("excluded_implausible_bwt_ga", partition_cols) %in% names(excl_csv))) {
    yr_check <- excl_csv |>
      filter(year >= 2003) |>
      mutate(diff = excluded_implausible_bwt_ga -
               (excluded_bwt_missing + excluded_bwt_low + excluded_bwt_high))
    check("Alexander = missing + low + high (partition identity)",
          all(yr_check$diff == 0))
  }
} else {
  cat("  (skipped: exclusion_counts.csv not yet generated)\n")
}

cat("\n")

#### Test 13: detailed decompositions exist and are additive ####
cat("Test 13: detailed decomposition (main 7-category + 5/16-cat sensitivities)\n")

# (label, detailed CSV, overall CSV): the main 7-category analysis plus the two
# sensitivity schemes. All must reconstruct the Kitagawa formula and sum to the
# scheme's total gap.
detailed_specs <- list(
  c("main 7cat", "data/analysis_data/detailed_decomposition.csv",
    "data/analysis_data/decomposition_results.csv"),
  c("5cat", "data/analysis_data/detailed_decomposition_5cat.csv",
    "data/analysis_data/decomposition_results_5cat.csv"),
  c("16cat", "data/analysis_data/detailed_decomposition_16cat.csv",
    "data/analysis_data/decomposition_results_16cat.csv")
)
for (spec in detailed_specs) {
  scheme <- spec[1]; f <- spec[2]; overall_f <- spec[3]
  # A missing detailed output is a failure, not a silent skip.
  f_exists <- file.exists(f)
  check(paste0(scheme, " detailed output file exists"), f_exists)
  if (!f_exists) next
  det <- read_csv(f, show_col_types = FALSE)
  for (col in c("year", "gestational_age_group",
                "prop_black", "prop_white",
                "mortality_black", "mortality_white",
                "composition_contribution", "rate_contribution",
                "total_contribution")) {
    check(paste0(scheme, " detailed has ", col), col %in% names(det))
  }

  # Reconstruct each row's contributions from the per-GA Kitagawa formula to
  # verify the formula, not just CSV serialization fidelity.
  det_chk <- det |>
    mutate(
      comp_recomp = (prop_black - prop_white) *
        (mortality_black + mortality_white) / 2,
      rate_recomp = (mortality_black - mortality_white) *
        (prop_black + prop_white) / 2
    )
  check(paste0(scheme, " composition matches Kitagawa formula (tol 1e-10)"),
        nrow(det_chk) > 0 &&
          max(abs(det_chk$comp_recomp - det_chk$composition_contribution)) < 1e-10)
  check(paste0(scheme, " rate matches Kitagawa formula (tol 1e-10)"),
        nrow(det_chk) > 0 &&
          max(abs(det_chk$rate_recomp - det_chk$rate_contribution)) < 1e-10)

  # Year-level: sum of per-GA contributions == total gap for that scheme.
  if (file.exists(overall_f)) {
    overall <- read_csv(overall_f, show_col_types = FALSE) |>
      select(year, total_difference)
    year_check <- det |>
      summarise(total_sum = sum(total_contribution), .by = year) |>
      inner_join(overall, by = "year") |>
      mutate(diff = abs(total_sum - total_difference))
    check(paste0(scheme, " per-GA sums = total gap (year-level, tol 1e-8)"),
          nrow(year_check) > 0 && max(year_check$diff) < 1e-8)
  }
}

cat("\n")

#### Test 14: Birthweight helper module (Alexander 1996) ####
cat("Test 14: birthweight_helpers.R (Alexander)\n")

# The helpers live in a dedicated, sourceable module (mirrors categorize_cause.R),
# so the test exercises the exact functions the pipeline uses.
source("scripts/birthweight_helpers.R")

# alexander_reason() single-source classifier (drives both the flag and the
# exclusion-count partition).
check("alexander_reason: normal term -> plausible",
      alexander_reason(3000L, 40L) == "plausible")
check("alexander_reason: NA bwt -> missing_bwt",
      alexander_reason(NA_integer_, 40L) == "missing_bwt")
check("alexander_reason: 124g at 20wk -> implausibly_low",
      alexander_reason(124L, 20L) == "implausibly_low")
check("alexander_reason: 1251g at 20wk -> implausibly_high",
      alexander_reason(1251L, 20L) == "implausibly_high")
check("alexander_reason: 6000g at >=38wk -> plausible (boundary)",
      alexander_reason(6000L, 40L) == "plausible")
# Extrapolated boundary rows: GA outside 20-47 is NOT excluded for GA; instead
# the nearest boundary bound is applied (<20 wk -> 20-wk bound; >=47 -> >=38).
check("alexander_reason: 500g at GA 18 (<20) -> plausible (20-wk bound)",
      alexander_reason(500L, 18L) == "plausible")
check("alexander_reason: 3000g at GA 18 (<20) -> implausibly_high (>1250)",
      alexander_reason(3000L, 18L) == "implausibly_high")
check("alexander_reason: 3000g at GA 50 (>=47) -> plausible (>=38 bound)",
      alexander_reason(3000L, 50L) == "plausible")
check("alexander_reason: 500g at GA 50 (>=47) -> implausibly_low (<1000)",
      alexander_reason(500L, 50L) == "implausibly_low")
# With extrapolated bounds, no non-NA gestational age is ever 'ga_out_of_range'.
ga_grid <- 15:55
check("alexander_reason: no in-data GA yields ga_out_of_range",
      !any(alexander_reason(3000L, ga_grid) == "ga_out_of_range"))
# alexander_reason only ever returns documented levels (matches the values the
# cleaning pipeline maps to exclusion-count columns).
grid <- expand.grid(bwt = c(NA, 100L, 500L, 3000L, 7000L),
                    weeks = c(NA, 19L, 21L, 32L, 40L, 45L))
check("alexander_reason returns only documented levels",
      all(alexander_reason(grid$bwt, grid$weeks) %in% alexander_reason_levels))

cat("\n")

#### Test 15: pipeline IMRs validate against published CDC/NCHS rates ####
cat("Test 15: pipeline IMRs vs published CDC reference (within tolerance)\n")

cdc_path <- "data/analysis_data/cdc_reference_imr.csv"
imr_path <- "data/analysis_data/imr_with_ci.csv"
cdc_exists <- file.exists(cdc_path)
imr_exists <- file.exists(imr_path)
check("cdc_reference_imr.csv exists", cdc_exists)
check("imr_with_ci.csv exists", imr_exists)

if (cdc_exists && imr_exists) {
  cdc_ref <- read_csv(cdc_path, show_col_types = FALSE)
  imr_wide <- read_csv(imr_path, show_col_types = FALSE) |>
    mutate(race = recode(race,
                         "Non-Hispanic Black" = "black",
                         "Non-Hispanic White" = "white")) |>
    select(year, race, imr) |>
    pivot_wider(names_from = race, values_from = imr)

  cmp <- inner_join(imr_wide, cdc_ref, by = "year")

  # Guard against a silently empty or mis-joined frame: all(<empty>) is TRUE,
  # so without this the tolerance checks below would pass vacuously if the join
  # produced no rows or the race columns failed to materialise (e.g. an upstream
  # race-label change slipping past recode()). Gate every downstream assertion
  # on cols_ok so a broken join fails the suite instead of greening it.
  cols_ok <- all(c("black", "white") %in% names(cmp)) && nrow(cmp) >= 20

  # The pipeline applies the Alexander birthweight/GA exclusion and restricts
  # to records with known gestational age, so its IMRs run slightly below the
  # published full-population rates. A 0.5-per-1,000 tolerance catches gross
  # errors (e.g. a misaligned fixed-width field) while tolerating that known
  # downward bias.
  tol <- 0.5
  check("CDC join yields both race columns and >= 20 overlapping years", cols_ok)
  check("Black IMR within 0.5/1,000 of CDC for all overlapping years",
        cols_ok && all(abs(cmp$black - cmp$cdc_imr_black) < tol))
  check("White IMR within 0.5/1,000 of CDC for all overlapping years",
        cols_ok && all(abs(cmp$white - cmp$cdc_imr_white) < tol))
  # Exclusions remove records, so pipeline rates should not sit above CDC on
  # average (checked for both races; individual years such as 2014 may exceed).
  check("pipeline IMRs not systematically above CDC (Black and White)",
        cols_ok &&
          mean(cmp$black - cmp$cdc_imr_black) <= 0 &&
          mean(cmp$white - cmp$cdc_imr_white) <= 0)
}

cat("\n")

#### Summary ####
cat("==============================================\n")
cat("Results:", pass, "passed,", fail, "failed\n")
cat("==============================================\n")

if (fail > 0) {
  stop(fail, " test(s) failed")
} else {
  cat("All tests passed.\n")
}
