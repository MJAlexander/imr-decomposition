# Birthweight / gestational-age plausibility helpers
#
# Shared by:
#   - scripts/01-clean_data.R (applies the exclusions during cleaning)
#   - scripts/04-test_cause_categories.R (unit tests)
#
# Implements the inclusion criteria of Alexander GR, Himes JH, Kaufman RB,
# Mor J, Kogan M. A United States National Reference for Fetal Growth.
# Obstet Gynecol 1996;87:163-168. Table 1 gives inclusive birthweight ranges
# by gestational age week.
#
# Requires tidyverse to be attached by the sourcing script (case_when is used
# unqualified, matching the scripts/categorize_cause.R convention).

#### Alexander 1996 ####

# Lower bound of plausible birthweight (g) by GA week. The Alexander 1996
# Table 1 rows cover 20-47 weeks; we extrapolate the boundary rows so no record
# is excluded purely for an out-of-range gestational age: weeks < 20 use the
# 20-week bound and weeks >= 47 use the >=38-week bound. Returns NA only when
# weeks itself is NA.
alexander_bwt_lo <- function(weeks) {
  case_when(
    weeks <= 21 ~ 125L,              # 20-week bound, extrapolated to weeks < 20
    weeks == 22 ~ 125L,
    weeks == 23 ~ 125L,
    weeks == 24 ~ 125L,
    weeks == 25 ~ 250L,
    weeks == 26 ~ 250L,
    weeks == 27 ~ 250L,
    weeks == 28 ~ 250L,
    weeks == 29 ~ 250L,
    weeks == 30 ~ 375L,
    weeks == 31 ~ 375L,
    weeks == 32 ~ 500L,
    weeks == 33 ~ 500L,
    weeks == 34 ~ 750L,
    weeks == 35 ~ 750L,
    weeks == 36 ~ 750L,
    weeks == 37 ~ 1000L,
    weeks >= 38 ~ 1000L,             # >=38-week bound, extrapolated to >= 47
    TRUE ~ NA_integer_
  )
}

# Upper bound of plausible birthweight (g) by GA week. Boundary rows are
# extrapolated as in alexander_bwt_lo. Returns NA only when weeks is NA.
alexander_bwt_hi <- function(weeks) {
  case_when(
    weeks <= 21 ~ 1250L,             # 20-week bound, extrapolated to weeks < 20
    weeks == 22 ~ 1375L,
    weeks == 23 ~ 1500L,
    weeks == 24 ~ 1625L,
    weeks == 25 ~ 1750L,
    weeks == 26 ~ 2000L,
    weeks == 27 ~ 2250L,
    weeks == 28 ~ 2500L,
    weeks == 29 ~ 2750L,
    weeks == 30 ~ 3000L,
    weeks == 31 ~ 3250L,
    weeks == 32 ~ 3500L,
    weeks == 33 ~ 3750L,
    weeks == 34 ~ 4000L,
    weeks == 35 ~ 4500L,
    weeks == 36 ~ 5000L,
    weeks == 37 ~ 5500L,
    weeks >= 38 ~ 6000L,             # >=38-week bound, extrapolated to >= 47
    TRUE ~ NA_integer_
  )
}

# Single source of truth for the Alexander classification. Returns one of the
# levels in alexander_reason_levels. The cleaning pipeline derives both the
# boolean exclusion flag (reason != "plausible") and the exclusion-count
# partition (one count per reason) from this one factor, so they cannot disagree.
alexander_reason_levels <- c("plausible", "missing_bwt", "ga_out_of_range",
                             "implausibly_low", "implausibly_high")

alexander_reason <- function(bwt, weeks) {
  lo <- alexander_bwt_lo(weeks)
  hi <- alexander_bwt_hi(weeks)
  case_when(
    is.na(bwt) ~ "missing_bwt",
    # Bounds are defined for every non-NA gestational age (boundary rows
    # extrapolated), so this only fires if weeks itself is NA. Such records are
    # already excluded earlier as unknown GA, so this level is effectively a
    # defensive guard and contributes 0 to the analysed exclusion counts.
    is.na(lo) ~ "ga_out_of_range",
    bwt < lo ~ "implausibly_low",
    bwt > hi ~ "implausibly_high",
    TRUE ~ "plausible"
  )
}
