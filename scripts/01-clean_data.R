#### Preamble ####
# Purpose: Process raw NCHS linked birth-infant death files (2003-2024)
#          to create analysis dataset for Kitagawa decomposition
# Author: Monica Alexander
# Date: 2026-01-31
# Contact: monica.alexander@utoronto.ca
# License: MIT
# Pre-requisites:
#   - Download NCHS Period Linked Birth-Infant Death files from:
#     https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm
#   - For 2003-2017: Place LinkPE03US.zip through LinkPE17US.zip in:
#     data/raw_data/linked_birth_death/
#   - For 2018-2024: Place period-cohort-linked files (e.g., 2024PE2023CO.zip) in:
#     data/raw_data/period_cohort_linked/
#   - Required R packages: tidyverse, archive
# Notes:
#   - Uses period linked files which contain both birth counts and infant deaths
#   - Analysis begins in 2003 (revised birth certificate). Pre-2003 LinkPE files
#     use distinct layouts that are not supported and are filtered out.
#   - Record layouts differ across years:
#     * 2003: 751-byte records (transition year, unique mother's age position)
#     * 2004-2013: 751-byte denominator, 1142-byte numerator records
#     * 2014-2024: 1346-byte denominator, 1743-byte numerator records
#   - Gestational age methodology changed in 2014 (LMP to OE estimate)
#   - Mother's age available for stratification (added 2024)

#### Workspace setup ####
library(tidyverse)
library(archive)

# Create output directories
dir.create("data/analysis_data", recursive = TRUE, showWarnings = FALSE)

#### Define file paths ####
# Period linked files (2003-2017). The analysis begins in 2003 (the revised
# birth certificate); pre-2003 LinkPE files use distinct fixed-width layouts
# that are not supported here and would be dropped downstream regardless, so
# they are filtered out at discovery rather than read.
linked_dir <- "data/raw_data/linked_birth_death"
linked_files <- list.files(linked_dir, pattern = "LinkPE\\d{2}US\\.zip$", full.names = TRUE)
linkpe_year <- function(f) {
  s <- as.integer(sub(".*LinkPE(\\d{2})US.*", "\\1", basename(f)))
  ifelse(s >= 95, 1900L + s, 2000L + s)
}
linked_files <- linked_files[linkpe_year(linked_files) >= 2003]

# Period-cohort linked files (2018-2024)
period_cohort_dir <- "data/raw_data/period_cohort_linked"
period_cohort_files <- list.files(period_cohort_dir, pattern = "\\d{4}PE\\d{4}CO\\.zip$", full.names = TRUE)

cat("Found", length(linked_files), "period linked files (2003-2017)\n")
cat("Found", length(period_cohort_files), "period-cohort linked files (2018-2024)\n")

#### Define record layouts by era ####

# Era definitions based on NCHS documentation
# Key variables needed:
# - Year of birth
# - Mother's race (need to identify Non-Hispanic Black and Non-Hispanic White)
# - Gestational age (weeks)
# - For numerator: death indicator or count

# 2014-2024 denominator file layout (based on 2014+ User Guides)
# Record length: 1346 bytes for denominator, 1743 bytes for numerator
# Format is consistent from 2014 through 2024 (2003 revised birth certificate)
# Verified positions through file inspection
layout_2014_2024_den <- fwf_cols(
  year = c(9, 12),               # Year of birth (positions 9-12)
  mrace = c(117, 117),           # Mother's single race (1=White, 2=Black, etc.)
  mhisp = c(115, 115),           # Mother's Hispanic origin (0=Non-Hispanic, 1-5=Hispanic)
  mage = c(75, 76),              # Mother's age in single years (MAGER)
  gest_weeks = c(490, 491),      # OE gestational age in weeks (17-47, 99=unknown)
  bwt = c(512, 515)              # Birth weight in grams (DBWT), verified by inspection
)

# 2014-2024 numerator file layout (same positions as denominator for key fields)
layout_2014_2024_num <- fwf_cols(
  year = c(9, 12),               # Year of birth
  mrace = c(117, 117),           # Mother's single race
  mhisp = c(115, 115),           # Mother's Hispanic origin
  mage = c(75, 76),              # Mother's age in single years (MAGER)
  gest_weeks = c(490, 491),      # OE gestational age in weeks
  bwt = c(512, 515),             # Birth weight in grams (DBWT)
  ucod = c(1368, 1371)           # Underlying cause of death (ICD-10, 4 chars)
)

# 2004-2013 file layout (based on 2003-2013 User Guide)
# Record length: ~751 bytes for denominator, ~1142 bytes for numerator
layout_2004_2013_den <- fwf_cols(
  year = c(15, 18),              # Year of birth
  mrace = c(143, 143),           # Mother's race (bridged, 1=White, 2=Black)
  mhisp = c(148, 148),           # Mother's Hispanic origin (0=Non-Hispanic, 1-5=Hispanic)
  mage = c(89, 90),              # Mother's age in years (verified position)
  gest_weeks = c(451, 452),      # Combined gestational age (LMP-based)
  bwt = c(467, 470)              # Birth weight in grams (DBWT), verified by inspection
)

layout_2004_2013_num <- fwf_cols(
  year = c(15, 18),              # Year of birth
  mrace = c(143, 143),           # Mother's race (bridged)
  mhisp = c(148, 148),           # Mother's Hispanic origin
  mage = c(89, 90),              # Mother's age in years
  gest_weeks = c(451, 452),      # Combined gestational age (LMP-based)
  bwt = c(467, 470),             # Birth weight in grams (DBWT)
  ucod = c(884, 887)             # Underlying cause of death (ICD-10, 4 chars)
)

# 2003 file layout (transition year - unique layout for mother's age)
# The 2003 revised birth certificate was introduced but mother's age position differs
# NOTE: Position 184-185 is the best available but the distribution is imperfect
#   (Under 20 ~4% vs expected ~10%, 40+ ~8% vs expected ~3%). This is a known
#   data quality limitation for the 2003 transition year.
layout_2003_den <- fwf_cols(
  year = c(15, 18),              # Year of birth
  mrace = c(143, 143),           # Mother's race (bridged, 1=White, 2=Black)
  mhisp = c(148, 148),           # Mother's Hispanic origin (0=Non-Hispanic, 1-5=Hispanic)
  mage = c(184, 185),            # Mother's age in years (2003-specific position)
  gest_weeks = c(451, 452),      # Combined gestational age (LMP-based)
  bwt = c(467, 470)              # Birth weight in grams (DBWT), verified by inspection
)

layout_2003_num <- fwf_cols(
  year = c(15, 18),              # Year of birth
  mrace = c(143, 143),           # Mother's race (bridged)
  mhisp = c(148, 148),           # Mother's Hispanic origin
  mage = c(184, 185),            # Mother's age in years (2003-specific position)
  gest_weeks = c(451, 452),      # Combined gestational age (LMP-based)
  bwt = c(467, 470),             # Birth weight in grams (DBWT)
  ucod = c(767, 770)             # Underlying cause of death (ICD-10, 4 chars)
)

#### Helper functions ####

# Determine which era a year belongs to (analysis starts in 2003; pre-2003
# files are filtered out at discovery and never reach here).
get_era <- function(year) {
  if (year >= 2014) return("2014_2024")
  if (year >= 2004) return("2004_2013")
  return("2003")
}

# Get the correct layout for a given year and file type
get_layout <- function(year, file_type = "den") {
  era <- get_era(year)
  if (era == "2014_2024") {
    if (file_type == "den") return(layout_2014_2024_den)
    return(layout_2014_2024_num)
  } else if (era == "2004_2013") {
    if (file_type == "den") return(layout_2004_2013_den)
    return(layout_2004_2013_num)
  } else {
    if (file_type == "den") return(layout_2003_den)
    return(layout_2003_num)
  }
}

# Standardize race coding across eras
# Goal: Extract *confirmed* Non-Hispanic White and Non-Hispanic Black only.
# A record is kept only when Hispanic origin is explicitly Non-Hispanic (0) and
# race is White (1) or Black (2). Records with Hispanic origin (1-5), unknown or
# missing Hispanic origin (9 or NA), or any other race code are excluded (NA).
# Unknown origin is excluded rather than assumed Non-Hispanic, so the analysis
# population does not include records whose ethnicity is not actually known.
# Race and Hispanic coding is consistent across all eras:
#   mrace: 1=White, 2=Black, others=excluded
#   mhisp: 0=Non-Hispanic, 1-5=Hispanic, 9=Unknown/Not stated
standardize_race <- function(mrace, mhisp) {
  case_when(
    mhisp == 0 & mrace == 1 ~ "Non-Hispanic White",
    mhisp == 0 & mrace == 2 ~ "Non-Hispanic Black",
    TRUE ~ NA_character_
  )
}

# Categorize gestational age into the 7 main-analysis groups
# Note: NCHS uses 99 for unknown gestational age - must exclude before other checks
# (5- and 16-group schemes are run as sensitivity analyses in 03-sensitivity_analyses.R)
categorize_gest_age <- function(weeks) {
  case_when(
    is.na(weeks) | weeks == 99 ~ NA_character_,
    weeks < 28 ~ "<28 weeks",
    weeks >= 28 & weeks <= 31 ~ "28-31 weeks",
    weeks >= 32 & weeks <= 33 ~ "32-33 weeks",
    weeks >= 34 & weeks <= 36 ~ "34-36 weeks",
    weeks >= 37 & weeks <= 38 ~ "37-38 weeks",
    weeks >= 39 & weeks <= 41 ~ "39-41 weeks",
    weeks >= 42 ~ "42+ weeks",
    TRUE ~ NA_character_
  )
}

# Birthweight / gestational-age plausibility helpers (Alexander 1996) live in a
# shared module so the cleaning script and the test suite use the same
# definitions (mirrors scripts/categorize_cause.R).
source("scripts/birthweight_helpers.R")

# Categorize mother's age into groups
# Note: NCHS uses 99 for unknown mother's age - must exclude before other checks
categorize_mother_age <- function(age) {
  case_when(
    age == 99 ~ NA_character_,
    age < 20 ~ "Under 20",
    age >= 20 & age <= 24 ~ "20-24",
    age >= 25 & age <= 29 ~ "25-29",
    age >= 30 & age <= 34 ~ "30-34",
    age >= 35 & age <= 39 ~ "35-39",
    age >= 40 ~ "40 and over",
    TRUE ~ NA_character_
  )
}

# Build a one-row exclusion-counts tibble for either a births_all or
# deaths_all data frame produced by process_linked_file. Counts records by
# sequential exclusion: race -> gestational age -> Alexander 1996 birthweight
# implausibility. The Alexander reason breakdown is read from the per-record
# bwt_ga_reason column (set once in the mutate via alexander_reason) so the
# counts cannot drift from the flag.
#
# Columns are consumed by eTable 1 (totals + Alexander reason breakdown).
# Counts that are exact arithmetic of other stored columns are not kept.
build_exclusion_tibble <- function(df) {
  ga_known <- !is.na(df$race) & !is.na(df$gestational_age_group)
  reason <- df$bwt_ga_reason
  tibble(
    year = df$year[1],
    total_records = nrow(df),
    excluded_race = sum(is.na(df$race)),
    excluded_race_unknown = sum(
      is.na(df$race) &
      !(df$mhisp %in% c(1, 2, 3, 4, 5)) &
      (is.na(df$mrace) | !(df$mrace %in% 1:8))
    ),
    # White/Black-race mothers excluded solely because Hispanic origin is
    # unknown or missing (mhisp is 9 or NA, not the explicit non-Hispanic 0).
    # These are a subset of excluded_race; reported in the Methods.
    excluded_hispanic_unknown = sum(
      (df$mrace %in% c(1, 2)) &
      !(df$mhisp %in% c(1, 2, 3, 4, 5)) &
      (is.na(df$mhisp) | df$mhisp != 0)
    ),
    excluded_ga = sum(!is.na(df$race) & is.na(df$gestational_age_group)),
    excluded_implausible_bwt_ga = sum(ga_known & df$implausible_bwt_ga),
    # Mutually exclusive partition of excluded_implausible_bwt_ga, read from
    # the single-source alexander_reason() factor. Gestational age is never an
    # exclusion reason on its own (boundary bounds are extrapolated), so there
    # is no ga_out_of_range count: missing + low + high == implausible total.
    excluded_bwt_missing = sum(ga_known & reason == "missing_bwt"),
    excluded_bwt_low = sum(ga_known & reason == "implausibly_low"),
    excluded_bwt_high = sum(ga_known & reason == "implausibly_high"),
    included = sum(ga_known & !df$implausible_bwt_ga)
  )
}

# Extract year from original LinkPE filename (2003-2017)
extract_year_from_linkpe <- function(filename) {
  match <- str_extract(basename(filename), "(?<=LinkPE)\\d{2}")
  if (is.na(match)) return(NA)
  year_suffix <- as.integer(match)
  if (year_suffix >= 95) {
    return(1900 + year_suffix)
  } else {
    return(2000 + year_suffix)
  }
}

# Extract year from period-cohort filename (2018-2024)
# Files named like "2024PE2023CO.zip" - we want the period year (first 4 digits)
extract_year_from_period_cohort <- function(filename) {
  match <- str_extract(basename(filename), "^\\d{4}")
  if (is.na(match)) return(NA)
  return(as.integer(match))
}

# Find the denominator and numerator files within a zip archive
# For period-cohort files, we need to select the PERIOD files (not cohort)
# Period files have the period year in the filename (e.g., VS20LINK for 2020, VS2021LINK for 2021)
find_files_in_archive <- function(zip_path, period_year = NULL) {
  files <- archive(zip_path)
  den_files <- files$path[str_detect(str_to_upper(files$path), "DEN")]
  num_files <- files$path[str_detect(str_to_upper(files$path), "NUM")]

 # If period_year is specified, filter for files matching that year
  if (!is.null(period_year)) {
    # Create patterns to match the period year (e.g., "VS20" for 2020, "VS2021" for 2021)
    year_2digit <- substr(as.character(period_year), 3, 4)
    year_4digit <- as.character(period_year)
    # Match either VS{2digit}LINK or VS{4digit}LINK at the start
    pattern <- paste0("^VS(", year_2digit, "|", year_4digit, ")LINK")

    den_match <- str_detect(den_files, regex(pattern, ignore_case = TRUE))
    num_match <- str_detect(num_files, regex(pattern, ignore_case = TRUE))

    if (any(den_match)) den_files <- den_files[den_match]
    if (any(num_match)) num_files <- num_files[num_match]
  }

  list(den = den_files, num = num_files)
}

#### Process each year's linked file ####

# Generic processing function that handles both file types
process_linked_file <- function(zip_path, file_type = "linkpe") {
  # Extract year based on file type
  if (file_type == "linkpe") {
    file_year <- extract_year_from_linkpe(zip_path)
  } else {
    file_year <- extract_year_from_period_cohort(zip_path)
  }

  if (is.na(file_year)) {
    warning("Could not extract year from: ", zip_path)
    return(NULL)
  }

  cat("Processing year", file_year, "...\n")

  # Find files in archive (pass year for period-cohort files to select correct period data)
  if (file_type == "period_cohort") {
    archive_files <- find_files_in_archive(zip_path, period_year = file_year)
  } else {
    archive_files <- find_files_in_archive(zip_path)
  }

  if (length(archive_files$den) == 0 || length(archive_files$num) == 0) {
    warning("Missing denominator or numerator file in: ", zip_path)
    return(NULL)
  }

  # Get appropriate layouts based on file year
  den_layout <- get_layout(file_year, "den")
  num_layout <- get_layout(file_year, "num")

  # Read denominator file (all births)
  den_data <- tryCatch({
    read_fwf(
      archive_read(zip_path, file = archive_files$den[1]),
      col_positions = den_layout,
      col_types = cols(.default = "c"),
      progress = FALSE
    )
  }, error = function(e) {
    warning("Error reading denominator file for year ", file_year, ": ", e$message)
    return(NULL)
  })

  if (is.null(den_data)) return(NULL)

  # Read numerator file (infant deaths)
  num_data <- tryCatch({
    read_fwf(
      archive_read(zip_path, file = archive_files$num[1]),
      col_positions = num_layout,
      col_types = cols(.default = "c"),
      progress = FALSE
    )
  }, error = function(e) {
    warning("Error reading numerator file for year ", file_year, ": ", e$message)
    return(NULL)
  })

  if (is.null(num_data)) return(NULL)

  # Handle the 2014 era transition in the numerator file
  # The 2014 numerator contains 2013-born death records that were reformatted to
  # the 2014 layout. All fields use the same positions EXCEPT mrace: 2013-born
  # records store bridged race at position 110 instead of detailed race at 117.
  # We re-read mrace from position 110 for those records.
  # (The 2003 file's 2002-born records use the same 2003 layout and need no fix.)
  if (file_year == 2014) {
    cat("  Note: Fixing mrace for 2013-born death records (position 110)\n")
    # Re-read just the bridged race field at position 110
    bridged_race_data <- tryCatch({
      read_fwf(
        archive_read(zip_path, file = archive_files$num[1]),
        col_positions = fwf_cols(birth_year = c(9, 12), mrace_bridged = c(110, 110)),
        col_types = cols(.default = "c"),
        progress = FALSE
      )
    }, error = function(e) {
      warning("Error re-reading bridged race for 2014: ", e$message)
      return(NULL)
    })

    if (!is.null(bridged_race_data)) {
      # For 2013-born records, replace mrace with the bridged race from position 110
      is_prior_year <- as.integer(bridged_race_data$birth_year) == 2013
      num_data$mrace[is_prior_year] <- bridged_race_data$mrace_bridged[is_prior_year]
      cat("    Corrected mrace for", sum(is_prior_year), "prior-year records\n")
    }
  }

  # All processed files are 2003+ (pre-2003 are filtered out at discovery), so
  # the Alexander birthweight check always applies. This flag is a defensive
  # guard should the file filter ever change.
  apply_alexander <- file_year >= 2003

  # Process denominator (births) - include mother's age and birth weight
  births_all <- den_data |>
    mutate(
      year = as.integer(year),
      mrace = as.integer(mrace),
      mhisp = as.integer(mhisp),
      mage = as.integer(mage),
      gest_weeks = as.integer(gest_weeks),
      # NCHS codes "not stated" birthweight as 9999; treat as missing.
      bwt = if ("bwt" %in% colnames(den_data)) na_if(as.integer(bwt), 9999L) else NA_integer_,
      race = standardize_race(mrace, mhisp),
      gestational_age_group = categorize_gest_age(gest_weeks),
      mother_age_group = categorize_mother_age(mage),
      # Single-source Alexander classification: one factor per record, from
      # which both the boolean flag and the exclusion-reason partition derive.
      bwt_ga_reason = if (apply_alexander)
        alexander_reason(bwt, gest_weeks) else "plausible",
      implausible_bwt_ga = bwt_ga_reason != "plausible"
    )

  # Process numerator (deaths) - include mother's age and birth weight
  # NOTE: The year field in death records is the BIRTH year of the infant, not the
  # death year. Period linked files contain all infant deaths occurring in the period
  # year, including infants born in the prior year (~12-13% of deaths). We override
  # year to file_year so all period deaths are attributed to the correct period,
  # matching the CDC's period infant mortality rate definition.
  deaths_all <- num_data |>
    mutate(
      year = file_year,
      mrace = as.integer(mrace),
      mhisp = as.integer(mhisp),
      mage = as.integer(mage),
      gest_weeks = as.integer(gest_weeks),
      # NCHS codes "not stated" birthweight as 9999; treat as missing.
      bwt = if ("bwt" %in% colnames(num_data)) na_if(as.integer(bwt), 9999L) else NA_integer_,
      race = standardize_race(mrace, mhisp),
      gestational_age_group = categorize_gest_age(gest_weeks),
      mother_age_group = categorize_mother_age(mage),
      # Single-source Alexander classification (see births block above).
      bwt_ga_reason = if (apply_alexander)
        alexander_reason(bwt, gest_weeks) else "plausible",
      implausible_bwt_ga = bwt_ga_reason != "plausible",
      ucod = str_trim(ucod)
    )

  # Count all exclusions before filtering. Exclusion categories are applied
  # sequentially: race -> gestational age -> Alexander 1996 implausibility.
  exclusion_births <- build_exclusion_tibble(births_all)
  exclusion_deaths <- build_exclusion_tibble(deaths_all)

  # Filter to known GA and plausible birthweight for main analysis
  births <- births_all |>
    filter(!is.na(race), !is.na(gestational_age_group), !implausible_bwt_ga) |>
    count(year, race, gestational_age_group, mother_age_group, name = "births")

  deaths <- deaths_all |>
    filter(!is.na(race), !is.na(gestational_age_group), !implausible_bwt_ga) |>
    count(year, race, gestational_age_group, mother_age_group, name = "deaths")

  # Aggregate deaths by cause AND gestational age (for cause-by-GA appendix)
  # Saves raw ICD-10 codes; categorization done in analysis script
  deaths_by_cause_ga <- deaths_all |>
    filter(!is.na(race), !is.na(gestational_age_group), !implausible_bwt_ga) |>
    count(year, race, gestational_age_group, ucod, name = "deaths")

  # Per-single-week aggregates (for the 5- and 16-category GA sensitivity
  # analyses in 03-sensitivity_analyses.R; the 7-category scheme is the main
  # analysis and uses analysis_data.csv directly)
  births_by_week <- births_all |>
    filter(!is.na(race), !is.na(gestational_age_group), !implausible_bwt_ga) |>
    count(year, race, gest_weeks, name = "births")

  deaths_by_week <- deaths_all |>
    filter(!is.na(race), !is.na(gestational_age_group), !implausible_bwt_ga) |>
    count(year, race, gest_weeks, name = "deaths")

  data_by_week <- births_by_week |>
    left_join(deaths_by_week, by = c("year", "race", "gest_weeks")) |>
    mutate(deaths = replace_na(deaths, 0))

  # Combine births and deaths
  result <- births |>
    left_join(deaths, by = c("year", "race", "gestational_age_group", "mother_age_group")) |>
    mutate(deaths = replace_na(deaths, 0))

  cat("  Births:", format(sum(result$births), big.mark = ","), "\n")
  cat("  Deaths:", format(sum(result$deaths), big.mark = ","), "\n")

  return(list(data = result,
              exclusions_births = exclusion_births,
              exclusions_deaths = exclusion_deaths,
              deaths_by_cause_ga = deaths_by_cause_ga,
              data_by_week = data_by_week))
}

#### Main processing loop ####

cat("\n=== Processing NCHS Linked Birth-Infant Death Files ===\n\n")

# Process all files
all_results <- list()
all_exclusions_births <- list()
all_exclusions_deaths <- list()
all_deaths_by_cause_ga <- list()
all_data_by_week <- list()

# Process original LinkPE files (2003-2017)
if (length(linked_files) > 0) {
  cat("Processing period linked files (2003-2017)...\n")
  for (zip_file in linked_files) {
    result <- process_linked_file(zip_file, file_type = "linkpe")
    if (!is.null(result)) {
      all_results[[length(all_results) + 1]] <- result$data
      all_exclusions_births[[length(all_exclusions_births) + 1]] <- result$exclusions_births
      all_exclusions_deaths[[length(all_exclusions_deaths) + 1]] <- result$exclusions_deaths
      all_deaths_by_cause_ga[[length(all_deaths_by_cause_ga) + 1]] <- result$deaths_by_cause_ga
      all_data_by_week[[length(all_data_by_week) + 1]] <- result$data_by_week
    }
  }
}

# Process period-cohort linked files (2018-2024)
if (length(period_cohort_files) > 0) {
  cat("\nProcessing period-cohort linked files (2018-2024)...\n")
  for (zip_file in period_cohort_files) {
    result <- process_linked_file(zip_file, file_type = "period_cohort")
    if (!is.null(result)) {
      all_results[[length(all_results) + 1]] <- result$data
      all_exclusions_births[[length(all_exclusions_births) + 1]] <- result$exclusions_births
      all_exclusions_deaths[[length(all_exclusions_deaths) + 1]] <- result$exclusions_deaths
      all_deaths_by_cause_ga[[length(all_deaths_by_cause_ga) + 1]] <- result$deaths_by_cause_ga
      all_data_by_week[[length(all_data_by_week) + 1]] <- result$data_by_week
    }
  }
}

# Combine all years
if (length(all_results) == 0) {
  stop("No data could be processed. Check that linked birth-death files are present.")
}

analysis_data_raw <- bind_rows(all_results)

# Combine and save exclusion counts
exclusion_counts <- bind_rows(
  bind_rows(all_exclusions_births) |> mutate(file_type = "births"),
  bind_rows(all_exclusions_deaths) |> mutate(file_type = "deaths")
) |>
  select(year, file_type, total_records, excluded_race, excluded_race_unknown,
         excluded_hispanic_unknown, excluded_ga, excluded_implausible_bwt_ga,
         excluded_bwt_missing, excluded_bwt_low, excluded_bwt_high, included)
write_csv(exclusion_counts, "data/analysis_data/exclusion_counts.csv")
cat("Exclusion counts saved: data/analysis_data/exclusion_counts.csv\n")

# Combine deaths by cause and gestational age
cause_ga_data <- bind_rows(all_deaths_by_cause_ga) |>
  summarise(deaths = sum(deaths), .by = c(year, race, gestational_age_group, ucod))

# Get births per year-race-GA for rate computation
births_by_year_race_ga <- analysis_data_raw |>
  summarise(births = sum(births), .by = c(year, race, gestational_age_group))

cause_ga_data <- cause_ga_data |>
  left_join(births_by_year_race_ga, by = c("year", "race", "gestational_age_group"))

write_csv(cause_ga_data, "data/analysis_data/cause_of_death_by_ga.csv")
cat("Cause-by-GA data saved: data/analysis_data/cause_of_death_by_ga.csv\n")

# Per-single-week aggregates (for the 5- and 16-category GA sensitivity
# analyses; the 7-category scheme is the main analysis)
data_by_week <- bind_rows(all_data_by_week) |>
  summarise(births = sum(births), deaths = sum(deaths),
            .by = c(year, race, gest_weeks)) |>
  arrange(year, race, gest_weeks)
write_csv(data_by_week, "data/analysis_data/analysis_data_by_week.csv")
cat("Per-week analysis data saved: data/analysis_data/analysis_data_by_week.csv\n")

cat("\n=== Processing Complete ===\n")
cat("Years processed:", n_distinct(analysis_data_raw$year), "\n")
cat("Year range:", min(analysis_data_raw$year), "-", max(analysis_data_raw$year), "\n")

#### Aggregate and calculate rates ####

# Order gestational age groups (7 main-analysis groups)
ga_levels <- c(
  "<28 weeks",
  "28-31 weeks",
  "32-33 weeks",
  "34-36 weeks",
  "37-38 weeks",
  "39-41 weeks",
  "42+ weeks"
)

# Order mother's age groups
mage_levels <- c(
  "Under 20",
  "20-24",
  "25-29",
  "30-34",
  "35-39",
  "40 and over"
)

# Save detailed data with mother's age for stratified analysis
analysis_data_by_age <- analysis_data_raw |>
  filter(!is.na(mother_age_group)) |>
  mutate(
    gestational_age_group = factor(gestational_age_group, levels = ga_levels),
    mother_age_group = factor(mother_age_group, levels = mage_levels)
  ) |>
  mutate(
    mortality_rate = deaths / births * 1000,
    .by = c(year, race, gestational_age_group, mother_age_group)
  ) |>
  arrange(year, race, gestational_age_group, mother_age_group)

# Main analysis: aggregate over mother's age
analysis_data <- analysis_data_raw |>
  summarise(
    births = sum(births),
    deaths = sum(deaths),
    .by = c(year, race, gestational_age_group)
  ) |>
  mutate(
    gestational_age_group = factor(gestational_age_group, levels = ga_levels)
  ) |>
  # Calculate mortality rate and birth proportion
  mutate(
    mortality_rate = deaths / births * 1000,
    birth_proportion = births / sum(births),
    .by = c(year, race)
  ) |>
  arrange(year, race, gestational_age_group)

#### Validation ####

# Check that all year-race combinations have all 7 gestational age groups
complete_check <- analysis_data |>
  summarise(n_groups = n_distinct(gestational_age_group), .by = c(year, race))

incomplete <- complete_check |>
  filter(n_groups < 7)

if (nrow(incomplete) > 0) {
  cat("\nWarning: Some year-race combinations have incomplete gestational age data:\n")
  print(incomplete)

  # Keep only complete years
  complete_years <- complete_check |>
    filter(n_groups == 7) |>
    summarise(n = n(), .by = year) |>
    filter(n == 2) |>
    pull(year)

  analysis_data <- analysis_data |>
    filter(year %in% complete_years)

  cat("\nKeeping only complete years:", paste(complete_years, collapse = ", "), "\n")
}

# Check proportions sum to 1
prop_check <- analysis_data |>
  summarise(total = sum(birth_proportion), .by = c(year, race))

if (any(abs(prop_check$total - 1) > 0.001)) {
  warning("Birth proportions do not sum to 1 for some year-race combinations")
}

#### Summary statistics ####

cat("\n=== Data Summary ===\n")
cat("Final years:", min(analysis_data$year), "-", max(analysis_data$year), "\n")
cat("Race groups:", paste(unique(analysis_data$race), collapse = ", "), "\n")
cat("Gestational age groups:", length(levels(analysis_data$gestational_age_group)), "\n")
cat("Total observations:", nrow(analysis_data), "\n\n")

# Summary by race
summary_by_race <- analysis_data |>
  summarise(
    total_births = sum(births),
    total_deaths = sum(deaths),
    imr = sum(deaths) / sum(births) * 1000,
    .by = race
  )

cat("Overall summary by race:\n")
print(summary_by_race)

# IMR by year and race
cat("\nInfant Mortality Rate by year (first and last years):\n")
imr_by_year <- analysis_data |>
  summarise(
    total_births = sum(births),
    total_deaths = sum(deaths),
    imr = sum(deaths) / sum(births) * 1000,
    .by = c(year, race)
  ) |>
  filter(year %in% c(min(year), max(year)))

print(imr_by_year)

#### Save cleaned data ####

# Save main analysis data (aggregated over mother's age)
write_csv(analysis_data, "data/analysis_data/analysis_data.csv")

# Save data stratified by mother's age for additional analyses
write_csv(analysis_data_by_age, "data/analysis_data/analysis_data_by_mother_age.csv")

cat("\n=== Output Files ===\n")
cat("Main analysis data: data/analysis_data/analysis_data.csv\n")
cat("By mother's age: data/analysis_data/analysis_data_by_mother_age.csv\n")
cat("Per-week analysis data: data/analysis_data/analysis_data_by_week.csv\n")
cat("Cause-by-GA data: data/analysis_data/cause_of_death_by_ga.csv\n")
cat("Exclusion counts: data/analysis_data/exclusion_counts.csv\n")
