# Categorize ICD-10 codes into cause-of-death groups
# Categories follow Wolf et al. (2025, JAMA Pediatrics) eTable 1
#
# This function is shared across:
#   - scripts/02-analyze_data.R
#   - scripts/03-sensitivity_analyses.R
#   - scripts/04-test_cause_categories.R
#
# Edit this file to update cause categories, then re-run
# 04-test_cause_categories.R to verify correctness.

categorize_cause <- function(ucod) {
  ucod_3 <- str_sub(ucod, 1, 3)
  first_char <- str_sub(ucod, 1, 1)
  case_when(
    first_char == "Q" ~ "Congenital malformations",
    ucod_3 == "P07" ~ "Short gestation/LBW",
    # SUID: SIDS + ill-defined/unspecified + accidental suffocation in bed
    ucod_3 == "R95" | ucod_3 == "R99" | ucod_3 == "W75" ~
      "Sudden unexpected infant death",
    # Unintentional injuries: V01-V99, W00-W74, W76-W99, X00-X59
    first_char == "V" ~ "Unintentional injuries",
    first_char == "W" & ucod_3 != "W75" ~ "Unintentional injuries",
    first_char == "X" & ucod_3 <= "X59" ~ "Unintentional injuries",
    ucod_3 == "P01" ~ "Maternal complications",
    ucod_3 == "P02" ~ "Placenta/cord/membranes",
    # Infection: A00-B99, P35-P39 (includes bacterial sepsis P36)
    first_char == "A" | first_char == "B" |
      (ucod_3 >= "P35" & ucod_3 <= "P39") ~ "Infection",
    ucod_3 == "P22" ~ "Respiratory distress",
    first_char == "I" ~ "Circulatory diseases",
    (ucod_3 >= "P50" & ucod_3 <= "P52") | ucod_3 == "P54" ~ "Hemorrhage",
    ucod_3 >= "P20" & ucod_3 <= "P21" ~ "Hypoxia and birth asphyxia",
    ucod_3 == "P77" ~ "Necrotizing enterocolitis",
    # Assault (homicide): U01, X85-Y09
    ucod_3 == "U01" ~ "Assault (homicide)",
    first_char == "X" & ucod_3 >= "X85" ~ "Assault (homicide)",
    first_char == "Y" & ucod_3 <= "Y09" ~ "Assault (homicide)",
    TRUE ~ "Other"
  )
}

# Same as categorize_cause(), but splits the combined "Sudden unexpected infant
# death" category into its three ICD-10 components: SIDS (R95), accidental
# suffocation/strangulation in bed (W75), and ill-defined/unspecified (R99).
categorize_cause_suid_split <- function(ucod) {
  ucod_3 <- str_sub(ucod, 1, 3)
  case_when(
    ucod_3 == "R95" ~ "SIDS (R95)",
    ucod_3 == "W75" ~ "Accidental suffocation in bed (W75)",
    ucod_3 == "R99" ~ "Ill-defined/unspecified (R99)",
    TRUE ~ categorize_cause(ucod)
  )
}
