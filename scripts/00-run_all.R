#### Preamble ####
# Purpose: Master script to run all analysis scripts in sequence
# Author: Monica Alexander
# Date: 2026-01-31
# Contact: monica.alexander@utoronto.ca
# License: MIT
# Pre-requisites: R packages listed in README
# Note: Run this script from the project root directory

#### Setup ####
# Record start time
start_time <- Sys.time()
cat("==============================================\n")
cat("Running full analysis pipeline\n")
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================\n\n")

# Set working directory to project root (if not already there)
if (!file.exists("decomposition.Rproj")) {
  stop("Please run this script from the project root directory")
}

#### Run scripts in sequence ####

cat("Step 1: Data cleaning\n")
cat("-----------------------------------\n")
source("scripts/01-clean_data.R")
cat("\n")

cat("Step 2: Kitagawa decomposition + Monte Carlo CIs\n")
cat("-----------------------------------\n")
source("scripts/02-analyze_data.R")
cat("\n")

cat("Step 3: Sensitivity analyses (7- & 16-cat GA, SUID split)\n")
cat("-----------------------------------\n")
source("scripts/03-sensitivity_analyses.R")
cat("\n")

cat("Step 4: Tests\n")
cat("-----------------------------------\n")
source("scripts/04-test_cause_categories.R")
cat("\n")

#### Summary ####
end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "secs")

cat("==============================================\n")
cat("Analysis pipeline complete\n")
cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total runtime:", round(runtime, 1), "seconds\n")
cat("==============================================\n\n")

cat("To generate the brief report, run:\n")
cat("  quarto render paper/paper_brief_report.qmd --to pdf\n")
cat("  (or paper/render_docx.sh for the Word/JAMA submission version)\n")
