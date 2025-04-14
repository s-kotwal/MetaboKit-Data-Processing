# *Please cite [paper] and github repository if being used.*

# Presented is a code used to automatically filter .tsv files "Report_RTseparated" or "Report_RTseparated_fill" that are auto-generated from processing Mass Spec .mzML files using MetaboKit. The final output is a file that removes duplicate annotations, filter compounds based blanks, abundance cutoffs, etc. 
# Must have blank and experimental samples in MetaboKit .tsv as 'pos-sample-it1' to indicate mode, sample name, and number of iterative injection.
# This code takes into account rows where if a blank has an annotation, it is included only if one of the experimental samples is 3x greater than the average of non-zero blanks.
# Need to define mode as "positive" or "negative", and .tsv file exclusively.

# Install necessary packages if not already installed
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# Load necessary libraries
library(dplyr)
library(openxlsx)

# Define mode
mode <- "positive"  # Change to "negative" if needed

# Prefix definitions based on mode
mode_prefix <- ifelse(mode == "positive", "pos", "neg")
blank_sample_prefix <- paste0(mode_prefix, ".blank") # Must have blank samples named 'pos-blank-it1' in .tsv file.
experimental_sample_prefixes <- c(paste0(mode_prefix, "sample.name")) # List experimental sample names. '.' are read as hyphens in .tsv files
num_injections <- 5  # Number of iterative injections

# Generate column names for blank and experimental samples
blank_cols <- paste0(blank_sample_prefix, ".it", 1:num_injections)
experimental_cols <- unlist(lapply(experimental_sample_prefixes, function(prefix) {
  paste0(prefix, ".it", 1:num_injections)
}))

# Read the TSV file
file_path <- "C:/your_file_path/Report_RTseparated_fill.tsv"
data <- read.delim(file_path, sep = "\t")

# Replace NA or empty values in the specified columns with 0
data <- data %>%
  mutate(across(all_of(c(blank_cols, experimental_cols)), ~ ifelse(is.na(.) | . == "", 0, .)))

# Filtering rows based on the criteria
blank_alignment <- data %>%
  rowwise() %>% # Apply function to each row
  filter({
    # Calculate the average of non-zero blank values
    blank_values <- c_across(all_of(blank_cols))
    non_zero_blanks <- blank_values[blank_values > 0]
    avg_blank <- if (length(non_zero_blanks) > 0) mean(non_zero_blanks) else 0
    
    # Check if at least one experimental value is at least 3 times greater than the average blank
    experimental_values <- c_across(all_of(experimental_cols))
    keep_row <- any(experimental_values >= 3 * avg_blank) || avg_blank == 0
    
    keep_row
  }) %>%
  ungroup()

# Create a function to get the maximum value and its index among the iterative injections
get_max_value_and_index <- function(data, cols) {
  max_values <- suppressWarnings(apply(data[, cols], 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))))
  max_indices <- suppressWarnings(apply(data[, cols], 1, function(x) ifelse(all(is.na(x)), NA, which.max(x))))
  list(max_values = max_values, max_indices = max_indices)
}

# Add columns for the aggregated values
aggregated_columns <- list()

# Aggregate blank and experimental columns
aggregated_columns[[blank_sample_prefix]] <- get_max_value_and_index(blank_alignment, blank_cols)$max_values
for (prefix in experimental_sample_prefixes) {
  aggregated_columns[[prefix]] <- get_max_value_and_index(blank_alignment, paste0(prefix, ".it", 1:num_injections))$max_values
}

# Add aggregated values to the dataframe
for (col_name in names(aggregated_columns)) {
  blank_alignment[[col_name]] <- aggregated_columns[[col_name]]
}

# Add columns for RT, score, mp, and S/N for each sample prefix using the max index to extract corresponding values
attributes <- c("RT_", "score_", "mp_", "S.N_")
for (prefix in experimental_sample_prefixes) {
  # Get the index of the maximum abundance value for the current experimental prefix
  max_indices <- get_max_value_and_index(blank_alignment, paste0(prefix, ".it", 1:num_injections))$max_indices
  
  # Extract values for each attribute based on the max_indices
  for (attr in attributes) {
    group_cols <- paste0(attr, prefix, ".it", 1:num_injections)
    if (all(group_cols %in% colnames(blank_alignment))) {
      # Add a new column with the selected value for the corresponding attribute
      blank_alignment[[paste0(attr, prefix)]] <- mapply(function(idx, row) {
        if (!is.na(idx)) {
          row[[group_cols[idx]]]
        } else {
          NA
        }
      }, max_indices, split(blank_alignment, seq(nrow(blank_alignment))))
    }
  }
}

# Explicitly remove the original iterative injection columns
all_iterative_columns <- grep("\\.it\\d+$", colnames(blank_alignment), value = TRUE)

# Remove the iterative columns
blank_alignment <- blank_alignment %>%
  select(-all_of(all_iterative_columns))

# Replace #NAME? and -Inf with NA
blank_alignment[blank_alignment == "#NAME?" | blank_alignment == "-Inf"] <- NA

# Remove duplicates in the "InChIKey" column, keeping the row with the highest value in the specified columns
blank_alignment <- blank_alignment %>%
  group_by(InChIKey) %>%
  filter(row_number(desc(pmax(!!!syms(c(blank_sample_prefix, experimental_sample_prefixes))))) == 1) %>%
  ungroup()

# Additional filtering based on "accession..HMDB." and "adduct"
filtered_data <- blank_alignment %>%
  mutate(accession..HMDB. = as.character(accession..HMDB.)) %>%
  filter(is.na(accession..HMDB.) | accession..HMDB. == "") %>%
  bind_rows(
    blank_alignment %>%
      filter(!is.na(accession..HMDB.) & accession..HMDB. != "") %>%
      dplyr::group_by_at(vars(accession..HMDB., adduct)) %>%
      filter({
        if (n() > 1) {
          if (max(peak_shape..median., na.rm = TRUE) != min(peak_shape..median., na.rm = TRUE)) {
            peak_shape..median. == max(peak_shape..median., na.rm = TRUE)
          } else {
            score_average <- rowMeans(pick(starts_with("score_")), na.rm = TRUE)
            if (all(is.na(score_average))) {
              TRUE
            } else if (max(score_average, na.rm = TRUE) != min(score_average, na.rm = TRUE)) {
              score_average == max(score_average, na.rm = TRUE)
            } else {
              sn_average <- rowMeans(pick(starts_with("S.N_")), na.rm = TRUE)
              if (all(is.na(sn_average))) {
                TRUE
              } else {
                sn_average == max(sn_average, na.rm = TRUE)
              }
            }
          }
        } else {
          TRUE
        }
      }) %>%
      ungroup()
  )




# Rename 'group' to 'ID' and create a new column with IDs starting from 1
filtered_data <- filtered_data %>%
  mutate(ID = row_number())

# Reorder columns and filter them
filtered_dataset <- filtered_data %>%
  select(ID, name, all_of(blank_sample_prefix), !!!syms(experimental_sample_prefixes),
         Median.RT, name..HMDB., formula, InChIKey, feature_m.z, adduct,
         super.class..HMDB., class..HMDB., sub.class..HMDB., accession..HMDB.,
         confidence.level, Min..RT, Median.RT, Max..RT,
         matching_peaks..median., X.lib_peaks..median., peak_shape..median., X.detected,
         starts_with("RT_"), starts_with("score_"), starts_with("mp_"), starts_with("S.N_"))

# Remove empty columns
filtered_dataset <- filtered_dataset %>%
  select_if(function(col) !all(is.na(col)))

# Create the "Notes" data
notes_header <- c(paste("All samples analyzed by", num_injections, "iterative injections"), "")
notes_data <- data.frame(
  Field = c("ID", "name", blank_sample_prefix, experimental_sample_prefixes, " ", "Median.RT",
            "name..HMDB.", "formula", "InChIKey", "feature.m.z", "adduct", "super.class..HMDB.", "class..HMDB.",
            "sub.class..HMDB.", "accession..HMDB.", "confidence.level", "Min..RT", "Max..RT", "matching_peaks..median.",
            "X.lib_peaks..median.", "peak_shape..median.", "X.detected", "RT_", "score_",
            "mp_", "S.N_"),
  Description = c("ID from filtered dataset", "Compound name", 
                  "Sample containing only extraction solvent", rep("Experimental samples in positive mode", length(experimental_sample_prefixes)), "", "RT = Retention time in minutes",
                  "Human Metabolome Database", "", "", "Observed m/z of compound", "Observed ion of compound", "", "",
                  "", "", "ID by MS1 or MS2 (MS/MS)", "", "", "", "", "", "",
                  "RT = Retention time in minutes", "Score of library and experimental spectra", "", "Signal to noise")
)

# Add header and notes data to the workbook
wb <- createWorkbook()
addWorksheet(wb, "MetaboKit_filtered_dataset")
writeData(wb, "MetaboKit_filtered_dataset", filtered_dataset)
addWorksheet(wb, "Notes")
writeData(wb, "Notes", notes_header, startCol = 1, startRow = 1)
writeData(wb, "Notes", notes_data, startCol = 1, startRow = 2)

# Create the "Abundance annotated" sheet
abundance_annotated_header <- paste("Highest abundance of compound in", num_injections, "iterative injections")
abundance_annotated_data <- filtered_dataset %>%
  select(ID, name, all_of(blank_sample_prefix), !!!syms(experimental_sample_prefixes))

addWorksheet(wb, "Abundance annotated")
writeData(wb, "Abundance annotated", abundance_annotated_header, startCol = 1, startRow = 1)
writeData(wb, "Abundance annotated", abundance_annotated_data, startCol = 1, startRow = 2)

# Create the "Unique" sheet
unique_header <- "Unique compound identified in sample"
unique_data <- abundance_annotated_data

# Create a copy to work on unique identification
unique_copy <- unique_data

# Initialize all entries as empty strings
unique_copy[, c(experimental_sample_prefixes)] <- ""

# Check for uniqueness
for (col in c(blank_sample_prefix, experimental_sample_prefixes)) {
  other_cols <- setdiff(c(experimental_sample_prefixes), col)
  
  unique_copy[[col]] <- ifelse(
    unique_data[[col]] > 0 & rowSums(unique_data[other_cols] > 0) == 0,
    "x",
    ""
  )
}

# Write the unique data to the workbook
addWorksheet(wb, "Unique")
writeData(wb, "Unique", unique_header, startCol = 1, startRow = 1)
writeData(wb, "Unique", unique_copy, startCol = 1, startRow = 2)

# Add a new sheet called "Raw Data" with the raw dataset from the original file
addWorksheet(wb, "Raw Data")
writeData(wb, "Raw Data", data)

# Save the workbook with the added "Raw Data" sheet as the last sheet
saveWorkbook(wb, paste0(mode_prefix, "_filename_.xlsx"), overwrite = TRUE)
