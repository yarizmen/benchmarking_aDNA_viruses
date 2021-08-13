#!/usr/bin/env Rscript

# Merge the outputs per virus into a single table.
# This script groups the counts from "<classsifer>_output/<sequence>/<sequence>_ReadspTaxon.tsv" 
# into a single table.

# Requirements:
# A list of the accsession numbers of the viral sequences of interest, 
# & "<classsifer>_output/<sequence>/<sequence>_ReadspTaxon.tsv" files as input.

# Actions:
# Merge the different output into a single table.
# If output is empty, creates empty table (first file) or empty rows (subsequent files).

# Running example:
# Rscript --vanilla merge_into_counts_table.R -o "Centrifuge_output/counts_table.tsv"


#### Libraries #####
library(data.table)
library(purrr)
library(stringr)
library(tidyr)
library(optparse)

#### Arguments ####
option_list = list(
  make_option(c("-o", "--output_file"), type = "character", default = NULL, 
              help = "output file name", 
              metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

output_file <- opt$output_file

#### Main ####

### First sample
# Get classifier name
classifier <- str_split(output_file, pattern = "/")[[1]][1] %>% str_split(., pattern = "_")
classifier <- classifier[[1]][1]

# Get sequences name
seqs_paths <- list.files(path = paste0(classifier, "_output/"), pattern = "*_ReadspTaxon\\.txt", 
                        recursive = T, full.names = F)

viral_sequences <- str_split(seqs_paths, pattern = "/") %>% 
                  map_depth(., 1, 2) %>% 
                  unlist %>% 
                  str_replace(., "_ReadspTaxon\\.txt", "")

sample <- paste0(classifier, "_output/", viral_sequences[1], "/", viral_sequences[1], 
                "_ReadspTaxon.txt")

# Check if file is not empty
if (file.info(sample)[1, "size"] != 0) {
  counts <- fread(input = sample)
  colnames(counts) <- c(viral_sequences[1], "taxID")
  
} else { # if empty, create an empty table
  counts <- data.table(0, 10239)
  colnames(counts) <- c(viral_sequences[1], "taxID")
}

### Rest of samples

for (i in c(2:length(viral_sequences))) {
  # Reads per taxon file
  sample <- paste0(classifier, "_output/", viral_sequences[i], "/", viral_sequences[i], 
                "_ReadspTaxon.txt")
  
  # Check if it is not empty
  if (file.info(sample)[1, "size"] != 0) {
    # Read table
    reads_per_taxon <- fread(input = sample)
    colnames(reads_per_taxon) <- c(viral_sequences[i], "taxID")
    # Merge with existing table
    counts <- merge(counts, reads_per_taxon, all = T, by = "taxID")
  } else {
    # Empty file: add row with zeros
    counts[, viral_sequences[i] := rep(0, nrow(counts))]
  }
}

# NAs to 0
counts[is.na(counts)] <- 0

write.table(counts, file = output_file, quote = F, col.names = T, row.names = F, sep = "\t")