#!/usr/bin/env Rscript

# Comparing outputs by molecule type.
# This script groups the counts from "<viral_sequence>_Correct_Incorrect.tsv" 
# (given a list of viral sequences), then computes the proportions of correctly 
# and incorrectly classified reads using "total_simulated_reads_accession.tsv", 
# and finally gets the means of those proportions.

# Requirements:
# A list of the accsession numbers of the viral sequences of interest, 
# "<viral_sequence>_Correct_Incorrect.tsv" 
# & "total_simulated_reads_accession.tsv" 
# files as input.

# Actions:
# Groups counts and stores them into data.table.
# Computes proportions.
# Computes means.

# Running example:
# Rscript --vanilla compare_output.R -v "c(\"NC_001338.1\", \"NC_008517.1\", 
# \"NC_008976.1\", \"NC_026238.1\")" -s "total_simulated_reads_accession.tsv"


#### Libraries #####
library(data.table)
library(tidyverse)
library(optparse)

#### Arguments ####
option_list = list(
  make_option(c("-v", "--viral_sequences"), type = "character", default = NULL, 
              help = "viral sequences' accession numbers", 
              metavar = "character"),
  make_option(c("-s", "--simulated_reads"), type = "character", 
              default = NULL,
              help = "table with accession number and number of simulated reads")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

viral_sequences <- eval(parse(text = (opt$viral_sequences)))
simulated_reads <- opt$simulated_reads
simulated_reads <- fread(input = simulated_reads)

#### Function to group counts ####
group_compute <- function(viral_seq, dt_core){
  counts <- fread(input = paste0(viral_seq, "/", viral_seq, 
                                 "_Correct_Incorrect.tsv"))
  grouped <- counts[, .(total = sum(V1)), by = V6]
  grouped <- dcast(melt(grouped, id.vars = "V6"), variable ~ V6)
  grouped[, variable := NULL]
  grouped[, accession := viral_seq]
  
  if (!("Unclassified" %in% colnames(grouped))) {
    merged <- merge(dt_core, grouped, all = T, 
                    by = colnames(grouped))
  } else {
    merged <- rbind(dt_core, list(accession = viral_seq, 
                                  Correct_species = 0, 
                                  Correct_higher = 0, 
                                  Incorrect = 0))
  }
  
  return(merged)
}

#### Main ####

# data.table skeleton
dt_core <- data.table(c(NA), c(NA), c(NA), c(NA))
colnames(dt_core) <- c("accession", "Correct_species", "Correct_higher", 
                       "Incorrect")
# Group counts
for (i in viral_sequences) {
  dt_core <- group_compute(viral_seq = i, dt_core = dt_core)
}

# Remove first row and change NAs to 0
dt_core <- dt_core[-1]
dt_core[is.na(dt_core)] <- 0

# Add total number of reads
dt_simreads <- simulated_reads[V1 %in% viral_sequences, ]
colnames(dt_simreads) <- c("accession", "total_reads")

dt_core <- merge(dt_core, dt_simreads, by = "accession")

# Proportions
dt_core[, prop_correct_sp := (Correct_species * 100)/total_reads]
dt_core[, prop_correct_ht := (Correct_higher * 100)/total_reads]
dt_core[, prop_incorrect := (Incorrect * 100)/total_reads]
dt_core[, unclassified := total_reads - (Correct_species + Correct_higher + Incorrect)]
dt_core[, prop_unclassified := (unclassified * 100/total_reads)]

print(dt_core[, .(mean(prop_correct_sp), mean(prop_correct_ht), 
                  mean(prop_incorrect), mean(prop_unclassified))])