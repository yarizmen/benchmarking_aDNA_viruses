#!/usr/bin/env Rscript

# Classifier comparison.
# Difference between Classifier_comparison.R & this script is that the statistics are computed 
# different. In the present script the statistic is computed per each virus, and then the mean of 
# all statistics is used.  
# I'll take advantage of the "<classifier>_RawNumbers.txt" tables.

# Actions:
# Requires "<classifier>_RawNumbers.txt" & "<classifier>_output/Reported_virus_<classifier>.txt"
# files as input.
# Computes stats (proportions [classified, correct species, correct higher taxa, incorrect],
# sensitivity, & precision) per virus, then computes mean.
# Computes number of identified virus & extra taxa reported
# Stores the stats into a table.
# All above per classifier.

# Running example:
# Rscript --vanilla Classifier_comparison_mean_of_means.R 
# -p "walle_virome/simulations_Feb_2020b/Lengths/length_30/"
# -c "c(\"Centrifuge\", \"Kraken2\", \"DIAMOND\", \"MetaPhlAn2\")"

#### Libraries #####
library(tidyr)
library(data.table)
library(optparse)


#### Arguments ####
option_list = list(
  make_option(c("-p", "--path"), type = "character", default = NULL, 
              help = "path to data", metavar = "character"),
  make_option(c("-c", "--classifiers"), type = "character", 
              default = "c(\"Centrifuge\", \"Kraken2\", \"DIAMOND\", \"MetaPhlAn2\")",
              help = "vector with classifier names")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

path2data <- opt$path
classifiers <- eval(parse(text = (opt$classifiers)))

#### Function to compute desired statistics ####
compute_statistics <- function(classifier, path2data) {
    raw_numbers <- fread(file = paste0(path2data, classifier, "_120v_RawNumbers.txt"))

    proportion_classified <- raw_numbers[, mean((100 * class)/TotalSimReads)]
    correct_species <- raw_numbers[, mean((100 * Correct_Sp)/TotalSimReads)]
    correct_higher <- raw_numbers[, mean((100 * Correct_HT)/TotalSimReads)]
    incorrect <- raw_numbers[, mean((100 * Incorrect)/TotalSimReads)]
    precision <- raw_numbers[, ((Correct_Sp/class) * 100) %>% replace_na(., 0) %>% mean]
    sensitivity <- raw_numbers[, mean((Correct_Sp/TotalSimReads) * 100)]
    
    reported_virus <- fread(file = paste0(path2data, classifier, "_output/Reported_virus_", 
                                          classifier, ".txt"), 
                            header = F, select = c(1))
    
    identified <- raw_numbers[, sum(Correct_Sp != 0)]
    extra_taxa <- reported_virus[1, V1] - identified
    
    stats <- c(proportion_classified, correct_species, correct_higher, incorrect, precision,
               sensitivity, identified, extra_taxa)
    names(stats) <- c("prop_classified", "correct_sp", "correct_ht", "incorrect", "precision",
                      "sensitivity", "identified", "extra_taxa")
    
    return(stats)
}


#### Main ####

stats_table <- sapply(X = classifiers, 
                      FUN = function(x) {compute_statistics(classifier = x, 
                                                            path2data = path2data)})

write.table(x = stats_table, file = paste0(path2data, "stats_means_classifiers.csv"), sep = ",")

