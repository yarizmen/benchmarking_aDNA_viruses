#!/usr/bin/env Rscript

# Dealing with MetaPhlAn2 labels. From labels to taxIDs.

# March 2nd 2021
# Version 2: optparse included

#### Libraries ####
library(data.table)
.libPaths("/home/yarizmen/R/x86_64-koji-linux-gnu-library/3.5/")
library(taxize)
library(optparse)
sessionInfo()

#### Arguments ####
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "File with MetaPhlAn2 names", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file with the taxIDs")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

input_MPA2_names <- opt$input
output_file <- opt$output

#### Main ####

MPA2_120v_names <- fread(input = input_MPA2_names, header = F, sep = ",")
message("There are ", length(MPA2_120v_names$V1), " taxa in LastTaxRank file")

#Look for taxIDs
MPA2_taxID <- c(1:length(MPA2_120v_names$V1))

MPA2_taxID <- sapply(MPA2_120v_names$V1, 
                    function(x) get_uid(x, key = "8f5a7ee70253e037c5641bab48b0d6d74908", 
                                        ask = T)[1], USE.NAMES = F, simplify = "vector")
Check <- which(is.na(MPA2_taxID))


if(!(length(Check) < 1)){
        message(length(Check), " taxa names must be resolved.")

        counter <- 1

        for (i in Check){
                message(MPA2_120v_names$V1[Check[counter]], " must be resolved.")
                message("Looking for possible names: ")

                To_resolve <- resolve(query = MPA2_120v_names$V1[Check[counter]])
                print(To_resolve)

                message("Possible names: ")

                print(To_resolve$gnr[, 3])

                Row_selection <- readline(prompt = "Type the number of the row with the name you want to use ")
                MPA2_120v_names$V1[Check[counter]] <- To_resolve$gnr[Row_selection, 3]
                MPA2_taxID[Check[counter]] <- get_uid(MPA2_120v_names$V1[Check[counter]], 
                                                    key = "8f5a7ee70253e037c5641bab48b0d6d74908")[1]

                counter <- counter + 1

        }

}

write.table(MPA2_taxID, file = output_file, quote = F, col.names = F, row.names = F)