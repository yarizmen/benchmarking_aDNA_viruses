#!/usr/bin/env Rscript

# Dealing with MetaPhlAn2 labels. From labels to taxIDs.

# March 2nd 2021
# Version 2: optparse included

#### Libraries ####
library(data.table)
library(taxize)
library(optparse)


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
                    function(x) get_uid(x, ask = F)[1], USE.NAMES = F, simplify = "vector")
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

                cat("Type the number of the row with the name you want to use ")
                Row_selection <- readLines(file("stdin"), n = 1)
                Row_selection <- as.numeric(Row_selection)

                MPA2_120v_names$V1[Check[counter]] <- To_resolve$gnr[Row_selection, 3][[1]]
                MPA2_taxID[Check[counter]] <- get_uid(MPA2_120v_names$V1[Check[counter]])[1]

                counter <- counter + 1

        }

}

write.table(MPA2_taxID, file = output_file, quote = F, col.names = F, row.names = F)
