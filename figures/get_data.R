folder="/Users/sneuensc/Sapfo/virome_yami/analyses_review"
setwd(folder)

library(data.table)
library(ggplot2)
library(taxonomizr)



## sum up (individual contigs)
classifiers <- c("Centrifuge", "Kraken2", "MALT", "DIAMOND", "MetaPhlAn2")
#classifiers <- c("Centrifuge", "Kraken2", "DIAMOND", "MetaPhlAn2")
#classifiers <- c("NCBI-acc2taxa", "MEGAN-acc2taxa","MEGAN-mapDB")

sets <- data.frame("set" = c(rep('deamination', 11), 
                             rep('lengths', 7),
                             rep('sequencing_errors', 6)),
                   "idx" = c(
                     c("deams_0", paste0("deams_point", c("05","1","15","2","25","3","35","4","45","5"))),
                     paste0("length_", c(30,40,50,60,90,120,150)),
                     c("seqerr_0", paste0("seqerr_minus", c(1,3,5,7,9)))
                     )
)

sets <- sets[sets$idx=="length_60",]

thres <- c(1,5,10,20,50)

sqlFile='accessionTaxa.sql'
desiredTaxa <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
rank.taxa <- c("Archaea","Bacteria","Eukaryota","Viruses","")

################################################################################
## remove double taxIDs
taxID.2049444 <- c("NC_040550.1","NC_040619.1","NC_040620.1","NC_040688.1","NC_040691.1","NC_040803.1","NC_040804.1","NC_040805.1","NC_040806.1")
taxID.944645 <- c("NC_043223.1","NC_043224.1","NC_043225.1","NC_043226.1","NC_043227.1","NC_043228.1") ## sum up
taxID.129875 <- c("AC_000005.1","NC_001460.1")
taxID.130310 <- c("AC_000006.1","NC_010956.1")
taxID.108098 <- c("NC_011202.1","NC_011203.1")
taxID.NC_022518 <- "NC_022518.1" ## Human endogenous retrovirus K113

## take one
for(s in classifiers){
  print(s)
  for(x in 1:nrow(sets)){
    set <- sets$set[x]
    idx <- sets$idx[x]
    file <- paste0(s,"/",set,"-",idx,"-count_table.tsv")
    d <- read.table(paste0(folder,"/",file), h=T)
    
    ## sum up the hits on contigs (and store sum in first and remove others (see below))
    d[,taxID.944645[1]] <- apply(d[, colnames(d) %in% taxID.944645,], 1, sum)
    
    ## take one
    d <- d[,! colnames(d) %in% c(taxID.2049444[-1],     ## take first
                                 taxID.944645[-1],      ## take first
                                 taxID.129875[1],       ## take second
                                 taxID.130310[1],       ## take second
                                 taxID.108098[2],       ## take first
                                 taxID.NC_022518)]
    
    write.table(d, paste0(folder,"/",s,"/",set,"-",idx,"-count_table2.tsv"), row.names = F, sep ="\t")
  }
}


## on the read count table
counts <- data.table(read.table("counts_orig.txt"))

## sum up and store it in first
counts[counts$V4 %in% taxID.944645[1],'V1'] <- counts[counts$V4 %in% taxID.944645,sum(V1), by=c('V2','V3')]$V1

## remove all but one
counts <- counts[! counts$V4 %in% c(taxID.2049444[-1],     ## take first
                             taxID.944645[-1],      ## take first
                             taxID.129875[1],       ## take second
                             taxID.130310[1],       ## take second
                             taxID.108098[2],       ## take first
                             taxID.NC_022518), ]    ## remove retrovirus

write.table(counts, paste0(folder, "/counts_orig2.txt"), row.names = F, col.names = F, sep ="\t")

################################################################################
## get classification categories
data <- NULL
for(s in classifiers){
  print(s)
  for(x in 1:nrow(sets)){
    set <- sets$set[x]
    idx <- sets$idx[x]
    file <- paste0(folder,"/",s,"/",set,"-",idx,"-count_table2.tsv")
    #print(file)
    d <- read.table(file, h=T)
    nb.taxa <- ncol(d)
    
    ## get higher ranks for all classified taxIDs
    tmp <- data.table(getTaxonomy(d$taxID, paste0(folder, "/", sqlFile), desiredTaxa=desiredTaxa))
    tmp[is.na(tmp)] = ""
    d <- cbind(d, tmp)

    sum.correct <- NULL
    sum.higher <- NULL
    sum.incorrect <- NULL
    sum.unclassified <- NULL
    extra.taxa <-NULL
    extra.rank<-NULL
    extra.taxa.reads <- NULL
    extra.taxa.taxa <- NULL
    higher.rank<-NULL
    for(i in 2:nb.taxa){
      taxa = colnames(d)[i]
      #print(taxa)
      
      ## get higher taxa
      true.taxa <- read.table(paste0(folder,"/Higher_taxa_ids2/", taxa, "_taxIDs.txt"), sep="\t")$V1

      sum.higher[i-1] <- sum(d[d$taxID %in% true.taxa[-1], i])
      sum.correct[i-1] <- sum(d[d$taxID %in% true.taxa[1], i])
      sum.incorrect[i-1] <- sum(d[!d$taxID %in% c(true.taxa, -2, 0), i])
      sum.unclassified[i-1] <- sum(d[d$taxID %in% c(-2,0), i])
      
      ## get it per superkingdom
      tmp.reads <- NULL
      tmp.taxa <- NULL
      for(r in 1:length(rank.taxa)){
        tmp <- d[!d$taxID %in% c(true.taxa, -2,0) & d$superkingdom == rank.taxa[r], i]
        
        ## get number of reads on extra taxa
        tmp.reads[r] <- sum(tmp)
        
        ## get number of extra taxa
        for(e in 1:length(thres)){
          tmp.taxa[length(tmp.taxa) + 1] <- sum(tmp >= e)
        }
      }
      extra.taxa.reads <- rbind(extra.taxa.reads, tmp.reads)
      extra.taxa.taxa <- rbind(extra.taxa.taxa, tmp.taxa)
      
      ## get extra taxa
      tmp <- NULL
      for(e in 1:length(thres)){
        tmp[e] <- sum(d[!d$taxID %in% c(true.taxa, -2,0), i] >= e)
      }
      extra.taxa <- rbind(extra.taxa, tmp)
    }
    
    ## combine the data
    colnames(extra.taxa.reads)<-paste0(rank.taxa, "_nb.reads")
    rownames(extra.taxa.reads)<-NULL
    
    names <- t(outer(paste0(rank.taxa, "_nb.taxa_"), thres, FUN = "paste0"))
    dim(names) <- NULL
    colnames(extra.taxa.taxa)<-names
    rownames(extra.taxa.taxa)<-NULL
    
    colnames(extra.taxa)<-paste(thres, "reads")
    rownames(extra.taxa)<-NULL
    
    tmp <- cbind(data.frame("Classifier"=s,
                            "Set"=set,
                            "Idx"=idx,
                            "GenBankID"=colnames(d)[2:nb.taxa], 
                            "Correct_Sp"=sum.correct, 
                            "Correct_HT"=sum.higher, 
                            "incorrect"=sum.incorrect, 
                            "unclassified"=sum.unclassified),
                  extra.taxa)
    tmp <- cbind(tmp, extra.taxa.reads)
    tmp <- cbind(tmp, extra.taxa.taxa)
    
    data <- rbind(data, tmp)
  }
}

## add total number of reads
counts <- read.table(paste0(folder,"/counts_orig2.txt"))
colnames(counts) <- c("Total", "Set", "Idx", "GenBankID")


## split length and the other sets to ge the total number of reads
data <- data.table(rbind(merge(data[data$Set == "lengths",], counts, all.x=T),
              merge(data[data$Set != "lengths",], counts[counts$Idx=="seqerr_minus9",c('GenBankID','Total')], all.x=T)))

################################################################################
## add unclassified
data$Unclassified1 <- data$Total - (data$Correct_Sp + data$Correct_HT + data$incorrect)

## get ratios
data$"Correct species" <- 100 * data$Correct_Sp / data$Total
data$"Correct higher" <- 100 * data$Correct_HT / data$Total
data$"Unclassified" <- 100 * data$Unclassified1 / data$Total
data$"Incorrect" <- 100 * data$incorrect / data$Total

data$Classifier <- factor(data$Classifier, levels=classifiers)

dataAll<-data
saveRDS(dataAll, paste0(folder, "/dataAll.rds"))
## dataAll<-readRDS(paste0(folder, "/dataAll.rds"))

################################################################################
################################################################################


## colours
colorcillos <- c(rgb(0, 158, 115, maxColorValue = 255),
                 rgb(181, 230, 0, maxColorValue = 255),
                 rgb(243, 155, 0, maxColorValue = 255),
                 rgb(213, 213, 213, maxColorValue = 255))

color.classifers <- c( "#d95f02", "#7570b3", "#666666", "#e7298a", "#a6761d")
if(length(classifiers)==4) color.classifers <- color.classifers[-3]

colors_bw <- c("black","gray")


## change legend size
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "cm"))
}


################################################################################
## write the spreadsheet
# library("xlsx")
# 
# counts <- read.table(paste0(folder,"/counts_orig2.txt"))
# colnames(counts) <- c("Total", "Set", "Idx", "GenBankID")
# 
# sets1 <- sets[sets$idx=="length_60",]
# for(s in classifiers){
#   print(s)
#   for(x in 1:nrow(sets1)){
#     set <- sets$set[x]
#     idx <- sets$idx[x]
#     file <- paste0(folder,"/",s,"/",set,"-",idx,"-count_table2.tsv")
#     d <- read.table(file, h=T)
#     
#     ## remove strange rows
#     d <- d[d$taxID>0,]
#     
#     tmp <- apply(d[,-1], 2, sum)
#     tmp <- data.frame("GenBankID" = names(tmp), "Classified" = tmp)
#     rownames(tmp1) <- tmp1$GenBankID
#     
#     tmp1 <- merge(tmp, counts[counts$Idx == "length_60",c('GenBankID','Total')], all.x=T)
#     tmp1$unclassified <- tmp1$Total - tmp1$Classified
#     
#     print(range(tmp1$unclassified))
#     
#     rownames(tmp1) <- tmp1$GenBankID
#     
#     tmp2 <- cbind(data.frame('taxID'= c(0,0,0,0)), t(tmp1))
#     
#     tmp3 <- rbind(tmp2[4,], d)
#     
#     file <- paste0(folder,"/table_S02_viral_classification.xlsx")
#     write.xlsx(tmp3, file, sheetName = s, col.names = TRUE, row.names = F, append = T)
#   }
# }