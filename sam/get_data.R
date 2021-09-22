folder="/Users/sneuensc/Sapfo/virome_yami/analyses_review"
setwd(folder)

library(data.table)
library(ggplot2)

classifiers <- c("Kraken2", "Centrifuge", "DIAMOND", "Malt", "MetaPhlAn2")
sets <- c("lengths","deamination","sequencing_errors")

sets <- data.frame("set" = c(rep('deamination', 11), 
                             rep('lengths', 7),
                             rep('sequencing_errors', 6)),
                   "idx" = c(
                     c("deams_0", paste0("deams_point", c("05","1","15","2","25","3","35","4","45","5"))),
                     paste0("length_", c(30,40,50,60,90,120,150)),
                     c("seqerr_0", paste0("seqerr_minus", c(1,3,5,7,9)))
                     )
)
#sets<-sets[1:2,]

thres <- c(1,5,10,20,50)

## get classification categories
data <- NULL
for(s in classifiers){
  print(s)
  for(x in 1:nrow(sets)){
    set <- sets$set[x]
    idx <- sets$idx[x]
    file <- paste0(s,"/",set,"-",idx,"-count_table.tsv")
    #print(file)
    d <- read.table(file, h=T)
    
    sum.correct <- NULL
    sum.higher <- NULL
    sum.incorrect <- NULL
    sum.unclassified <- NULL
    extra.taxa <-NULL
    for(i in 2:ncol(d)){
      taxa = colnames(d)[i]
      #print(taxa)
      
      ## get higher taxa
      true.taxa <- read.table(paste0("Higher_taxa_ids/", taxa, "_taxIDs.txt"), sep="\t")$V1
      
      sum.higher[i-1] <- sum(d[d$taxID %in% true.taxa[-1], i])
      sum.correct[i-1] <- sum(d[d$taxID %in% true.taxa[1], i])
      sum.incorrect[i-1] <- sum(d[!d$taxID %in% c(true.taxa, -2,0), i])
      sum.unclassified[i-1] <- sum(d[d$taxID %in% c(-2,0), i])
      
      ## get extra taxa
      tmp <- NULL
      for(e in 1:length(thres)){
        tmp[e] <- sum(d[!d$taxID %in% c(true.taxa, -2,0), i] >= e)
      }
      extra.taxa <- rbind(extra.taxa, tmp)
    }
    
    ## combine the data
    colnames(extra.taxa)<-paste(thres, "reads")
    rownames(extra.taxa)<-NULL
    tmp <- cbind(data.frame("Classifier"=s,
                            "Set"=set,
                            "Idx"=idx,
                            "GenBankID"=colnames(d)[-1], 
                            "Correct_Sp"=sum.correct, 
                            "Correct_HT"=sum.higher, 
                            "Incorrect"=sum.incorrect, 
                            "unclassified"=sum.unclassified), extra.taxa)
    
    data <- rbind(data, tmp)
  }
}

## add total number of reads
counts <- read.table("counts_orig.txt")
colnames(counts) <- c("Total", "Set", "Idx", "GenBankID")

## split length and the other sets to ge the total number of reads
data <- rbind(merge(data[data$Set == "lengths",], counts, all.x=T),
              merge(data[data$Set != "lengths",], counts[counts$Idx=="seqerr_minus9",c('GenBankID','Total')], all.x=T))


## add unclassified
data$Unclassified1 <- data$Total - (data$Correct_Sp + data$Correct_HT + data$Incorrect)

## get ratios
data$"Correct species" <- 100 * data$Correct_Sp / data$Total
data$"Correct higher" <- 100 * data$Correct_HT / data$Total
data$"Unclassified" <- 100 * data$Unclassified1 / data$Total
data$"Incorrect" <- 100 * data$Incorrect / data$Total

data$Classifier <- factor(data$Classifier, levels=classifiers)

dataAll<-data.table(data)



## number of unclassified per classifier
data <- data.table(data)
data[,list("unclassified"=mean(unclassified)), by=Classifier]

## colours
colorcillos <- c(rgb(0, 158, 115, maxColorValue = 255),
                 rgb(181, 230, 0, maxColorValue = 255),
                 rgb(243, 155, 0, maxColorValue = 255),
                 rgb(213, 213, 213, maxColorValue = 255))

color.classifers <- c( "#56B4E9","#0072B2", "#D55E00", "#CC79A7", "orange ")

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