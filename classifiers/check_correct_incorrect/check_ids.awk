# Check if the classifications made were correct (at species or higher taxa) or incorrect. If 
# correct, indicate how far in the taxonomy they are respect to the sequence taxID.

# Required:
# File with taxIDs in a single column.
# {sample}_ReadsTaxon_NamesRanks.tsv classifications table from "run_classifiers.smk"

# Use:
# awk -f check_ids.awk example_taxIDs.txt example_ReadsTaxon_NamesRanks.tsv

BEGIN{FS="\t"; level=0}

{

if(NR==FNR){
    taxIDs[$1]=level; 
    level+=1;
    next
}

if($2 in taxIDs){
    if(taxIDs[$2]==0){
        print $0"\tCorrect_species\t"taxIDs[$2]
    }else{
        print $0"\tCorrect_higher\t"taxIDs[$2]
    }     
}else{
        print $0"\tIncorrect\t"0
}

}

