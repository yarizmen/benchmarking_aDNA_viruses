# Check for correct (both at species level and at a higher taxa) and incorrect reads in old output.

configfile: 'config_check_correct_incorrect.yaml'

rule all:
    input:
        expand("Centrifuge_output/{sample}/{sample}_Correct_Incorrect.tsv", 
                sample=config['samples']),
        expand("Kraken2_output/{sample}/{sample}_Correct_Incorrect.tsv", 
                sample=config['samples']),
        expand("DIAMOND_output/{sample}/{sample}_Correct_Incorrect.tsv", 
                sample=config['samples']),
        expand("MetaPhlAn2_output/{sample}/{sample}_Correct_Incorrect.tsv", 
                sample=config['samples'])
        
rule check_ids:
    input:
        taxids = ancient("Higher_taxa_ids/{sample}_taxIDs.txt"),
        classifications = "{classifier}_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv"
    output:
        "{classifier}_output/{sample}/{sample}_Correct_Incorrect.tsv"
    params:
        runtime = "01:30:00",
        check_ids = config['path2check_ids']
    resources:
        memory = 2000
    log:
        "logs/check_correct/{classifier}/{sample}_check_ids.log"
    shell:
        '''
        awk -f {params.check_ids} {input.taxids} {input.classifications} > {output}
        '''

