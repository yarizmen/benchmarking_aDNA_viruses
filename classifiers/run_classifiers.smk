# Run Centrifuge, Kraken2, DIAMOND, and MetaPhlAn2 for a given sample

configfile: 'config_run_classifiers.yaml'

rule all:
    input:
        expand("Centrifuge_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv", 
                sample=config['samples']),
        expand("Kraken2_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv", 
                sample=config['samples']),
        expand("DIAMOND_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv", 
                sample=config['samples']),
        expand("MetaPhlAn2_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv", 
                sample=config['samples'])

rule get_fasta:
    input:
        "samples/fastq_files/{sample}.fq.gz"
    output:
        "samples/fasta_files/{sample}.fna"
    resources:
        memory = 5000
    params:
        runtime = "04:00:00"
    log:
        "logs/{sample}/{sample}_get_fasta.log"
    shell:
        '''
        module add UHTS/Analysis/seqtk/1.2;
        seqtk seq -a {input} > {output}
        '''

rule run_Centrifuge:
    input:
        "samples/fasta_files/{sample}.fna"
    output:
        report_file = "Centrifuge_output/{sample}/{sample}_report_file.tsv",
        classifications = "Centrifuge_output/{sample}/{sample}_classifications.txt"
    params:
        runtime = "10:00:00",
        Centrifuge = config['path2Centrifuge'],
        database = config['path2Centrifuge_db'],
        assign = 1
    resources:
        memory = 30000
    threads: 5
    benchmark:
        "benchmarks/Centrifuge/{sample}/{sample}.benchmark.txt"
    log:
        "logs/Centrifuge/{sample}/{sample}_run_Centrifuge.log"
    shell:
        '''
        {params.Centrifuge} -x {params.database} -f -U {input} -p {threads} -k {params.assign} \
        --report-file {output.report_file} -S {output.classifications}
        '''

rule reads_per_taxon_Centrifuge:
    input:
        "Centrifuge_output/{sample}/{sample}_classifications.txt"
    output:
        "Centrifuge_output/{sample}/{sample}_ReadspTaxon.txt"
    params:
        runtime = "02:00:00",
        empty_ReadsTaxon_NamesRanks = "Centrifuge_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv",
        empty_Correct_Incorrect = "Centrifuge_output/{sample}/{sample}_Correct_Incorrect.tsv"
    resources:
        memory = 3000
    log:
        "logs/Centrifuge/{sample}/{sample}_reads_per_taxon_Centrifuge.log"
    shell:
        '''
        cut -f3 {input} | grep -v 'taxID' | sort | uniq -c | awk '{{print $1"\t"$2}}' > {output}

        if ! [ -s {output} ]; then 
                awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0}}' > {params.empty_ReadsTaxon_NamesRanks};
                awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0"\\tUnclassified\\t"0}}' > \
                {params.empty_Correct_Incorrect}
        fi
        '''

rule run_Kraken2:
    input:
        "samples/fasta_files/{sample}.fna"
    output:
        K2out = "Kraken2_output/{sample}/{sample}_K2out.txt",
        report_file = "Kraken2_output/{sample}/{sample}_K2report.txt"
    params:
        runtime = "10:00:00",
        Kraken2 = config['path2Kraken2'],
        database = config['path2Kraken2_db']
    resources:
        memory = 48000
    threads: 5
    benchmark:
        "benchmarks/Kraken2/{sample}/{sample}.benchmark.txt"
    log:
        "logs/Kraken2/{sample}/{sample}_run_Centrifuge.log"
    shell:
        '''
        {params.Kraken2} --db {params.database} {input} --threads {threads} \
        --report {output.report_file} --use-names > {output.K2out}
        '''

rule reads_per_taxon_Kraken2:
    input:
        "Kraken2_output/{sample}/{sample}_K2out.txt"
    output:
        "Kraken2_output/{sample}/{sample}_ReadspTaxon.txt"
    params:
        runtime = "02:00:00",
        empty_ReadsTaxon_NamesRanks = "Kraken2_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv",
        empty_Correct_Incorrect = "Kraken2_output/{sample}/{sample}_Correct_Incorrect.tsv"
    resources:
        memory = 3000
    log:
        "logs/Kraken2/{sample}/{sample}_reads_per_taxon_Kraken2.log"
    shell:
        '''
        set +e

        grep '^C' {input} | cut -f3 | grep -oP '\(taxid \d+\)' | egrep -o '[0-9]+' | sort | uniq -c | \
        awk '{{print $1"\t"$2}}' > {output}

        exitcode=${{PIPESTATUS[0]}}
        if [[ exitcode -eq 1 ]]; then
            awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0}}' > {params.empty_ReadsTaxon_NamesRanks};
            awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0"\\tUnclassified\\t"0}}' > \
            {params.empty_Correct_Incorrect}
            exit 0
        fi
        '''

rule run_DIAMOND:
    input:
        "samples/fasta_files/{sample}.fna"
    output:
        "DIAMOND_output/{sample}/{sample}_tax.txt"
    params:
        runtime = "10:00:00",
        DIAMOND = config['path2DIAMOND'],
        mode = 'blastx',
        database = config['path2DIAMOND_db']
    resources:
        memory = 24000
        #constraint = 'AVX512' 
    threads: 5
    benchmark:
        "benchmarks/DIAMOND/{sample}/{sample}.benchmark.txt"
    log:
        "logs/DIAMOND/{sample}/{sample}_run_DIAMOND.log"
    shell:
        '''
        {params.DIAMOND} {params.mode} -p {threads} -d {params.database} -q {input} -o {output} \
        -f 102
        '''

rule reads_per_taxon_DIAMOND:
    input:
        "DIAMOND_output/{sample}/{sample}_tax.txt"
    output:
        "DIAMOND_output/{sample}/{sample}_ReadspTaxon.txt"
    params:
        runtime = "02:00:00",
        empty_ReadsTaxon_NamesRanks = "DIAMOND_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv",
        empty_Correct_Incorrect = "DIAMOND_output/{sample}/{sample}_Correct_Incorrect.tsv"
    resources:
        memory = 3000
    log:
        "logs/DIAMOND/{sample}/{sample}_reads_per_taxon_DIAMOND.log"
    shell:
        '''
        set +e
        
        cut -f2 {input} | grep -wv '0' | sort | uniq -c | awk '{{print $1"\t"$2}}' > {output}
        
        exitcode=${{PIPESTATUS[1]}}
        if [[ exitcode -eq 1 ]]; then
            awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0}}' > {params.empty_ReadsTaxon_NamesRanks};
            awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0"\\tUnclassified\\t"0}}' > \
            {params.empty_Correct_Incorrect}
            exit 0
        fi
        '''

rule run_MetaPhlAn2:
    input:
        "samples/fasta_files/{sample}.fna"
    output:
        bt2o = "MetaPhlAn2_output/{sample}/{sample}.bowtie2.bz2",
        MPA2_out = "MetaPhlAn2_output/{sample}/{sample}_out.txt"
    params:
        runtime = "10:00:00",
        MetaPhlAn2 = config['path2MetaPhlAn2'],
        i_type = 'fasta',
        o_type = 'reads_map'
    resources:
        memory = 12000
    threads: 5
    benchmark:
        "benchmarks/MetaPhlAn2/{sample}/{sample}.benchmark.txt"
    log:
        "logs/MetaPhlAn2/{sample}/{sample}_run_MetaPhlAn2.log"
    shell:
        '''
        module add UHTS/Aligner/bowtie2/2.3.4.1;

        {params.MetaPhlAn2} {input} --bowtie2out {output.bt2o} --nproc {threads} \
        --input_type {params.i_type} -t {params.o_type} -o {output.MPA2_out}
        '''

rule reads_per_taxon_MetaPhlAn2:
    input:
        "MetaPhlAn2_output/{sample}/{sample}_out.txt"
    output:
        counts = temp("MetaPhlAn2_output/{sample}/{sample}_counts.txt"),
        MPA2_names = temp("MetaPhlAn2_output/{sample}/LastTaxRank")
    params:
        runtime = "01:00:00",
        empty_ReadsTaxon_NamesRanks = "MetaPhlAn2_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv",
        empty_Correct_Incorrect = "MetaPhlAn2_output/{sample}/{sample}_Correct_Incorrect.tsv"
    resources:
        memory = 1500
    log:
        "logs/MetaPhlAn2/{sample}/{sample}_reads_per_taxon_MetaPhlAn2.log"
    shell:
        '''
        set +e

        grep -v '#' {input} | grep -v '^$' | cut -f2 | sed 's/|t__.\+//g' | sort | uniq -c | \
        awk '{{print $1"\t"$2}}' | cut -f1 > {output.counts}


        exitcode=${{PIPESTATUS[1]}}
        if [[ exitcode -eq 1 ]]; then
            awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0}}' > {params.empty_ReadsTaxon_NamesRanks};
            awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0"\\tUnclassified\\t"0}}' > \
            {params.empty_Correct_Incorrect}
            touch {output.MPA2_names}
            exit 0
        else
            grep -v '#' {input} | grep -v '^$' | cut -f2 | sed 's/|t__.\+//g' | sort | uniq -c | \
            awk '{{print $1"\t"$2}}' | cut -f2 | grep -oP '\|\w__[^\|]+$' | \
            sed 's/^|\w__//g' | sed 's/_/ /g' > {output.MPA2_names}
        fi

        '''

rule get_taxID_from_MPA2_names:
    input:
        counts = ancient("MetaPhlAn2_output/{sample}/{sample}_counts.txt"),
        MPA2_names = ancient("MetaPhlAn2_output/{sample}/LastTaxRank")
    output:
        taxIDs = temp("MetaPhlAn2_output/{sample}/MPA2_taxID"),
        ReadspTaxon = "MetaPhlAn2_output/{sample}/{sample}_ReadspTaxon.txt"
    params:
        runtime = "01:00:00",
        path2Rscript = config['path2Rscript']
    resources:
        memory = 3000
    log:
        "logs/MetaPhlAn2/{sample}/{sample}_get_taxID_from_MPA2_names.log"
    shell:
        '''
        module add R/3.6.1;

        Rscript --vanilla {params.path2Rscript} -i {input.MPA2_names} -o {output.taxIDs}

        paste {input.counts} {output.taxIDs} > {output.ReadspTaxon}
        '''

rule get_name_rank:
    input:
        taxon_id = ancient("{classifier}_output/{sample}/{sample}_ReadspTaxon.txt")
    output:
        name_rank = temp("{classifier}_output/{sample}/{sample}_NamesRanks.tsv"),
        join = temp("{classifier}_output/{sample}/{sample}_joined.tsv")
    params:
        api_key = config['api_key'],
        runtime = "04:00:00"
    resources:
        memory = 2000
    log:
        "logs/{classifier}/{sample}/{sample}_get_name_rank.log"
    shell:
        '''
        for i in `cut -f2 {input.taxon_id}`; do
        efetch -db taxonomy -id ${{i}} -format xml | \
        xtract -pattern Taxon -element TaxId ScientificName Rank >> {output.name_rank};
        done

        paste {input.taxon_id} {output.name_rank} > {output.join}
        '''

rule reformat_extra_taxids:
    input:
        "{classifier}_output/{sample}/{sample}_joined.tsv"
    output:
        ReadsTaxon="{classifier}_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv",
        moreTaxID="{classifier}_output/{sample}/{sample}_extra_TaxIDs_retrieved.tsv"
    params:
        runtime = "01:00:00"
    resources:
        memory = 1000
    log:
        "logs/{classifier}/{sample}/{sample}_reformat_extra_taxids.log"
    shell:
        '''
        awk 'BEGIN{{FS="\\t"; OFS="\\t"}} \
        {{if (NF>5) {{print $0 > "{output.moreTaxID}"; print $1,$2,$4,$5,$6 > \
        "{output.ReadsTaxon}"}} \
        else {{ print$0 > "{output.ReadsTaxon}"}} }}' {input}

        touch {output.moreTaxID}
        '''

# rule reads_per_taxon_MetaPhlAn2:
#     input:
#         "MetaPhlAn2_output/{sample}/{sample}_out.txt"
#     output:
#         counts = temp("MetaPhlAn2_output/{sample}/{sample}_ReadspTaxon_counts.txt"),
#         MPA2_names = temp("MetaPhlAn2_output/{sample}/{sample}_MPA2_names.txt"),
#         name_rank_id = temp("MetaPhlAn2_output/{sample}/{sample}_name_rank_id.txt"),
#         counts_name_rank_id = "MetaPhlAn2_output/{sample}/{sample}_ReadspTaxon_name_rank_id.txt"
#     params:
#         runtime = "02:00:00",
#         api_key = config['api_key']
#     resources:
#         memory = 3000
#     log:
#         "logs/DIAMOND/{sample}_reads_per_taxon_DIAMOND.log"
#     shell:
#         '''
#         grep -v '#' {input} | grep -v '^$' | cut -f2 | sort | uniq -c | \
#         awk '{{print $1"\t"$2}}' | cut -f1 > {output.counts}

#         grep -v '#' {input} | grep -v '^$' | cut -f2 | sort | uniq -c | awk '{{print $1"\t"$2}}' | \
#         cut -f2 | sed 's/|t__.\+//g' | grep -oP '\|\w__[^\|]+$' | sed 's/^|\w__//g' | \
#         sed 's/_/ /g' > {output.MPA2_names}

#         cat {output.MPA2_names} | while IFS= read -r line; do
#         esearch -api_key {params.api_key} -db taxonomy \
#         -query "${{line}}[Scientific Name]" < /dev/null | efetch -format xml | \
#         xtract -pattern Taxon -element TaxId ScientificName Rank;
#         done >> {output.name_rank_id};

#         paste {output.counts} {output.name_rank_id} > {output.counts_name_rank_id}
#         '''
