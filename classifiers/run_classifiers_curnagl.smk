# Run Centrifuge, Kraken2, DIAMOND, MetaPhlAn2, and MALT for a given sample in curnagl. 

configfile: 'config_run_classifiers_curnagl.yaml'

rule all:
    input:
        #"MALT_output/rma6/count_table.tsv",
        "Centrifuge_output/count_table.tsv",
        "Kraken2_output/count_table.tsv",
        "DIAMOND_output/count_table.tsv",
        "MetaPhlAn2_output/count_table.tsv"

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
        database = config['path2Centrifuge_db'],
        assign = 1
    resources:
        memory = 50000
    threads: 5
    benchmark:
        "benchmarks/Centrifuge/{sample}/{sample}.benchmark.txt"
    log:
        "logs/Centrifuge/{sample}/{sample}_run_Centrifuge.log"
    shell:
        '''
        centrifuge -x {params.database} -f -U {input} -p {threads} -k {params.assign} \
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
        database = config['path2Kraken2_db']
    resources:
        memory = 80000
    threads: 5
    benchmark:
        "benchmarks/Kraken2/{sample}/{sample}.benchmark.txt"
    log:
        "logs/Kraken2/{sample}/{sample}_run_Kraken2.log"
    shell:
        '''
        kraken2 --db {params.database} {input} --threads {threads} \
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
        runtime = "40:00:00",
        mode = 'blastx',
        database = config['path2DIAMOND_db']
        # seed_shape = config['seed_shape']
    resources:
        memory = 26000 
    threads: 5
    benchmark:
        "benchmarks/DIAMOND/{sample}/{sample}.benchmark.txt"
    log:
        "logs/DIAMOND/{sample}/{sample}_run_DIAMOND.log"
    shell:
        '''
        diamond {params.mode} -p {threads} -d {params.database} -q {input} \
        -o {output} -f 102
        '''
        # --shape-mask {params.seed_shape}

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
        mpa2 = config['path2MetaPhlAn2'],
        runtime = "10:00:00",
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
        {params.mpa2} {input} --bowtie2out {output.bt2o} --nproc {threads} \
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
        Rscript --vanilla {params.path2Rscript} -i {input.MPA2_names} -o {output.taxIDs}

        paste {input.counts} {output.taxIDs} > {output.ReadspTaxon}
        '''

rule join_into_counts_table:
    input:
        expand("{classifier}_output/{sample}/{sample}_ReadspTaxon.txt", 
                sample=config['samples'], classifier="{classifier}")
    output:
        "{classifier}_output/count_table.tsv"
    params:
        runtime = "01:00:00",
        path2merge_script = config['path2merge_script']
    resources:
        memory = 5000
    log:
        "logs/{classifier}_output/join_into_counts_table.log"
    shell:
        '''
        Rscript --vanilla {params.path2merge_script} -o {output}
        '''


## run for each fastq file
## the output for each set has to be in a seperate folder
rule run_malt:
    input:
        fastq = "samples/fastq_files/{sample}.fq.gz",
        db = "/work/FAC/FBM/DBC/amalaspi/popgen/sneuensc/virome/malt/index_nucl"
    output:
        "MALT_output/rma6/{sample}.rma6"
    log:
        "logs/MALT/rma6/{sample}_run_malt.log"
    params:
        runtime = "01:00:00"
    resources:
        memory = 40000
    conda:
        "/work/FAC/FBM/DBC/amalaspi/popgen/sneuensc/conda/envs/hops"
    shell:
        """
        malt-run -J-Xmx40000m -d {input.db} -i {input.fastq} -o {output} -m BlastN
        """
 
 ## run on all runs of a set (all rma6 of a set should be located in the same folder)
rule extract_malt:
    input:
        rma6 = expand("MALT_output/rma6/{sample}.rma6", sample=config['samples'])
    output:
        "MALT_output/rma6/count_table.tsv"
    params:
        runtime = "05:00:00"
    resources:
        memory = 10000
    log:
        "logs/MALT/extract_malt.log"
    conda:
        "/work/FAC/FBM/DBC/amalaspi/popgen/sneuensc/conda/envs/megan"
    shell:
        """
        bash /work/FAC/FBM/DBC/amalaspi/popgen/sneuensc/virome/malt/rma-tabuliser/rma-tabuliser -d MALT_output/rma6 -u
        sed -i 's:_simulated::g' MALT_output/rma6/count_table.tsv
        """