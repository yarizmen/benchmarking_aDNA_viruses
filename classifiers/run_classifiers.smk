# Run Centrifuge, Kraken2, DIAMOND, and MetaPhlAn2 for a given sample

configfile: 'config_run_classifiers.yaml'

rule all:
    input:
        expand("Centrifuge_output/{sample}/{sample}_classifications.txt", sample=config['samples']),
        expand("Kraken2_output/{sample}/{sample}_K2out.txt", sample=config['samples']),
        expand("DIAMOND_output/{sample}/{sample}_tax.txt", sample=config['samples']),
        expand("MetaPhlAn2_output/{sample}/{sample}_out.txt", sample=config['samples'])

rule get_fasta:
    input:
        "samples/fastq_files/{sample}.fastq"
    output:
        "samples/fasta_files/{sample}.fna"
    resources:
        memory = 5000
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
        memory = 24000
    threads: 5
    benchmark:
        "benchmarks/Centrifuge/{sample}.benchmark.txt"
    log:
        "logs/Centrifuge/{sample}_run_Centrifuge.log"
    shell:
        '''
        {params.Centrifuge} -x {params.database} -f -U {input} -p {threads} -k {params.assign} \
        --report-file {output.report_file} -S {output.classifications}
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
        memory = 36000
    threads: 5
    benchmark:
        "benchmarks/Kraken2/{sample}.benchmark.txt"
    log:
        "logs/Kraken2/{sample}_run_Centrifuge.log"
    shell:
        '''
        {params.Kraken2} --db {params.database} {input} --threads {threads} \
        --report {output.report_file} --use-names > {output.K2out}
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
        "benchmarks/DIAMOND/{sample}.benchmark.txt"
    log:
        "logs/DIAMOND/{sample}_run_DIAMOND.log"
    shell:
        '''
        {params.DIAMOND} {params.mode} -p {threads} -d {params.database} -q {input} -o {output} \
        -f 102
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
        "benchmarks/MetaPhlAn2/{sample}.benchmark.txt"
    log:
        "logs/MetaPhlAn2/{sample}_run_MetaPhlAn2.log"
    shell:
        '''
        module add UHTS/Aligner/bowtie2/2.3.4.1;

        {params.MetaPhlAn2} {input} --bowtie2out {output.bt2o} --nproc {threads} \
        --input_type {params.i_type} -t {params.o_type} -o {output.MPA2_out}
        '''