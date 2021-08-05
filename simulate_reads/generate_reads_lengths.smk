# Generate simulated reads from a given sample

configfile: 'config_generate_reads.yaml'

rule all:
    input:
        expand("samples/fastq_files/{fasta}.fq.gz", fasta=config['source_fastas'])

rule run_ART:
    input:
        "source_fasta/{fasta}.fa"
    output:
        "samples/fastq_files/{fasta}.fq"
    params:
        runtime = "10:00:00",
        seq_technology = config['seq_technology'],
        length = config['length'],
        coverage = config['coverage'],
        prefix = "samples/fastq_files/{fasta}",
        seed = config['seed']
    resources:
        memory = 8000
    log:
        "logs/ART/{fasta}/{fasta}_run_ART.log"
    shell:
        '''
        module add Phylogeny/art/2.5.8;

        art_illumina -ss {params.seq_technology} -i {input} -l {params.length} -f {params.coverage}\
        -rs {params.seed} -o {params.prefix}
        '''

rule gzip_fastq:
    input:
        "samples/fastq_files/{fasta}.fq"
    output:
        "samples/fastq_files/{fasta}.fq.gz"
    params:
        runtime = "08:00:00"
    resources:
        memory = 2000
    log:
        "logs/ART/{fasta}/{fasta}_gzip_fastq.log"
    shell:
        '''
        gzip {input}
        '''
