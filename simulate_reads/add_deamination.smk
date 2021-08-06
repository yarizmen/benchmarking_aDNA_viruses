# Add deamination damage to reads

configfile: 'config_add_deamination.yaml'

rule all:
    input:
        expand("samples/fasta_files/{fasta}.fna", fasta=config['source_fastas'])

rule run_gargammel:
    input:
        "non_deaminated_reads/{fasta}.fna"
    output:
        "samples/fasta_files/{fasta}.fna"
    params:
        nick = 0.03,
        lovends = 0.25,
        deamd = 0.01,
        deams = config['deams']
    resources:
        memory = 8000
    log:
        logs/gargammel/{fasta}/{fasta}_run_gargammel.log
    shell:
        '''
        deamSim -damage {params.nick},{params.lovends},{params.deamd},{params.deams} {input} > \
        {output}
        '''