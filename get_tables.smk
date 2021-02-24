# Get tables with mean of the stats (proportions, sensitivity and precision)

configfile: 'config_get_tables.yaml'

rule all:
    input:
        expand("{condition}_{value}/stats_means_classifiers.csv", condition=config['condition'], 
            value=config['condition_values'])

rule run_Classifier_comparison:
    input: 
        expand("{condition}_{value}/{classifier}_120v_RawNumbers.txt", 
            value="{value}", condition=config['condition'], classifier=config['classifiers']), 
        expand("{condition}_{value}/{classifier}_output/Reported_virus_{classifier}.txt", 
            value="{value}", condition=config['condition'], classifier=config['classifiers'])
    output:
        "{condition}_{value}/stats_means_classifiers.csv"
    params:
        runtime = "01:00:00",
        script = config['script_Classifier_comparison'],
        path = "{condition}_{value}/",
        classifiers = config['classifiers_R']
    resources:
        memory = 2000
    log:
        "logs/{condition}_{value}/{condition}_{value}.log"
    shell:
        '''
        module add R/3.6.1;

        Rscript --vanilla {params.script} -p {params.path} -c {params.classifiers}
        '''
