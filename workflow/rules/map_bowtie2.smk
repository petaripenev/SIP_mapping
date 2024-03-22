
rule generate_bowtie2_index:
    input:
        scaffold_path = "results/scaffolds/{sample}_scaffold.fa" if config["single_scaffold_file"] == 'f' \
                   else f"results/scaffolds/{config['output_prefix']}_scaffold.fa" if config["single_scaffold_file"] == 't' \
                   else SINGLE_SCAFFOLD_ERROR,
    output:
        index=directory("results/{sample}/ref_index/{sample}") if config["single_scaffold_file"] == 'f' \
                   else directory(f"results/{config['output_prefix']}_ref_index") if config["single_scaffold_file"] == 't' \
                   else SINGLE_SCAFFOLD_ERROR,
    log:
        "results/{sample}/{sample}_bowtie2_index.log" if config["single_scaffold_file"] == 'f' \
        else "results/bowtie2_index.log" if config["single_scaffold_file"] == 't' \
        else SINGLE_SCAFFOLD_ERROR,
    threads: config["max_threads"]
    params:
        bt2_index_base="results/{sample}/ref_index/{sample}/{sample}_scaffold.fa" if config["single_scaffold_file"] == 'f' \
                 else f"results/{config['output_prefix']}_ref_index/scaffold_ix.fa" if config["single_scaffold_file"] == 't' \
                 else SINGLE_SCAFFOLD_ERROR,
    shell:
        """
        /shared/software/bowtie2/2.5.2/bowtie2-build \
        {input.scaffold_path} \
        {params.bt2_index_base} \
        --threads {threads} \
        --large-index \
        &> {log}
        """

rule map_fr_reads:
    input:
        index_folder = "results/{sample}/ref_index/{sample}" if config["single_scaffold_file"] == 'f' \
                   else "results/{output_prefix}_ref_index" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        fw_reads_path = "results/reads/{n}/{sample_map}_{n}_{isotope}_fw.fq.gz",
        rv_reads_path = "results/reads/{n}/{sample_map}_{n}_{isotope}_rv.fq.gz",
    output:
        mapping = "results/maps/{n}/{sample}-v-{sample_map}_{n}_{isotope}_mapped.sort.bam" if config["single_scaffold_file"] == 'f' \
             else "results/maps/{output_prefix}/{n}/{sample_map}_{n}_{isotope}_mapped.sort.bam" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    wildcard_constraints:
        output_prefix=config["output_prefix"],
        n=fraction_loop_regex,
        isotope=config["isotope_regex"],
    log:
        bowtie2 = "logs/maps/{n}/{sample}-v-{sample_map}_{n}_{isotope}_mapped.log" if config["single_scaffold_file"] == 'f' \
        else "logs/maps/{output_prefix}/{n}/{sample_map}_{n}_{isotope}_mapped.log" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        shrinksam = "logs/maps/{n}/{sample}-v-{sample_map}_{n}_{isotope}_shrinksam.log" if config["single_scaffold_file"] == 'f' \
        else "logs/maps/{output_prefix}/{n}/{sample_map}_{n}_{isotope}_shrinksam.log" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    threads: config["max_threads"]
    params:
        bt2_index_base = "results/{sample}/ref_index/{sample}/{sample}_scaffold.fa" if config["single_scaffold_file"] == 'f' \
                   else f"results/{config['output_prefix']}_ref_index/scaffold_ix.fa" if config["single_scaffold_file"] == 't' \
                   else SINGLE_SCAFFOLD_ERROR,
    benchmark:
        "benchmarks/maps/{n}/{sample}-v-{sample_map}_{n}_{isotope}_mapped.tsv" if config["single_scaffold_file"] == 'f' \
        else "benchmarks/maps/{output_prefix}/{n}/{sample_map}_{n}_{isotope}_mapped.tsv" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    # conda:
    #     "../envs/bowtie2.yaml"
    shell:
        """
        /shared/software/bowtie2/2.5.2/bowtie2 \
        -x {params.bt2_index_base} \
        -1 {input.fw_reads_path} \
        -2 {input.rv_reads_path} \
        --threads {threads} \
        --sensitive \
        2> {log.bowtie2} \
        | shrinksam -v \
        2> {log.shrinksam} \
        > {output.mapping}
        """


rule map_unfr_reads:
    input:
        index_folder = "results/{sample}/ref_index/{sample}" if config["single_scaffold_file"] == 'f' \
                   else "results/{output_prefix}_ref_index" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        unfr_fw = "results/reads/t{inc_time}/{sample_map}_t{inc_time}_fw.fq.gz",
        unfr_rv = "results/reads/t{inc_time}/{sample_map}_t{inc_time}_rv.fq.gz",
    output:
        mapping = "results/maps/t{inc_time}/{sample}-v-{sample_map}_t{inc_time}_mapped.sort.bam" if config["single_scaffold_file"] == 'f' \
             else "results/maps/{output_prefix}/t{inc_time}/{sample_map}_t{inc_time}_mapped.sort.bam" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    log:
        bowtie2="logs/maps/t{inc_time}/{sample}-v-{sample_map}_t{inc_time}_mapping.log" if config["single_scaffold_file"] == 'f' \
        else "logs/maps/{output_prefix}/t{inc_time}/{sample_map}_t{inc_time}_mapping.log" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        shrinksam="logs/maps/t{inc_time}/{sample}-v-{sample_map}_t{inc_time}_shrinksam.log" if config["single_scaffold_file"] == 'f' \
        else "logs/maps/{output_prefix}/t{inc_time}/{sample_map}_t{inc_time}_shrinksam.log" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    params:
        bt2_index_base = "results/{sample}/ref_index/{sample}/{sample}_scaffold.fa" if config["single_scaffold_file"] == 'f' \
                   else f"results/{config['output_prefix']}_ref_index/scaffold_ix.fa" if config["single_scaffold_file"] == 't' \
                   else SINGLE_SCAFFOLD_ERROR,
    threads: config["max_threads"]
    # conda:
    #     "../envs/bowtie2.yaml"
    shell:
        """
        /shared/software/bowtie2/2.5.2/bowtie2 \
        -x {params.bt2_index_base} \
        -1 {input.unfr_fw} \
        -2 {input.unfr_rv} \
        --threads {threads} \
        --sensitive \
        2> {log.bowtie2} \
        | shrinksam -v \
        2> {log.shrinksam} \
        > {output.mapping}
        """