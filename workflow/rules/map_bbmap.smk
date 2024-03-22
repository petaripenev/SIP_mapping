
rule generate_bbmap_index:
    input:
        scaffold_path = "results/scaffolds/{sample}_scaffold.fa" if config["single_scaffold_file"] == 'f' \
                   else f"results/scaffolds/{config['output_prefix']}_scaffold.fa" if config["single_scaffold_file"] == 't' \
                   else SINGLE_SCAFFOLD_ERROR,
    output:
        index=directory("results/{sample}/ref_index/{sample}") if config["single_scaffold_file"] == 'f' \
                   else directory(f"results/{config['output_prefix']}_ref_index") if config["single_scaffold_file"] == 't' \
                   else SINGLE_SCAFFOLD_ERROR,
    log:
        "results/{sample}/{sample}_bbmap_index.log" if config["single_scaffold_file"] == 'f' \
        else "results/bbmap_index.log" if config["single_scaffold_file"] == 't' \
        else SINGLE_SCAFFOLD_ERROR,
    threads: config["max_threads"]
    benchmark:
        "benchmarks/index/{sample}_bbmap_index.tsv" if config["single_scaffold_file"] == 'f' \
        else "benchmarks/index/bbmap_index.tsv" if config["single_scaffold_file"] == 't' \
        else SINGLE_SCAFFOLD_ERROR,
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbmap.sh \
        reference={input.scaffold_path} \
        path={output.index} \
        threads={threads} \
        &> {log}
        """

rule map_fr_reads:
    input:
        index = "results/{sample}/ref_index/{sample}" if config["single_scaffold_file"] == 'f' \
           else "results/{output_prefix}_ref_index" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        fw_reads_path = "results/reads/{n}/{sample_map}_{n}_{isotope}_fw.fq.gz",
        rv_reads_path = "results/reads/{n}/{sample_map}_{n}_{isotope}_rv.fq.gz",
    output:
        coverage_stats = "results/maps/{n}/covstats_{sample}-v-{sample_map}_{n}_{isotope}.tsv" if config["single_scaffold_file"] == 'f' \
                    else "results/maps/{output_prefix}/{n}/covstats_{sample_map}_{n}_{isotope}.tsv" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        coverage_hist = "results/maps/{n}/covhist_{sample}-v-{sample_map}_{n}_{isotope}.tsv" if config["single_scaffold_file"] == 'f' \
                   else "results/maps/{output_prefix}/{n}/covhist_{sample_map}_{n}_{isotope}.tsv" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        mapping = "results/maps/{n}/{sample}-v-{sample_map}_{n}_{isotope}_mapped.sort.bam" if config["single_scaffold_file"] == 'f' \
             else "results/maps/{output_prefix}/{n}/{sample_map}_{n}_{isotope}_mapped.sort.bam" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    wildcard_constraints:
        output_prefix=config["output_prefix"],
        n=fraction_loop_regex,
        isotope=config["isotope_regex"],
    log:
        "logs/maps/{n}/{sample}-v-{sample_map}_{n}_{isotope}_mapped.log" if config["single_scaffold_file"] == 'f' \
        else "logs/maps/{output_prefix}/{n}/{sample_map}_{n}_{isotope}_mapped.log" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    params:
        ambiguous="random",
        minid=config["min_map_id"],
        perfectmode=config["perfectmode"],
    threads: config["max_threads"]
    benchmark:
        "benchmarks/maps/{n}/{sample}-v-{sample_map}_{n}_{isotope}_mapped.tsv" if config["single_scaffold_file"] == 'f' \
        else "benchmarks/maps/{output_prefix}/{n}/{sample_map}_{n}_{isotope}_mapped.tsv" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbmap.sh \
        in1={input.fw_reads_path} \
        in2={input.rv_reads_path} \
        path={input.index} \
        trimreaddescriptions=t \
        perfectmode={params.perfectmode} \
        ambiguous={params.ambiguous} \
        threads={threads} \
        covstats={output.coverage_stats} \
        covhist={output.coverage_hist} \
        outm=stdout.sam 2> {log} |\
        sambam > {output.mapping}
        """

rule map_unfr_reads:
    input:
        index = "results/{sample}/ref_index/{sample}" if config["single_scaffold_file"] == 'f' \
           else "results/{output_prefix}_ref_index" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        unfr_fw = "results/reads/t{inc_time}/{sample_map}_t{inc_time}_fw.fq.gz",
        unfr_rv = "results/reads/t{inc_time}/{sample_map}_t{inc_time}_rv.fq.gz",
    output:
        coverage_stats = "results/maps/t{inc_time}/covstats_{sample}-v-{sample_map}.txt" if config["single_scaffold_file"] == 'f' \
                    else "results/maps/{output_prefix}/t{inc_time}/covstats_{sample_map}.txt" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        coverage_hist = "results/maps/t{inc_time}/covhist_{sample}-v-{sample_map}.txt" if config["single_scaffold_file"] == 'f' \
                   else "results/maps/{output_prefix}/t{inc_time}/covhist_{sample_map}.txt" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
        mapping = "results/maps/t{inc_time}/{sample}-v-{sample_map}_t{inc_time}_mapped.sort.bam" if config["single_scaffold_file"] == 'f' \
             else "results/maps/{output_prefix}/t{inc_time}/{sample_map}_t{inc_time}_mapped.sort.bam" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    log:
        "logs/maps/t{inc_time}/{sample}-v-{sample_map}_t{inc_time}_mapping.log" if config["single_scaffold_file"] == 'f' \
        else "logs/maps/{output_prefix}/t{inc_time}/{sample_map}_t{inc_time}_mapping.log" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    wildcard_constraints:
        output_prefix=config["output_prefix"],
        inc_time=config["incubation_times"],
    params:
        ambiguous="random",
        minid=config["min_map_id"],
        perfectmode=config["perfectmode"],
    threads: config["max_threads"]
    benchmark:
        "benchmarks/maps/t{inc_time}/{sample}-v-{sample_map}_t{inc_time}_mapping.tsv" if config["single_scaffold_file"] == 'f' \
        else "benchmarks/maps/{output_prefix}/t{inc_time}/{sample_map}_t{inc_time}_mapping.tsv" if config["single_scaffold_file"] == 't' else SINGLE_SCAFFOLD_ERROR,
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbmap.sh \
        in1={input.unfr_fw} \
        in2={input.unfr_rv} \
        path={input.index} \
        trimreaddescriptions=t \
        perfectmode={params.perfectmode} \
        ambiguous={params.ambiguous} \
        threads={threads} \
        covstats={output.coverage_stats} \
        covhist={output.coverage_hist} \
        outm=stdout.sam 2> {log} |\
        sambam > {output.mapping}
        """
