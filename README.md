# Map fractions and T0 reads to assemblies

## Input files

Populate [`samples.tsv`](config/samples.tsv) with the following fields:
* `sample_name`: Name of the sample
* `scaffolds_path`: path to assembly
* `fraction_1_16O_fw`, `fraction_1_18O_fw`, `fraction_2_16O_fw` ... : paths to 16O and 18O fraction forward reads, file names should end with `PE.1.fastq.gz`. Changing the isotope name can be done here, for example if you are working with 12C/13C change the headers to `fraction_1_12C_fw`, `fraction_1_13C_fw`, `fraction_2_12C_fw` ...
* `t0_fw`: paths to forward T0 reads, file names should end with `PE.1.fastq.gz`. Same renaming will be done to get the reverse reads as the fractions.

This will assume the reverse reads can be found in the same paths as the forward reads. The following substitution will be done to get the file names:

```
"PE.1.fastq.gz" -> "PE.2.fastq.gz"
```

## Configuration

The following parameters can be set in [`config.yaml`](config/config.yaml):
* `samples`: path to [`samples.tsv`](config/samples.tsv)
* `max_threads`: number of maximum threads to use. Cluster submission automatically lowers to the number of available cores, so leave the setting to 10000
* `perfectmode` : `t` or `f`. if set to `t`, the mapping will be done in perfect mode, otherwise in sensitive mode. Perfect mode is faster but less sensitive.
* `min_map_id`: minimum mapping identity to keep a read. Default is 0.99.
* `coverm_perc_id`: minimum identity to count a read towards coverage calculation. Default is 0.99.
* `single_scaffold_mode`: `t` or `f`. If set to `t`, the pipeline will assume that all entries in column `scaffolds_path` of [`samples.tsv`](config/samples.tsv) point to the same scaffold. This mode is meant to be used if you have a specific gene of interest, have isolated all of its matches within the project, and clustered them down to a representative set. If set to `f`, the pipeline will perform an all versus all mapping between the scaffolds and reads.
* `isotope_regex`: regular expression to match the isotope name in the fraction file names. Default is `"16O|18O"`. If you are working with 12C/13C, change this to `"12C|13C"`. Should match with the fraction column names in [`samples.tsv`](config/samples.tsv).
* `isotopes`: list of isotope names to use in the output. Default is `["16O", "18O"]`. If you are working with 12C/13C, change this to `["12C", "13C"]`. Should match with the fraction column names in [`samples.tsv`](config/samples.tsv).


## Run with
```
snakemake --use-conda --cluster-config config/cluster.yaml --default-resources partition=standard --cluster "sbatch -J WYfrmap -p {cluster.p} -o {cluster.o}" -j8 --rerun-incomplete --latency-wait 60 --cluster-cancel scancel
```