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

## Run with
```
snakemake --use-conda --cluster-config config/cluster.yaml --default-resources partition=standard --cluster "sbatch -J WYfrmap -p {cluster.p} -o {cluster.o}" -j8 --rerun-incomplete --latency-wait 60 --cluster-cancel scancel
```