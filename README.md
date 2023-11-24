# Map fractions and T0 reads to assemblies

## Input files

Populate [`samples.tsv`](config/samples.tsv) with the following fields:
* `sample_name`: Name of the sample
* `scaffolds_path`: path to assembly
* `t0_fw` and `t0_rv`: paths to T0 reads
* `fraction_1_16O_fw`, `fraction_1_16O_rv`, `fraction_2_16O_fw` ... : paths to 16O fraction reads

This will assume the 18O fractions can be found in the same paths as the 18O fractions, but with the 16O replaced with 18O.

## Run with
```
snakemake --use-conda --cluster-config config/cluster.yaml --default-resources partition=standard --cluster "sbatch -J WYfrmap -p {cluster.p} -o {cluster.o}" -j8 --rerun-incomplete --latency-wait 60 --cluster-cancel scancel
```