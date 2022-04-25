# BUSCO - Benchmarking Universal Single-Copy Orthologs
- Comparison of assembly againsted expected gene sets  

## Resoucces
[BUSCO citation](https://academic.oup.com/bioinformatics/article/31/19/3210/211866)  
[BUSCO homepage](https://busco.ezlab.org/)  
[BUSCO user manual](https://busco.ezlab.org/busco_userguide.html) 

## Results

```shell
# BUSCO version is: 4.0.6
# The lineage dataset is: metazoa_odb10 (Creation date: 2019-11-20, number of species: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/REFS/Past/Past_genome_filtered_v1_Genewiz.fasta
# BUSCO was run in mode: genome

        ***** Results: *****

        C:69.3%[S:60.8%,D:8.5%],F:4.6%,M:26.1%,n:954
        661     Complete BUSCOs (C)
        580     Complete and single-copy BUSCOs (S)
        81      Complete and duplicated BUSCOs (D)
        44      Fragmented BUSCOs (F)
        249     Missing BUSCOs (M)
        954     Total BUSCO groups searched
```

## Run BUSCO on bluewaves
1. All BUSCO related files and scripts are in /data/putnamlab/shared/busco/
2. Metazoa_odb10 database used for corals and other files it requires to run are downloaded and stored here.
3. Follow directions in readme.txt that houses all the instructions to run it. Pasted below: 

```shell

How to run busco
----------------

sbatch -o ~/%u-%x.%j.out -e ~/%u-%x.%j.err /data/putnamlab/shared/busco/scripts/run-busco.sh

This will run /data/putnamlab/REFS/Past/Past_genome_filtered_v1_Genewiz.fasta as the query
and /data/putnamlab/shared/busco/downloads/lineages/metazoa_odb10 as the db to compare

The log of this will be written to /data/putnamlab/<user>/<user>-busco-<jobid>.out and .err

Passing different query/db_to_compare
-------------------------------------

To use a different query (genome/transcriptome) or db to compare do
# Different query: Use --export and the FULL path of the query

sbatch -o ~/%u-%x.%j.out -e ~/%u-%x.%j.err \
       --export query=<FULL PATH_to_your_query.fa>  \
       /data/putnamlab/shared/busco/scripts/run-busco.sh

# this will compare <your_query.fa> against
# /data/putnamlab/shared/busco/shared/downloads/lineages/metazoa_odb10

# Different db to compare: Use --export and FULL path to your db
# Note: make sure the db you want to use is untar'd. Currently the only db that is available is metazoa_odb10                                                      # If you need to use another db the location for all dbs is: /data/putnamlab/shared/busco/shared/downloads/lineages/

sbatch -o ~/%u-%x.%j.out -e ~/%u-%x.%j.err \
       --export query=<FULL_PATH_to_your.fa>,db_to_compare=/data/putnamlab/shared/busco/shared/downloads/lineages/insecta_odb10 \
       /data/putnamlab/shared/busco/scripts/run-busco.sh

```

