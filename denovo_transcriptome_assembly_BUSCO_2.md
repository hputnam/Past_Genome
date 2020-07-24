# Porites astreoides transcriptome assembly

## Table of Contents
  - [1. Quality & Trim](#1-Assess-quality-of-reads-using-FastQC-and-trim-using-fastq)
  - [2. *De novo* Transcriptome Assembly](#2-De-novo-transcriptome-assembly-with-Trinity)
    - [2.1 What is Trinity?](#21-What-is-Trinity)
  - [3. Assembly completeness](#3-Assembly-completeness)
    - [3.1 Methods and tools to test](#31-Methods-and-tools-to-test)
      - [BUSCO](#BUSCO)
      - [Orthofinder](#Orthofinder)
      - [Alignment](#Alignment-to-related-species-transcriptomes)
  - [4. Annotate the assembly](#4-Annotate-the-assembly)

-------

# Example pipelines for RNAseq and *de novo* transcriptome assembly

* et al. ####* <URL here>

Citation:

![assembly_pipeline_example](https://samgurr.github.io/SamJGurr_Lab_Notebook/images/assembly_pipeline_example.JPG "assembly_pipeline_example")

* et al. ####* <URL here>

Citation:

![assembly_pipeline_example_2](https://samgurr.github.io/SamJGurr_Lab_Notebook/images/assembly_pipeline_example_2.JPG "assembly_pipeline_example_2")


-------

# 1. Assess quality of reads using FastQC and trim using fastq

-------

# 2. De novo transcriptome assembly with Trinity

##  2.1 What is Trinity?

*Trinity wiki page here* https://github.com/trinityrnaseq/trinityrnaseq/wiki

*Trinity publication* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3571712/

*Trinity downstream analysis publication* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/

- Trinity is used fro de novo assembly of RNA-seq data from the illumina platform
- *de novo* in bioinformatics referes to genome assembly based on read data **without** use of a reference
- *The following is a summary of the github wiki above...*
  - assembles the RNA into full-length transcripts or dominant isoforms of genes (contigs) and reports  only the unique portion of the alternativel spliced isoforms (*Inchworm*)
  - takes the output and runs *de Bruijn* graphs to form gene clusters (*Chrysalis*)
  - processes the *de Bruijn* graphs and reports the full-length transcripts for all isoforms of each gene (*Butterfly*)

- What is a *de Bruijn* graph?

*Read about de Bruijn graphs here (image below)* https://homolog.us/Tutorials/book4/p2.1.html

- below is an image of a *de Bruijn* graph
  - this simplified example shows a contructed sequence by splitting into k-mers (k=7 or 7-mers) that connect to form the putative sequence assembly with an overlap of 6 nt or k-1

![debruijn_graph](https://samgurr.github.io/SamJGurr_Lab_Notebook/images/debruijn_graph.JPG "debruijn_graph")

#### Trinity shell script  *(not tested!)*:

##### PHASE 1 Trinity

      #!/bin/bash
      #SBATCH -t 30:00:00
      #SBATCH --nodes=1 --ntasks-per-node=1
      #SBATCH --export=NONE
      #SBATCH --account=putnamlab
      #SBATCH --output=trinity_out
      #SBATCH --job-name=trinity-phase1
      #SBATCH --mem=220G             # maximum memory available to Trinity
      #SBATCH --partition=bigmem    # based on memory requirements
      #SBATCH --hint=nomultithread  # disable hyper-threading
      #SBATCH -D /data/putnamlab/<your folder>

      module load <trininty package here: i.e. Trinity/2.8.4-foss-2016b>

      # run trinity, stop before phase 2
      srun Trinity --seqType <insert file type; i.e. fq> --left <path to forward.fq> --right <reverse.fq> --CPU 6 --max_memory 220G

# 3. Assembly completeness

##### Why assess assembly 'completeness' in your workflow?

- important to verify the quality of the assembly (i.e. gaps, inconsistancies in contig connection and direction) using ortholog markers before moving forward in downstream annotation and analysis

## 3.1 Methods and tools to test

#### BUSCO ( Benchmarking Universal Single-Copy Orthologs)
- uses highly conserved single-copy orthologs; evolutionary informed expectations of gene content.
- appears that youu can focus a BUSCO analysis to orthologs related to your target taxa.

- image below shows a BUSCO analysis comparing the crayfish targetted for tde novo transcriptome assembly to 44 other arthropod species assemblies and a single vertebrate assembly:

*Theissinger et al. 2016* https://www.sciencedirect.com/science/article/abs/pii/S1874778716300137

Citation: Theissinger, K., Falckenhayn, C., Blande, D., Toljamo, A., Gutekunst, J., Makkonen, J., ... & Kokko, H. (2016). De Novo assembly and annotation of the freshwater crayfish Astacus astacus transcriptome. Marine Genomics, 28, 7-10.

![busco_analysis](https://samgurr.github.io/SamJGurr_Lab_Notebook/images/busco_analysis.JPG "busco_analysis")


#### Orthofinder
- uses all-v-all BLAST to assess completeness of your assembly

*Paper comparing BUSCO and orthofinder* https://genome.cshlp.org/content/early/2019/06/21/gr.243212.118.full.pdf

Altenhoff, A. M., Levy, J., Zarowiecki, M., Tomiczek, B., Vesztrocy, A. W., Dalquen, D. A., ... & Dessimoz, C. (2019). OMA standalone: orthology inference among public and custom genomes and transcriptomes. Genome research, 29(7), 1152-1163.

- authors found Orthofinder detected more othologous groups than other methods...

![assembly_completeness_fig](https://samgurr.github.io/SamJGurr_Lab_Notebook/images/assembly_completeness_fig.JPG "assembly_completeness_fig")


#### Alignment to related species transcriptome(s)
 - Publicly available transcriptomes available on NCBI:
      - *Porites astreoides*
      -  Related taxa



-------

# 4. Annotate the assembly
