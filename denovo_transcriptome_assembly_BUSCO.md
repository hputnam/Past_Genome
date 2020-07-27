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

* [Holzer and Marz, 2019](https://watermark.silverchair.com/giz039.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAArowggK2BgkqhkiG9w0BBwagggKnMIICowIBADCCApwGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM1XrEyX8vopR0xjmZAgEQgIICbbyfMmevYUVAJ64h-hRStbmfLG5V-NwzIso9KTu4334Vzoebw6jn-4sS6I7Qc65XkIPw3W1IkCdyyuOfMpDf5_o9UBG7eGLv74VZ2arf2TMvUFL34fGF9mDJYApAGOpBKMcrf23ddaJ2VKHYEWuSwMC81OHU2awRWYc02Uo0chb3gsvdZYjfbeQ04sGHjq3JqQafiWoS8_Z7fTRo7F7D1QBfb0wgKP1AOgfPct1HIRZ--YeyFuEqUMBQUll22DfAsOm0qFv0ywczInxqrxADzEyFodCZRA60V_CctPAsNG-1vgBxq61_8ld_bPGCkCSkXFXzT0K53OfoOh9sOly1J-MIZC-GubllQoHFMk2obgTErq-C0QhL-YI-ZlwKE69vrDgDXEze20VUY1sspJLGNtxXIx4tzHOQTlTx7NLB_0e6QoCCEBaBXNI4g3M4mIa4UpqtXaCdgVUgbvNaefJe_os1B2-aoMAbtihQBGgyHARszcDnfOZB52yMWxk1a8dpSIq_EVg_eEWNsESQZQeCQVhuTxl3Ngh0Vp226mEtF7irAdZW3Snnh9oZz8Y7-cNDFWxsjb1sS_pJEEj6ViRDbEf2cOgmVq1n1jdRMDoyoxb5r8F6wL28lNIrISkMv0_hnv1UNSSQb2yEmYHJNw7EXOmDTnKe7dATNuAb_DP3rVdvTDGPgPgh0FQeRp7Ed7zgqJbGMDU8K6_OmJZUwYOiDOSOAf1RVEDUuqb-nw0oI3P1nPvDHedO39m1I3V_D8beFExAptXDAmQQRknfKJLwAUQ9ceEZVJVeUr6i78bdN-wj8mgjeTh7KK4QUwclIQ)

Citation:

![assembly_pipeline_example](https://samgurr.github.io/SamJGurr_Lab_Notebook/images/assembly_pipeline_example.JPG "assembly_pipeline_example")

* [Carruthers et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5759245/)

Citation:

![assembly_pipeline_example_2](https://samgurr.github.io/SamJGurr_Lab_Notebook/images/assembly_pipeline_example_2.JPG "assembly_pipeline_example_2")


-------

# 1. Assess quality of reads using FastQC and trim using fastp

Current files:
- Sample2_R1.fastq.gz
- Sample2_R2.fastq.gz
- Sample3_R1.fastq.gz
- Sample3_R2.fastq.gz  
- Sample4_R1.fastq.gz
- Sample4_R2.fastq.gz
- Sample5_R1.fastq.gz
- Sample5_R2.fastq.gz
- Past.md5                       

## 1.1 Checking the files downloaded correctly (cksum)

Uploading files from an external hard drive to Bluewaves
```
scp -r /Volumes/KWong_Hardrive_2/Porites_Ref_Transcriptomes/ kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kwong/P.astreoides_transcriptome_genome/raw_files/
cd data/putnamlab/kwong/P.astreoides_transcriptome_genome/raw_files/Porites_Ref_Transcriptomes/
```

For small jobs in Bluewave, turn on "interactive mode"

```
interactive
```

Create a seccondary md5 file with uploaded files to compare to original md5

```
md5sum *.fastq.gz > Past2.md5
```

Use 'cksum' to compare md5 files

```
cksum Past.md5 Past2.md5
```

The output should look like this:
```
1158401968 432 Past.md5
282010279 432 Past2.md5
```

The first number is the cksum ID and the second number is the size of the file (in bytes).
Since both files are the same, this suggests the transfer of files did not corrupt any of the files.


Don't forget to exit interactive mode!
```
exit
```


## 1.2 Quality control of raw sequencing reads (FastQC and MultiQC)
https://github.com/s-andrews/FastQC
https://github.com/ewels/MultiQC

Make a scripts folder to stay organized.
```
mkdir scripts
cd scripts
```

Making a shell script to assess the quality of the reads using FastQC and summarizing the results using FastQC.

```
nano fastqc_raw.sh
```

```
#!/bin/bash

#SBATCH --job-name="fastqc_raw_P.ast"
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kwong/P.astreoides_transcriptome_genome/raw_files/Porites_Ref_Transcriptomes/

module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-3.6.6

fastqc *fastq.gz .
echo -n "Finished FastQC:"

multiqc .
echo -n "Finished multiqc:"

```

Running the shell script.

```
sbatch fastqc_raw.sh
```

**ERROR**

* Sample 2 read 1 may be corrupt - Fastqc was not able to run.
* Also need to learn how to export into a different folder


Export the MultiQC report to view. It should be an HTML file.
```
scp -r  kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kwong/P.astreoides_transcriptome_genome/fastqc_results/multiqc_report.html /Users/kevinwong/Documents/MyProjects/Past_Genome/output/

cat */summary.txt > /Users/kevinwong/Documents/MyProjects/P.astreoides_Transcriptome_Assembly/Output/fastqc_summaries.txt
grep FAIL fastqc_summaries.txt | wc -l
```


## 1.4 Quality control using fastp

https://github.com/OpenGene/fastp


```nano fastp.sh```

```
#!/bin/bash

#SBATCH --job-name="fastp_raw_P.ast"
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kwong/P.astreoides_transcriptome_genome/raw_files/Porites_Ref_Transcriptomes/

module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-3.6.6

for filename in *.fastq.gz
> do
> fastp --in1 ${filename} --out1 ${filename}.clean \
> --cut_front 20 --cut_tail 20 --cut_window_size 5 \
> --cut_mean_quality 15 -q 15 -w 16 --trim_front1 14 \
> --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
> done
```

Count the reads before and after trimming to compare the reduction in size

Raw reads:
```
zgrep -c *.fastq.gz
```

Cleaned reads:
```
zgrep -c *.fastq.cleaned.gz
```


## 1.5 Quality control of trimmed sequencing reads (FASTQC)

Move trimmed reads to a new folder

```
mkdir cleaned_reads
mv *fastq.cleaned.gz* cleaned_reads/
```

Make a shell script to assess the quality of the trimmed reads

```
nano fastqc_clean.sh
```

```
#!/bin/bash

#SBATCH --job-name="fastqc_clean_P.ast"
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kwong/P.astreoides_transcriptome_genome/raw_files/Porites_Ref_Transcriptomes/cleaned_reads

module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-3.6.6

fastqc *fastq.clean.gz .
echo -n "Finished FastQC:"

multiqc .
echo -n "Finished multiqc:"
```

Export the MultiQC report to view. It should be an HTML file.
```
scp -r  kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kwong/P.astreoides_transcriptome_genome/cleaned_reads/multiqc_report.html /Users/kevinwong/Documents/MyProjects/Past_Genome/output/

cat */summary.txt > /Users/kevinwong/Documents/MyProjects/P.astreoides_Transcriptome_Assembly/Output/fastqc_summaries.txt
grep FAIL fastqc_summaries.txt | wc -l
```


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

- *Commands and overview for running BUSCO here*: https://busco.ezlab.org/busco_userguide.html
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
