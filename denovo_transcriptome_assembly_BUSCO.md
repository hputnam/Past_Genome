# Porites astreoides transcriptome assembly

## Table of Contents
[1. Quality & Trim](#1-Assess-quality-of-reads-using-FastQC-and-trim-using-fastp)

[2. *De novo* Transcriptome Assembly with Trinity](#2-De-novo-transcriptome-assembly-with-Trinity)

[3. Assembly completeness with BUSCO](#3-Assembly-completeness)

[4. Separate Host and Symbiont with psymtrans](#4-Separate-Host-and-Symbiont)

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

For small jobs in Bluewave, turn on "interactive mode"

```
interactive
```

Create a seccondary md5 file with uploaded files to compare to original md5

```
md5sum *.fastq.gz > Past_20200729.md5
```

Use 'cksum' to compare md5 files

```
cksum Past.md5 Past_20200729.md5
```

The output should look like this:
```
1158401968 432 Past.md5
282010279 432 Past2.md5
```

The first number is the cksum ID and the second number is the size of the file (in bytes).
Since both files are the same, this suggests the transfer of files did not corrupt any of the files.

```
less Past.md5
```

```
c9251a9b98768edd205724ed12c1c2a8  Sample2_R1.fastq.gz
21139f87ce4d0f4b293059a0020122e1  Sample2_R2.fastq.gz
173692a008fef82c0aecd55890eb985c  Sample3_R1.fastq.gz
962cf69dfcd1384c7f7ce3a792454ab6  Sample3_R2.fastq.gz
adf57023019e4c07fc8e5627db0faad7  Sample4_R1.fastq.gz
680ac6e4b9c4fead8f2f7a195b14dd79  Sample4_R2.fastq.gz
63c622d7c0e7093b537d0f3efd3c48c4  Sample5_R1.fastq.gz
4df8c11d45958e9081aead7d40a9386a  Sample5_R2.fastq.gz
```

```
less Past_20200729.md5
```

```
c9251a9b98768edd205724ed12c1c2a8  Sample2_R1.fastq.gz
21139f87ce4d0f4b293059a0020122e1  Sample2_R2.fastq.gz
173692a008fef82c0aecd55890eb985c  Sample3_R1.fastq.gz
962cf69dfcd1384c7f7ce3a792454ab6  Sample3_R2.fastq.gz
adf57023019e4c07fc8e5627db0faad7  Sample4_R1.fastq.gz
680ac6e4b9c4fead8f2f7a195b14dd79  Sample4_R2.fastq.gz
63c622d7c0e7093b537d0f3efd3c48c4  Sample5_R1.fastq.gz
4df8c11d45958e9081aead7d40a9386a  Sample5_R2.fastq.gz
```

Looks the same.

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

#SBATCH --job-name="fast_multi_qc_raw"
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/raw_reads

module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-2.7.15

fastqc *fastq.gz .
echo -n "Finished FastQC:"

multiqc .
echo -n "Finished multiqc:"

```


## 1.4 Quality control using fastp

https://github.com/OpenGene/fastp


```nano fastp.sh```

# loop fastp through left and right sequences; used the tutorial here: https://bash.programmingpedia.net/en/knowledge-base/11215088/bash-shell-script-two-variables-in-for-loop

```
#!/bin/bash

#SBATCH --job-name="fastp_multiqc"
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --output=../cleaned_reads/"%x_out.%j"
#SBATCH --error=../cleaned_reads/"%x_err.%j"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/raw_reads

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-2.7.15

# Make an array of sequences to trim
fwd_array=($(ls *R1.fastq.gz))
rev_array=($(ls *R2.fastq.gz))

#fastp loop
for ((i = 0; i < ${#fwd_array[@]} && i < ${#rev_array[@]}; i++)); do
         fastp --in1 ${fwd_array[i]} --in2 ${rev_array[i]} --out1 ../cleaned_reads/clean.${fwd_array[i]} --out2 ../cleaned_reads/clean.${rev_array[i]} --cut_front 20 --cut_tail 20 --cut_window_size 5 --trim_front1 3 cut_mean_quality 30 -q 30 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --json ../cleaned_reads/fastp.json  --html ../cleaned_reads/fastp.html
done

echo "Finished fastp. " $(date)

# Quality Assessment of Trimmed Reads

cd ../cleaned_reads

fastqc *fastq.gz .

echo "Finished FastQC. " $(date)

multiqc . #Compile MultiQC report from FastQC files

echo "Cleaned MultiQC report generated." $(date)

```

Count the reads before and after trimming to compare the reduction in size

### Raw reads:
#### original raw reads from /data/putnamlab/KITT/hputnam/20191010_Past_ubertrans
```
zgrep -c "@C" *.fastq.gz
```
Sample2_R1.fastq.gz:9869330
Sample2_R2.fastq.gz:9869330
Sample3_R1.fastq.gz:8471009
Sample3_R2.fastq.gz:8471009
Sample4_R1.fastq.gz:8244565
Sample4_R2.fastq.gz:8244565
Sample5_R1.fastq.gz:8920953
Sample5_R2.fastq.gz:8920953

### Cleaned reads:
#### titled 'clean' after the fastp job above
```
zgrep -c "@C" *.fastq.cleaned.gz
```
clean.Sample2_R1.fastq.gz:8726675
clean.Sample2_R2.fastq.gz:8726675
clean.Sample3_R1.fastq.gz:7255348
clean.Sample3_R2.fastq.gz:7255348
clean.Sample4_R1.fastq.gz:7237699
clean.Sample4_R2.fastq.gz:7237699
clean.Sample5_R1.fastq.gz:7673755
clean.Sample5_R2.fastq.gz:7673755


# Concatenate all R1 and R2 files together

cat *R1.fastq.gz > clean.Sample_ALL_R1.fastq.gz

cat *R2.fastq.gz > clean.Sample_ALL_R2.fastq.gz

# Summary table of counts

| Sample # | Read | Raw counts | Cleaned counts |
|:--------:|:----:|:----------:|:--------------:|
|  Sample2 |  R1  |   2394838  |     2095768    |
|  Sample2 |  R2  |   2450560  |     2022997    |
|  Sample3 |  R1  |   2076974  |     1659411    |
|  Sample3 |  R2  |   2184301  |     1784818    |
|  Sample4 |  R1  |   1981144  |     1711541    |
|  Sample4 |  R2  |   2050661  |     1705137    |
|  Sample5 |  R1  |   2124524  |     1686454    |
|  Sample5 |  R2  |   2259048  |     1902057    |
|  All_R1  |  R1  |            |     7153174    |
|  All_R2  |  R2  |            |     7415009    |

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


#### Trinity shell script:

```
nano trinity.sh
```

***note***: run this on andromeda

```
#!/bin/bash

#SBATCH --job-name="Trinity_clean_P.ast"
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH --mem=500GB
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/trinity_2

module load Trinity/2.9.1-foss-2019b-Python-3.7.4

Trinity --seqType fq --left ../cleaned_reads/clean.Sample_ALL_R1.fastq.gz  --right ../cleaned_reads/clean.Sample_ALL_R2.fastq.gz --CPU 20 --max_memory 500G --full_cleanup

```

### Trinity Statistics

```
nano trinity_stats.sh
```

```
#!/bin/bash

#SBATCH --job-name="Trinity_stats"
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH --mem=500GB
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/trinity_2

module load Trinity/2.9.1-foss-2019b-Python-3.7.4

perl /opt/software/Trinity/2.9.1-foss-2019b-Python-3.7.4/trinityrnaseq-v2.9.1/util/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity_stats.txt

```

```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  431885
Total trinity transcripts:	687749
Percent GC: 39.50

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 3573
        Contig N20: 2432
        Contig N30: 1824
        Contig N40: 1410
        Contig N50: 1094

        Median contig length: 454
        Average contig: 742.94
        Total assembled bases: 510954452


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 3394
        Contig N20: 2341
        Contig N30: 1758
        Contig N40: 1357
        Contig N50: 1049

        Median contig length: 433
        Average contig: 713.17
        Total assembled bases: 308007563

```


# 3. Assembly completeness

#### BUSCO ( Benchmarking Universal Single-Copy Orthologs)

- *Commands and overview for running BUSCO here*: https://busco.ezlab.org/busco_userguide.html
- uses highly conserved single-copy orthologs; evolutionary informed expectations of gene content.
- appears that youu can focus a BUSCO analysis to orthologs related to your target taxa.

- image below shows a BUSCO analysis comparing the crayfish targetted for tde novo transcriptome assembly to 44 other arthropod species assemblies and a single vertebrate assembly:

*Theissinger et al. 2016* https://www.sciencedirect.com/science/article/abs/pii/S1874778716300137

Citation: Theissinger, K., Falckenhayn, C., Blande, D., Toljamo, A., Gutekunst, J., Makkonen, J., ... & Kokko, H. (2016). De Novo assembly and annotation of the freshwater crayfish Astacus astacus transcriptome. Marine Genomics, 28, 7-10.


***note***: run this on bluewaves

```
sbatch -o ~/%u-%x.%j.out -e ~/%u-%x.%j.err --mail-type=BEGIN,END,FAIL --mail-user=kevin_wong1@uri.edu \
--export query=/data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/trinity_2/trinity_out_dir.Trinity.fasta  \
/data/putnamlab/kevin_wong1/scripts/run-busco-transcriptome.sh
```

### BUSCO Results

```

# BUSCO version is: 4.0.6
# The lineage dataset is: metazoa_odb10 (Creation date: 2019-11-20, number of species: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/trinity_2/trinity_$
# BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:21.4%[S:13.9%,D:7.5%],F:36.7%,M:41.9%,n:954
        205     Complete BUSCOs (C)
        133     Complete and single-copy BUSCOs (S)
        72	Complete and duplicated BUSCOs (D)
        350     Fragmented BUSCOs (F)
        399     Missing BUSCOs (M)
        954     Total BUSCO groups searched

```

# 4. Separate Host and Symbiont with psymtrans

https://github.com/sylvainforet/psytrans
https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/6-Symbiont-RNAseq-Analysis/a-QC-Align-Assemble/3_execute_psymtrans.sh


Input
- Host : Porites lutea (/data/putnamlab/kevin_wong1/REFS/Plutea/plut2v1.1.proteins.fasta)
- Symbiont: A4/A4a (/data/putnamlab/kevin_wong1/REFS/symA/syma_aug_37.aa.longest.fa)
  - http://sampgr.org.cn/index.php/download (Symbiodinium sp. protein)

***Note***: run on bluewaves

```
  mkdir tempDir
```

```
  psymtrans.sh
```

  ```
  #!/bin/bash
  #SBATCH --job-name="psytrans"
  #SBATCH -t 336:00:00
  #SBATCH --export=NONE
  #SBATCH --nodes=1 --ntasks-per-node=20
  #SBATCH --exclusive
  #SBATCH --output="%x_out.%j"
  #SBATCH --error="%x_err.%j"
  #SBATCH --mail-type=BEGIN,END,FAIL
  #SBATCH --mail-user=kevin_wong1@uri.edu
  #SBATCH -D /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/psytrans
  #SBATCH --mem=500GB
  #SBATCH -q putnamlab

  module load BLAST+/2.8.1-foss-2018b
  module load LIBSVM/3.23-foss-2018b

  echo "Separating Host and Symbiont" $(date)

  python psytrans.py -A ../../REFS/Plutea/plut2v1.1.proteins.fasta -B ../../REFS/symA/syma_aug_37.aa.longest.fa -p 20 -t tempDir ../trinity_2/trinity_out_dir.Trinity.fasta

  echo "Mission complete!" $(date)
  ```

Transcript Counts:
- Holobiont: 1375498
- Host: 1310924
- Symbiont: 64574


### BUSCO on just psytrans output of species 1 (host)

  *note: run this on bluewaves

  ```
  sbatch -o ~/%u-%x.%j.out -e ~/%u-%x.%j.err --mail-type=BEGIN,END,FAIL --mail-user=kevin_wong1@uri.edu \
  --export query=/data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/psytrans/species1_trinity_out_dir.Trinity.fasta  \
  /data/putnamlab/kevin_wong1/scripts/run-busco-transcriptome.sh
  ```
