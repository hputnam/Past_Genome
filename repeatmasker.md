# Repeat masking 

Repeat masking is required before annotation prediciton. This will be done using [RepeatMasker](http://www.repeatmasker.org/) in two steps. 

## Step 1: Finding repeats

RepeatMasker has a database of repeats from model organisms that one can directly make use of if working with model organisms. For non model organisms, one can find repeats and make their own custom database and use that as a reference to mask repeats in the genome fasta. It is recommended to use [RepeatModeler](https://blaxter-lab-documentation.readthedocs.io/en/latest/repeatmodeler.html) to find repeats in conjenction with RepeatMasker that can then be used to mask the repeats obtained from RepeatModeler. 

RepeatModeler is a de novo transposable element (TE) family identification and modeling package. At the heart of RepeatModeler are three de-novo repeat finding programs ( RECON, RepeatScout and LtrHarvest/Ltr_retriever ) which employ complementary computational methods for identifying repeat element boundaries and family relationships from sequence data.

RepeatModeler assists in automating the runs of the various algorithms given a genomic database, clustering redundant results, refining and classifying the families and producing a high quality library of TE families suitable for use with RepeatMasker and ultimately for submission to the Dfam database (http://dfam.org). 

Below is the script used to run Repeatmodeler. First we load RepeatModeler/1.0.8 and other required software. Then we build a database of our genome file in step 1 of the script. In the next two steps we use RepeatModeler and a script RepeatClassifier (a script within RepeatModeler) to find and classify repeats in our genome. It generates a consensus file that we then can use with RepeatMasker to mask the repeats identified here. 

```shell
#!/bin/bash

#SBATCH --job-name="repeatmasker-classifier"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G
#SBATCH --output="repeatm-%u-%x-%j"
#SBATCH --account=putnamlab
#SBATCH --export=NONE

echo "START" $(date)

#module load bio/TRF/4.09-linux64
#module load bio/RepeatModeler/1.0.8
module load Python/3.7.4-GCCcore-8.3.0

cd /data/putnamlab/tejashree/scripts/

infile=/data/putnamlab/REFS/Past/Past_genome_filtered_v1_Genewiz.fasta

RM=/data/putnamlab/tejashree/RepeatModeler

/opt/RepeatModeler/BuildDatabase -name Pastdb -engine ncbi ${infile}
$RM/RepeatModeler -engine ncbi -pa 20 -database Pastdb >& past.out
$RM/RepeatClassifier -consensi /data/putnamlab/tejashree/scripts/RM_27512.FriSep251947022020/consensi.fa >& past.classifer.out

```

## RepeatMasker 

Below is the script used to mask repeats. 

```shell

#!/bin/bash

#SBATCH --job-name="repeatmasker-classifier"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G
#SBATCH --output="repeatm-%u-%x-%j"
##SBATCH --account=putnamlab
#SBATCH --export=NONE

echo "START" $(date)

module load RepeatMasker/4.1.1-gompi-2019b-HMMER
module load Python/3.7.4-GCCcore-8.3.0

cd /data/putnamlab/tejashree/scripts/

infile=/data/putnamlab/REFS/Past/Past_genome_filtered_v1_Genewiz.fasta

RM=/data/putnamlab/tejashree/RepeatMasker

export TRF_PRGM=/data/putnamlab/tejashree/trf/bin/trf
export RMBLAST_DIR=/data/putnamlab/tejashree/rmblast/bin
$RM/RepeatMasker -engine ncbi -lib consensi.fa.classified -parallel 20 -dir . ${infile}

echo "DONE" $(date)
```
## Summary table showing percentages of types of repeats masked. 

Below is the summary file as output by RepeatMasker that shows total percentage of masked bases and the breakdown of types of repeat regions in the genome. This statistic is an important result in terms of the genome statistics that we report along with size, scaffold information etc.  

```shell
==================================================
file name: Past_genome_filtered_v1_Genewiz.fasta
sequences:          3051
total length:  677753397 bp  (677753397 bp excl N/X-runs)
GC level:         39.12 %
bases masked:  282185769 bp ( 41.64 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
SINEs:             1333       213564 bp    0.03 %
      ALUs            0            0 bp    0.00 %
      MIRs            0            0 bp    0.00 %

LINEs:            61765     26196412 bp    3.87 %
      LINE1        1004       834309 bp    0.12 %
      LINE2       22896     10809986 bp    1.59 %
      L3/CR1       3219      1234155 bp    0.18 %

LTR elements:      4722      4859871 bp    0.72 %
      ERVL            0            0 bp    0.00 %
      ERVL-MaLRs      0            0 bp    0.00 %
      ERV_classI      0            0 bp    0.00 %
      ERV_classII     0            0 bp    0.00 %

DNA elements:     29287     15555224 bp    2.30 %
     hAT-Charlie      0            0 bp    0.00 %
     TcMar-Tigger     0            0 bp    0.00 %

Unclassified:    931599    224488960 bp   33.12 %

Total interspersed repeats:271314031 bp   40.03 %


Small RNA:            0            0 bp    0.00 %
Satellites:           0            0 bp    0.00 %
Simple repeats:  131197     10089380 bp    1.49 %
Low complexity:   16915       853045 bp    0.13 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element


The query species was assumed to be homo
RepeatMasker Combined Database: Dfam_Consensus-20170127

run with rmblastn version 2.10.0+
The query was compared to classified sequences in "consensi.fa.classified"


```
