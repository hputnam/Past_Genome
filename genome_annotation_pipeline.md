# *Porites astreoides* Genome Annotation

## <span style="color:blue">**Table of Contents**</span>

1. [MAKER round 1](#1.-MAKER-Round-1)
2. [SNAP round 1](#2.-SNAP-Round-1)
3. [MAKER round 2](#3.-MAKER-Round-2)
4. [SNAP round 2](#4.-SNAP-Round-2)
5. [MAKER round 3](#5.-MAKER-Round-3)
6. [BUSCO](#6.-BUSCO)
7. [Functional Annotation](#7.-Functional-Annotation)

## 1. MAKER Round 1

MAKER is a genome annotation pipeline developed by the [Yandell lab](http://www.yandell-lab.org/software/maker.html), [Campbell et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4286374/).

A tutorial is available with details and explanations at [MAKER Tutorial](http://gmod.org/wiki/MAKER_Tutorial).


#### What does MAKER do?
- Identifies and masks out repeat elements
- Aligns ESTs to the genome
- Aligns proteins to the genome
- Produces ab initio gene predictions
- Synthesizes these data into final annotations
- Produces evidence-based quality values for downstream annotation management

#### Additional Information
- [Maker Information](https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Main_Page)
- [Maker Tutorial 2018](https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018)
- [Yandell tutorial](https://biohpc.cornell.edu/doc/annotation_2018_exercises1.pdf)
- [Fuller et al. 2020 A. millepora methods](https://przeworskilab.com/wp-content/uploads/Acropora_millepora_methods.pdf)

#### Files needed for MAKER round 1

To annotate a genome using Maker, you need the following files:

- The genome sequence in fasta format
- Assembled RNA-seq transcripts of the species in fasta format
- Protein sequences in fasta format, usually from closely related species or from a curated sequence database like UniProt/SwissProt.
Maker will align the transcript and protein sequences on the genome sequence to determine gene positions.

**Genome file:**
- REFERENCE: *Porites astreoides* genome assembly from GeneWiz [(link)](https://github.com/hputnam/Past_Genome/blob/master/De-novo_genome_30-323686303_GENEWIZ_Bioinformatics_Report.pdf)
- PATH: /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta

**Transcriptome file:**
- REFERENCE: *Porites astreoides* transcriptome (Florida Keys) from [Kenkel et al. 2013](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.12390) [(FASTA link)](https://matzlab.weebly.com/data--code.html)
- PATH: /data/putnamlab/kevin_wong1/Past_Genome/refs/Kenkel2013_past_transcriptome.fasta

**Protein file:**
- REFERENCE: *Porites lutea* from [Robbins et al. 2019](https://www.nature.com/articles/s41564-019-0532-4) [(FASTA link)](http://plut.reefgenomics.org/download/)
- PATH: /data/putnamlab/kevin_wong1/Past_Genome/refs/plut2v1.1.proteins.fasta

#### MAKER control files
The first thing is to create new control files for your MAKER run. The control files tell MAKER how to run and where to find additional software such as different gene callers.

In the desired folder for your MAKER run type:

`maker -CTL`

This will create three files:

- **maker_bopts.ctl** containing settings for BLAST and Exonerate.
- **masker_exe.ctl** with all the paths to different executables used by MAKER on your system.
- **maker_opts.ctl** is the file controlling MAKERs running behavior.

These files contain paths to files that are used by `maker` and our choices for analysis. Edit the maker_opts.ctl file to specify the genome assembly sequence, experimental alignment evidence and which gene finding method to use.

#### Modifications to `maker_opts.ctl`

For MAKER to run, modify the following with the appropriate paths:
- genome=**PATH_TO_GENOME**
- est=**PATH_TO_TRANSCRIPTOME**
- protein=**PATH_TO_PROTEIN**
- maxdnalength= (optional, default is 100000)


```
#-----Genome (these are always required)
genome=/data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=/data/putnamlab/kevin_wong1/Past_Genome/refs/Kenkel2013_past_transcriptome.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=/data/putnamlab/kevin_wong1/Past_Genome/refs/plut2v1.1.proteins.fasta  #protein sequence file in fasta format (i.e. from mutiple or$
protein_gff=  #aligned protein homology evidence from an external GFF3 file

maxdnalength=300000 #previously 1000000
```

`nano maker_rnd1.1.sh`

```
#!/bin/bash
#SBATCH --job-name="MAKER_RND1"
#SBATCH -t 500:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd1.1
#SBATCH --mem=250GB

module load maker/3.01.03

maker -cpus $SLURM_CPUS_ON_NODE -base Rnd1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "Mission complete." $(date)

```

*Took about 2 weeks*


#### Checking completion

`less -S Rnd1_master_datastore_index.log `

All contigs were completed successfully (Started and finished)

#### Merging files into one GFF3 and fasta file

`nano maker_rnd1_merge.sh`

```
#!/bin/bash
#SBATCH --job-name="MAKER_RND1.1_merge"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd1.1/Rnd1.maker.output
#SBATCH --mem=100GB

module load maker/3.01.03

fasta_merge -d Rnd1_master_datastore_index.log

gff3_merge -d Rnd1_master_datastore_index.log

maker2zff -l 50 -x 0.5 Rnd1.all.gff

echo "Mission complete." $(date)
```

## 2. SNAP Round 1

https://github.com/KorfLab/SNAP

`mkdir snap_1`

`nano snap_1.sh`

```
#!/bin/bash
#SBATCH --job-name="snap_1"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/snap_1
#SBATCH --mem=100GB

module load SNAP/2013-11-29-GCC-8.3.0

fathom -categorize 1000 ../maker_rnd1.1/Rnd1.maker.output/genome.ann ../maker_rnd1.1/Rnd1.maker.output/genome.dna

fathom -export 1000 -plus uni.ann uni.dna

forge export.ann export.dna

hmm-assembler.pl past . > ../past1.hmm

echo "Mission complete." $(date)
```

## 3. MAKER Round 2

`mkdir maker_rnd2`

`module load maker/3.01.03`

`maker -CTL`

Modifications to **maker_opts.ctl**

```
#-----Genome (these are always required)
genome=/data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)

#-----Re-annotation Using MAKER Derived GFF3
maker_gff=/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd1.1/Rnd1.maker.output/Rnd1.all.gff #MAKER derived GFF3 file
est_pass=1 #use ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=1 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=1 #use repeats in maker_gff: 1 = yes, 0 = no

#-----Gene Prediction
snaphmm=/data/putnamlab/kevin_wong1/Past_Genome/past1.hmm #SNAP HMM file
augustus_species=Pastreoides #Augustus gene prediction species model

pred_stats=1 #report AED and QI statistics for all predictions as well as models

```

`nano maker_rnd2.sh`

```
#!/bin/bash
#SBATCH --job-name="MAKER_RND2"
#SBATCH -t 500:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd2
#SBATCH --mem=250GB

module load maker/3.01.03

maker -cpus $SLURM_CPUS_ON_NODE -base Rnd2 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "Mission complete." $(date)

```

## 4. SNAP Round 2

## 5. MAKER Round 3

## 6. BUSCO

## 7. Functional Annotation
