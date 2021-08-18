# *Porites astreoides* Genome Annotation

## <span style="color:blue">**Table of Contents**</span>

1. [MAKER round 1](#1.-MAKER-Round-1)
2. [MAKER round 2](#2.-MAKER-round-2)
3. [BUSCO](#3.-BUSCO)
4. [Functional Annotation](#4.-Functional-Annotation)

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

#### Shell script: maker_rnd1.sh

```
#!/bin/bash
#SBATCH --job-name="MAKER_RND1"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome
#SBATCH --mem=100GB

module load maker/3.01.03

maker -base Rnd1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "Mission complete." $(date)

```

Started: August 17, 2021 at 12pm

Ended:

## 2. MAKER Round 2

## 3. BUSCO

## 4. Functional Annotation
