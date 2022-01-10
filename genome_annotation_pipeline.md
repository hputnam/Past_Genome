# *Porites astreoides* Genome Annotation

## </span>**Table of Contents**</span>

### Structural Annotation

*Initial MAKER Analysis*

1. [MAKER Round 1](#1.-MAKER-Round-1)

*Training gene prediction software*

2. [SNAP Round 1](#2.-SNAP-Round-1)

3. [AUGUSTUS Round 1](#3.-AUGUSTUS-Round-1)

*MAKER with Ab Initio gene predictors*

4. [MAKER Round 2](#4.-MAKER-Round-2)

*Iteratively running MAKER to improve annotation*

5. [SNAP Round 2](#4.-SNAP-Round-2)

6. [AUGUSTUS Round 2](#6.-AUGUSTUS-Round-2)

7. [MAKER Round 3](#7.-MAKER-Round-3)

8. [BUSCO](#8.-BUSCO)


### Functional Annotation



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
- [darencard tutorial](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2)

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


```bash
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

```bash
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

gff3_merge -s -d Rnd1_master_datastore_index.log > Rnd1.all.gff

#GFF w/o the sequences

gff3_merge -n -s -d Rnd1_master_datastore_index.log > Rnd1.all.noseq.gff

maker2zff -l 50 -x 0.5 Rnd1.all.gff

echo "Mission complete." $(date)
```

## 2. SNAP Round 1

https://github.com/KorfLab/SNAP

`mkdir snap_1`

`nano snap_1.sh`

```bash
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

## 3. AUGUSTUS Round 1

Here, I followed the AUGUSTUS training from the [darencard tutorial](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2). We are training AUGUSTUS using BUSCO.

First, we have to create training sequences from the MAKER round 1 output.

`nano transcripts1000.sh`

```bash
#!/bin/bash
#SBATCH --job-name="trans1000"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/aug_training_1
#SBATCH --mem=100GB

module load BEDTools/2.30.0-GCC-10.2.0

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' ../maker_rnd1.1/Rnd1.maker.output/Rnd1.all.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi ../past_filtered_assembly.fasta -bed - -fo Rnd1.all.maker.transcripts1000.fasta

```

### Setting up AUGUSTUS config

Since I do not have permissions to write to the server AUGUSTUS config path, I have to copy the config that to my own working directory, then export the config path in my script.

```bash
cd /opt/software/AUGUSTUS/3.4.0-foss-2020b/

cp -r config /data/putnamlab/kevin_wong1/Past_Genome/aug_training_1/
```
### Running AUGUSTUS training script

`nano aug_training1.sh`

```bash
#!/bin/bash
#SBATCH --job-name="aug_training1"
#SBATCH -t 150:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/aug_training_1
#SBATCH --mem=500GB
#SBATCH --exclusive

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#locating AUGUSTUS config path
export AUGUSTUS_CONFIG_PATH=/data/putnamlab/kevin_wong1/Past_Genome/aug_training_1/config

#run BUSCO
busco \
--config config.ini \
--in Rnd1.all.maker.transcripts1000.fasta \
--out past_rnd1_maker \
-l busco_downloads/metazoa_odb10 \
-m genome \
-f \
--long \
--augustus \
--augustus_parameters='--progress=true' \
--offline

echo "BUSCO Mission complete!" $(date)
```

`sbatch /data/putnamlab/kevin_wong1/Past_Genome/aug_training_1/aug_training1.sh`

### Renaming and moving AUGUSTUS training outputs into the species folder

The directory where the AUGUSTUS outputs are:

`ls /data/putnamlab/kevin_wong1/Past_Genome/aug_training_1/past_rnd1_maker/run_metazoa_odb10/augustus_output/retraining_parameters/BUSCO_past_rnd1_maker`

```bash
BUSCO_past_rnd1_maker_exon_probs.pbl    
BUSCO_past_rnd1_maker_intron_probs.pbl  
BUSCO_past_rnd1_maker_metapars.cgp.cfg  
BUSCO_past_rnd1_maker_parameters.cfg        
BUSCO_past_rnd1_maker_weightmatrix.txt
BUSCO_past_rnd1_maker_igenic_probs.pbl  
BUSCO_past_rnd1_maker_metapars.cfg      
BUSCO_past_rnd1_maker_metapars.utr.cfg  
BUSCO_past_rnd1_maker_parameters.cfg.orig1
```

The files need to be in a similar format as this in the AUGUSTUS species config path:

```bash
nematostella_vectensis_exon_probs.pbl    
nematostella_vectensis_intron_probs.pbl  
nematostella_vectensis_weightmatrix.txt
nematostella_vectensis_igenic_probs.pbl  
nematostella_vectensis_parameters.cfg    
README.TXT
```

Making a new species folder in the AUGUSTUS config path:

`mkdir /data/putnamlab/kevin_wong1/Past_Genome/aug_training_1/config/species/porites_astreoides1`

Copying all files over to the new species folder:

`cp ../../../past_rnd1_maker/run_metazoa_odb10/augustus_output/retraining_parameters/BUSCO_past_rnd1_maker/B* .`

Renaming files in the new species folder:

``` bash
mv BUSCO_past_rnd1_maker_exon_probs.pbl porites_astreoides1_exon_probs.pbl  
mv BUSCO_past_rnd1_maker_intron_probs.pbl porites_astreoides1_intron_probs.pbl
mv BUSCO_past_rnd1_maker_metapars.cgp.cfg porites_astreoides1_metapars.cgp.cfg
mv BUSCO_past_rnd1_maker_parameters.cfg porites_astreoides1_parameters.cfg       
mv BUSCO_past_rnd1_maker_weightmatrix.txt porites_astreoides1_weightmatrix.txt
mv BUSCO_past_rnd1_maker_igenic_probs.pbl porites_astreoides1_igenic_probs.pbl
mv BUSCO_past_rnd1_maker_metapars.cfg porites_astreoides1_metapars.cfg      
mv BUSCO_past_rnd1_maker_metapars.utr.cfg porites_astreoides1_metapars.utr.cfg
mv BUSCO_past_rnd1_maker_parameters.cfg.orig1 porites_astreoides1_parameters.cfg.orig1
```

Also coping in the original files with names because there seems to be an issue with MAKER finding these files.

`cp ../../../past_rnd1_maker/run_metazoa_odb10/augustus_output/retraining_parameters/BUSCO_past_rnd1_maker/B* .`

Now these files can be accessed by MAKER

## 4. MAKER Round 2

### Need to create gff files from the first round of MAKER to input into the second round.

`/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd1.1/Rnd1.maker.output`

```bash
# transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' Rnd1.all.noseq.gff > Rnd1.all.maker.est2genome.gff
# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' Rnd1.all.noseq.gff > Rnd1.all.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' Rnd1.all.noseq.gff > Rnd1.all.maker.repeats.gff
```

### Making a new maker folder with the opts files

```bash
mkdir ../../maker_rnd2.1

module load maker/3.01.03

maker -CTL
```

### Modifications to **maker_opts.ctl**

```bash
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
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd1.1/Rnd1.maker.output/Rnd1.all.maker.est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd1.1/Rnd1.maker.output/Rnd1.all.maker.protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd1.1/Rnd1.maker.output/Rnd1.all.maker.repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=/data/putnamlab/kevin_wong1/Past_Genome/past1.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=porites_astreoides1 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=1 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

### Adding AUGUSTUS path to the maker_exe.ctl file

`nano maker_exe.ctl`

```bash
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/opt/software/BLAST+/2.9.0-gompi-2019b/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/opt/software/BLAST+/2.9.0-gompi-2019b/bin/blastn #location of NCBI+ blastn executable
blastx=/opt/software/BLAST+/2.9.0-gompi-2019b/bin/blastx #location of NCBI+ blastx executable
tblastx=/opt/software/BLAST+/2.9.0-gompi-2019b/bin/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
prerapsearch= #location of prerapsearch executable
rapsearch= #location of rapsearch executable
RepeatMasker=/opt/software/RepeatMasker/4.0.9-p2-gompi-2019b-HMMER/RepeatMasker #location of RepeatMasker executable
exonerate=/opt/software/Exonerate/2.4.0-GCC-8.3.0/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/opt/software/SNAP/2013-11-29-GCC-8.3.0/bin/snap #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=/opt/software/AUGUSTUS/3.4.0-foss-2020b/bin/augustus #location of augustus executable
fgenesh= #location of fgenesh executable
evm= #location of EvidenceModeler executable
tRNAscan-SE= #location of trnascan executable
snoscan= #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)

```

### Running MAKER Round 2

`nano maker_rnd2.1.sh`

```
#!/bin/bash
#SBATCH --job-name="MAKER_RND2.1"
#SBATCH -t 500:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd2.1
#SBATCH --mem=120GB
#SBATCH --exclusive

module load maker/3.01.03-foss-2020b

#locating AUGUSTUS config path
export AUGUSTUS_CONFIG_PATH=/data/putnamlab/kevin_wong1/Past_Genome/aug_training_1/config

#Running maker with MPI
maker -cpus $SLURM_CPUS_ON_NODE -base Rnd2.1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "Mission complete." $(date)
```

#### Checking completion

`less -S Rnd2.1_master_datastore_index.log `

All contigs were completed successfully (Started and finished)

#### Merging files into one GFF3 and fasta file

`nano maker_rnd2.1_merge.sh`

```bash
#!/bin/bash
#SBATCH --job-name="MAKER_RND2.1_merge"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd2.1/Rnd2.1.maker.output
#SBATCH --mem=100GB

module load maker/3.01.03

fasta_merge -d Rnd2.1_master_datastore_index.log

gff3_merge -s -d Rnd2.1_master_datastore_index.log > Rnd2.1.all.gff

#GFF w/o the sequences

gff3_merge -n -s -d Rnd2.1_master_datastore_index.log > Rnd2.1.all.noseq.gff

maker2zff -l 50 -x 0.5 Rnd2.1.all.gff

echo "Mission complete." $(date)
```

## 5. SNAP Round 2

https://github.com/KorfLab/SNAP

`mkdir snap_2`

`nano snap_2.sh`

```bash
#!/bin/bash
#SBATCH --job-name="snap_2"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/snap_2
#SBATCH --mem=100GB

module load SNAP/2013-11-29-GCC-8.3.0

fathom -categorize 1000 ../maker_rnd2.1/Rnd2.1.maker.output/genome.ann ../maker_rnd2.1/Rnd2.1.maker.output/genome.dna

fathom -export 1000 -plus uni.ann uni.dna

forge export.ann export.dna

hmm-assembler.pl past . > ./past2.hmm

echo "Mission complete." $(date)
```

## 6. AUGUSTUS Round 2

Here, I followed the AUGUSTUS training from the [darencard tutorial](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2). We are training AUGUSTUS using BUSCO.

First, we have to create training sequences from the MAKER round 2 output.

`mkdir aug_training_2`

`nano transcripts1000.sh`

```bash
#!/bin/bash
#SBATCH --job-name="trans1000"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2
#SBATCH --mem=100GB

module load BEDTools/2.30.0-GCC-10.2.0

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' ../maker_rnd2.1/Rnd2.1.maker.output/Rnd2.1.all.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi ../past_filtered_assembly.fasta -bed - -fo Rnd2.1.all.maker.transcripts1000.fasta

```

### Setting up AUGUSTUS config

Since I do not have permissions to write to the server AUGUSTUS config path, I have to copy the config that to my own working directory, then export the config path in my script.

```bash
cd /opt/software/AUGUSTUS/3.4.0-foss-2020b/

cp -r config /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/

cp /data/putnamlab/kevin_wong1/Past_Genome/aug_training_1/config.ini /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/

cp -r /data/putnamlab/kevin_wong1/Past_Genome/aug_training_1/busco_downloads /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/
```
### Running AUGUSTUS training script

`nano aug_training2.sh`

```bash
#!/bin/bash
#SBATCH --job-name="aug_training2"
#SBATCH -t 150:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2
#SBATCH --mem=500GB
#SBATCH --exclusive

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#locating AUGUSTUS config path
export AUGUSTUS_CONFIG_PATH=/data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/config

#run BUSCO
busco \
--config config.ini \
--in Rnd2.1.all.maker.transcripts1000.fasta \
--out past_rnd2_maker \
-l busco_downloads/metazoa_odb10 \
-m genome \
-f \
--long \
--augustus \
--augustus_parameters='--progress=true' \
--offline

echo "BUSCO Mission complete!" $(date)
```

`sbatch /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/aug_training2.sh`


### Renaming and moving AUGUSTUS training outputs into the species folder

The directory where the AUGUSTUS outputs are:

`ls /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/past_rnd2_maker/run_metazoa_odb10/augustus_output/retraining_parameters/BUSCO_past_rnd2_maker`

```bash
BUSCO_past_rnd2_maker_exon_probs.pbl    
BUSCO_past_rnd2_maker_intron_probs.pbl  
BUSCO_past_rnd2_maker_metapars.cgp.cfg  
BUSCO_past_rnd2_maker_parameters.cfg        
BUSCO_past_rnd2_maker_weightmatrix.txt
BUSCO_past_rnd2_maker_igenic_probs.pbl  
BUSCO_past_rnd2_maker_metapars.cfg      
BUSCO_past_rnd2_maker_metapars.utr.cfg  
BUSCO_past_rnd2_maker_parameters.cfg.orig1
```

Making a new species folder in the AUGUSTUS config path:

`mkdir /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/config/species/porites_astreoides2`

Copying all files over to the new species folder:

`cp ../../../past_rnd2_maker/run_metazoa_odb10/augustus_output/retraining_parameters/BUSCO_past_rnd2_maker/B* .`

Renaming and copying files in the new species folder:

``` bash
cp BUSCO_past_rnd2_maker_exon_probs.pbl porites_astreoides12_exon_probs.pbl  
cp BUSCO_past_rnd2_maker_intron_probs.pbl porites_astreoides2_intron_probs.pbl
cp BUSCO_past_rnd2_maker_metapars.cgp.cfg porites_astreoides2_metapars.cgp.cfg
cp BUSCO_past_rnd2_maker_parameters.cfg porites_astreoides2_parameters.cfg       
cp BUSCO_past_rnd2_maker_weightmatrix.txt porites_astreoides2_weightmatrix.txt
cp BUSCO_past_rnd2_maker_igenic_probs.pbl porites_astreoides2_igenic_probs.pbl
cp BUSCO_past_rnd2_maker_metapars.cfg porites_astreoides2_metapars.cfg      
cp BUSCO_past_rnd2_maker_metapars.utr.cfg porites_astreoides2_metapars.utr.cfg
cp BUSCO_past_rnd2_maker_parameters.cfg.orig1 porites_astreoides2_parameters.cfg.orig1
```

Now these files can be accessed by MAKER


## 7. MAKER Round 3

### Need to create gff files from the first round of MAKER to input into the second round.

`/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd2.1/Rnd2.1.maker.output`

```bash
# transcript alignments
awk '{ if ($2 ~ "est2genome") print $0 }' Rnd2.1.all.noseq.gff > Rnd2.1.all.maker.est2genome.gff
# protein alignments
awk '{ if ($2 ~ "protein2genome") print $0 }' Rnd2.1.all.noseq.gff > Rnd2.1.all.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' Rnd2.1.all.noseq.gff > Rnd2.1.all.maker.repeats.gff
```

### Making a new maker folder with the opts files

```bash
mkdir ../../maker_rnd3

module load maker/3.01.03

maker -CTL
```

### Modifications to **maker_opts.ctl**

```bash
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
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd2.1/Rnd2.1.maker.output/Rnd2.1.all.maker.est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd2.1/Rnd2.1.maker.output/Rnd2.1.all.maker.protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=/data/putnamlab/kevin_wong1/Past_Genome/maker_rnd2.1/Rnd2.1.maker.output/Rnd2.1.all.maker.repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=/data/putnamlab/kevin_wong1/Past_Genome/snap_2/past2.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=porites_astreoides2 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=1 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

### Adding AUGUSTUS path to the maker_exe.ctl file

`nano maker_exe.ctl`

```bash
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/opt/software/BLAST+/2.9.0-gompi-2019b/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/opt/software/BLAST+/2.9.0-gompi-2019b/bin/blastn #location of NCBI+ blastn executable
blastx=/opt/software/BLAST+/2.9.0-gompi-2019b/bin/blastx #location of NCBI+ blastx executable
tblastx=/opt/software/BLAST+/2.9.0-gompi-2019b/bin/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
prerapsearch= #location of prerapsearch executable
rapsearch= #location of rapsearch executable
RepeatMasker=/opt/software/RepeatMasker/4.0.9-p2-gompi-2019b-HMMER/RepeatMasker #location of RepeatMasker executable
exonerate=/opt/software/Exonerate/2.4.0-GCC-8.3.0/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/opt/software/SNAP/2013-11-29-GCC-8.3.0/bin/snap #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=/opt/software/AUGUSTUS/3.4.0-foss-2020b/bin/augustus #location of augustus executable
fgenesh= #location of fgenesh executable
evm= #location of EvidenceModeler executable
tRNAscan-SE= #location of trnascan executable
snoscan= #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)

```

### Running MAKER Round 3

`nano maker_rnd3.sh`

```
#!/bin/bash
#SBATCH --job-name="MAKER_RND3"
#SBATCH -t 500:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd3
#SBATCH --mem=120GB
#SBATCH --exclusive

module load maker/3.01.03-foss-2020b

#locating AUGUSTUS config path
export AUGUSTUS_CONFIG_PATH=/data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/config

#Running maker
maker -cpus $SLURM_CPUS_ON_NODE -base Rnd3 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "Mission complete." $(date)
```

#### Checking completion

`less -S Rnd2.1_master_datastore_index.log `

All contigs were completed successfully (Started and finished)

#### Merging files into one GFF3 and fasta file

`nano maker_rnd3_merge.sh`

```bash
#!/bin/bash
#SBATCH --job-name="MAKER_RND3_merge"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd3/Rnd3.maker.output
#SBATCH --mem=100GB

module load maker/3.01.03

fasta_merge -d Rnd3_master_datastore_index.log

gff3_merge -s -d Rnd3_master_datastore_index.log > Rnd3.all.gff

#GFF w/o the sequences

gff3_merge -n -s -d Rnd3_master_datastore_index.log > Rnd3.all.noseq.gff

maker2zff -l 50 -x 0.5 Rnd3.all.gff

echo "Mission complete." $(date)
```


## Downstream processing and homology inference

`mkdir past_struc_annotations_v1`

```bash
cp Rnd3.all.gff ../../past_struc_annotations_v1/
cp Rnd3.all.noseq.gff ../../past_struc_annotations_v1/
cp Rnd3.all.maker.proteins.fasta ../../past_struc_annotations_v1/
cp Rnd3.all.maker.transcripts.fasta ../../past_struc_annotations_v1/
```

`nano maker3_postprocessing.sh`

```bash
#!/bin/bash
#SBATCH --job-name="MAKER_RND3_proc"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1
#SBATCH --mem=100GB

module load maker/3.01.03

# create naming table (there are additional options for naming beyond defaults)
maker_map_ids --prefix Pastreoides --justify 5 Rnd3.all.gff > Rnd3.all.maker.name.map

# replace names in GFF files
map_gff_ids Rnd3.all.maker.name.map Rnd3.all.gff
map_gff_ids Rnd3.all.maker.name.map Rnd3.all.noseq.gff

# replace names in FASTA headers
map_fasta_ids Rnd3.all.maker.name.map Rnd3.all.maker.transcripts.fasta
map_fasta_ids Rnd3.all.maker.name.map Rnd3.all.maker.proteins.fasta

# renaming final files

mv Rnd3.all.gff Pastreoides_all_v1.gff
mv Rnd3.all.noseq.gff Pastreoides_noseq_v1.gff
mv Rnd3.all.maker.transcripts.fasta Pastreoides_transcripts_v1.fasta
mv Rnd3.all.maker.proteins.fasta Pastreoides_protiens_v1.fasta

echo "Mission complete." $(date)
```

## Assessing Structural Annotation

### Counting number of gene models after each round

#### Round 1

```bash
cd /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd1.1/Rnd1.maker.output

cat Rnd1.all.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
```

* Gene models: 58308
* Average Gene length: 2371.32

#### Round 2

```bash
cd /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd2.1/Rnd2.1.maker.output

cat Rnd2.1.all.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
```

* Gene models: 68481
* Average Gene length: 4059.7

#### Round 3

```bash
cd /data/putnamlab/kevin_wong1/Past_Genome/maker_rnd3/Rnd3.maker.output

cat Rnd3.all.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
```

* Gene models: 64636
* Average Gene length: 4320.31

### BUSCO

#### Transcripts

`nano busco_transcripts.sh`

```
#!/bin/bash
#SBATCH --job-name="BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1
#SBATCH --mem=50GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4

#run BUSCO
busco --config /data/putnamlab/kevin_wong1/busco_downloads/config.ini \
-m transcriptome \
-i Pastreoides_transcripts_v1.fasta \
-o Past_transcripts_v1_BUSCO \
-l /data/putnamlab/kevin_wong1/busco_downloads/metazoa_odb10 \
--offline

echo "BUSCO Mission complete!" $(date)
```
