# *Porites astreoides* Genome Annotation

## </span>**Table of Contents**</span>

### Structural Annotation

*Initial MAKER Analysis*

1. [MAKER Round 1](https://github.com/hputnam/Past_Genome/blob/master/genome_annotation_pipeline.md#1-maker-round-1)

*Training gene prediction software*

2. [SNAP Round 1](https://github.com/hputnam/Past_Genome/blob/master/genome_annotation_pipeline.md#2-snap-round-1)

3. [AUGUSTUS Round 1](https://github.com/hputnam/Past_Genome/blob/master/genome_annotation_pipeline.md#3-augustus-round-1)

*MAKER with Ab Initio gene predictors*

4. [MAKER Round 2](https://github.com/hputnam/Past_Genome/blob/master/genome_annotation_pipeline.md#4-maker-round-2)

*Iteratively running MAKER to improve annotation*

5. [SNAP Round 2](https://github.com/hputnam/Past_Genome/blob/master/genome_annotation_pipeline.md#5-snap-round-2)

6. [AUGUSTUS Round 2](https://github.com/hputnam/Past_Genome/blob/master/genome_annotation_pipeline.md#6-augustus-round-2)

7. [MAKER Round 3](https://github.com/hputnam/Past_Genome/blob/master/genome_annotation_pipeline.md#7-maker-round-3)

*Assessing completeness*

8. [BUSCO](#8.-BUSCO)


### Functional Annotation

1. BLAST the protein sequences against Swiss-Prot


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
mv Rnd3.all.maker.proteins.fasta Pastreoides_proteins_v1.fasta

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

```bash
cat Pastreoides_all_v1.gff | awk '{ if ($3 == "mRNA") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
```

* Number of mRNA Transcripts: 64636
* Average transcript length: 4320.31

```bash
grep -c ">" Pastreoides_transcripts_v1.fasta
```

* Number of transcripts in fasta file: 64636

```bash
grep -c ">" Pastreoides_proteins_v1.fasta
```

* Number of proteins in fasta file: 64636


## 8. BUSCO

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

BUSCO Output
```
--------------------------------------------------
|Results from dataset metazoa_odb10               |
--------------------------------------------------
|C:78.1%[S:69.2%,D:8.9%],F:10.9%,M:11.0%,n:954    |
|745    Complete BUSCOs (C)                       |
|660    Complete and single-copy BUSCOs (S)       |
|85     Complete and duplicated BUSCOs (D)        |
|104    Fragmented BUSCOs (F)                     |
|105    Missing BUSCOs (M)                        |
|954    Total BUSCO groups searched               |
--------------------------------------------------
```

```bash
#!/bin/bash
#SBATCH --job-name="anno_eval"
#SBATCH -t 150:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1
#SBATCH --mem=500GB
#SBATCH --exclusive

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#locating AUGUSTUS config path
export AUGUSTUS_CONFIG_PATH=/data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/config

#run BUSCO
busco \
--config /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/config.ini \
--in Pastreoides_transcripts_v1.fasta \
-o annotation_eval \
-l /data/putnamlab/kevin_wong1/Past_Genome/aug_training_2/busco_downloads/metazoa_odb10 \
-m transcriptome \
-f \
--long \
--augustus \
--augustus_parameters='--progress=true' \
--offline

echo "BUSCO Mission complete!" $(date)
```

# Functional Annotation

* Resources: Functional Annotation pipeline by D. Becker-Polinski [link](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-12-08-Molecular-Underpinnings-Functional-Annotation-Pipeline.md)

## Align query protein sequences against databases

1) BLAST the protein sequences against Swiss-Prot

`nano swissprot_blast.sh`

```bash
#!/bin/bash
#SBATCH --job-name="swissprot-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="swissprot_blastp_out_error"
#SBATCH --output="swissprot_blastp_out"
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1  
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database" $(date)

blastp -max_target_seqs 5 \
-num_threads 20 \
-db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 \
-query /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_proteins_v1.fasta \
-evalue 1e-5 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_vs_sprot_1e-5_max5.out

echo "STOP" $(date)

```

i) Get the best hit for each Gene Model (protein) Swiss-Prot

```
#Sort by 1. query name, 2. bitscore, 3. evalue, 4. protein identity, and extract the best line for each query (bitscore more important than evalue, evalue more important than nucleotide identity).

cat PastGeneModels_vs_sprot_1e-5_max5.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PastGeneModels_vs_sprot_1e-5_besthit.out

wc -l PastGeneModels_vs_sprot_1e-5_max5.out #239,922
wc -l PastGeneModels_vs_sprot_1e-5_besthit.out #30,487
```

ii) Select the gene model proteins without hits in Swiss-Prot


first use awk to print a list of all the Gene Model names from besthits.out
```
awk '{print $1}' PastGeneModels_vs_sprot_1e-5_besthit.out > list_of_Pastgenemodelproteins_sprot.txt

wc -l list_of_Pastgenemodelproteins_sprot.txt #30,487
```

Then exclude these Gene Model names from your original fasta/.faa/protein file.
* needed to load the module that has the script with the -exclude command in it

#first, loaded the newest module for kentUtils/416-foss-2020b

`module load kentUtils/416-foss-2020b`

#second use module show command to see paths to certain scripts and softwares in the module

`module show kentUtils/416-foss-2020b`

```
-------------------------------------------------------------------
/opt/modules/all/kentUtils/416-foss-2020b:

module-whatis     Description: LiftOver, Blat and other utilities
module-whatis     Homepage: https://hgdownload.soe.ucsc.edu/
module-whatis     URL: https://hgdownload.soe.ucsc.edu/
conflict     kentUtils
prepend-path     CMAKE_PREFIX_PATH /opt/software/kentUtils/416-foss-2020b
prepend-path     PATH /opt/software/kentUtils/416-foss-2020b/bin
setenv         EBROOTKENTUTILS /opt/software/kentUtils/416-foss-2020b
setenv         EBVERSIONKENTUTILS 416
setenv         EBDEVELKENTUTILS /opt/software/kentUtils/416-foss-2020b/easybuild/kentUtils-416-foss-2020b-easybuild-devel
-------------------------------------------------------------------
```

I selected to the prepend-path /opt/software/kentUtils/416-foss-2020b/bin to see if it took me to the 'faSomeRecords' script which it did

`/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords`

I then ran the -exclude command to exclude the blasted Gene Models from the .faa file

```
/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords -exclude /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_proteins_v1.fasta list_of_Pastgenemodelproteins_sprot.txt Past_proteins_names_v1.0.faa.prot4trembl
```

Checking the number of Gene Models:

`grep -c ">" Past_proteins_names_v1.0.faa.prot4trembl #34,149`

* use this file to blast against trembl

Downloading the .xml file for Blast2Go

`nano swissprot_blast_xml.sh`

```bash
#!/bin/bash
#SBATCH --job-name="swissprot-blastp-protein-xml"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="xml_blastp_out_error"
#SBATCH --output="xml_blastp_out"
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1  

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database with xml format out" $(date)
blastp -max_target_seqs 5 \
-num_threads 20 \
-db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 \
-query /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_proteins_v1.fasta \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_maxhit.xml

echo "STOP" $(date)
```

2) BLAST the remaining protein sequences against Trembl

`nano trembl_blastp.sh`

```bash
#!/bin/bash
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="trembl_blastp_out_error"
#SBATCH --output="trembl_blastp_out"
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against trembl database" $(date)
blastp -max_target_seqs 5 \
-num_threads 20 \
-db /data/putnamlab/shared/databases/trembl_db/trembl_20211022 \
-query Past_proteins_names_v1.0.faa.prot4trembl \
-evalue 1e-5 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_vs_trembl_1e-5_max5.out

echo "STOP" $(date)
```

`nano trembl_blastp_hit1.sh`

```bash
#!/bin/bash
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="trembl_hit1_blastp_out_error"
#SBATCH --output="trembl_hit1_blastp_out"
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1
#SBATCH -c 36

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against trembl database" $(date)
blastp -max_target_seqs 1 \
-num_threads 20 \
-db /data/putnamlab/shared/databases/trembl_db/trembl_20211022 \
-query Past_proteins_names_v1.0.faa.prot4trembl \
-evalue 1e-5 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_vs_trembl_1e-5_max1.out

echo "STOP" $(date)
```

Get the best hit for each Gene Model (protein) Trembl
```
cat PastGeneModels_vs_trembl_1e-5_max1.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PastGeneModels_vs_trembl_1e-5_besthit.out

wc -l PastGeneModels_vs_trembl_1e-5_besthit.out #24525
```

Select the gene model proteins without hits in Trembl
```bash
#first use awk to print a list of all the Gene Model names from besthits.out

awk '{print $1}' PastGeneModels_vs_trembl_1e-5_besthit.out > list_of_Pastgenemodelproteins_trembl.txt

#load the newest module for kentUtils/416-foss-2020b

module load kentUtils/416-foss-2020b

#then exclude these Gene Model names from your original fasta/.faa/protein file

/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords -exclude Past_proteins_names_v1.0.faa.prot4trembl list_of_Pastgenemodelproteins_trembl.txt Past_proteins_names_v1.0.faa.prot4nr


#check the number of Gene Models

grep -c ">" Past_proteins_names_v1.0.faa.prot4nr #9624

#using this file to blast against nr database

```

`nano trembl_blastp_xml.sh`

```bash
#!/bin/bash
#SBATCH --job-name="xml-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="xml_blastp_out_terror"
#SBATCH --output="xml_blastp_tout"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against trembl database for xml" $(date)
blastp -max_target_seqs 1 \
-num_threads 20 \
-db /data/putnamlab/shared/databases/trembl_db/trembl_20211022 \
-query Past_proteins_names_v1.0.faa.prot4trembl \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out Past_protein_blastp_trembl.xml

echo "STOP" $(date)
```

3) BLAST the remaining protein sequences against nr

`nano ncbi_blastp_out.sh`

```bash
#!/bin/bash
#SBATCH --job-name="ncbi-blastp-protein-out"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="ncbi_blastp_out_error"
#SBATCH --output="ncbi_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against ncbi database" $(date)
blastp -max_target_seqs 5 \
-num_threads $SLURM_CPUS_ON_NODE \
-db /data/shared/ncbi-nr/nr \
-query Past_proteins_names_v1.0.faa.prot4nr \
-evalue 1e-5 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_ncbi_max5hits.out

echo "STOP" $(date)
```


`nano ncbi_blastp.sh`

```bash
#!/bin/bash
#SBATCH --job-name="ncbi-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="ncbi_blastp_out_error"
#SBATCH --output="ncbi_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against ncbi database" $(date)
blastp -max_target_seqs 5 \
-num_threads $SLURM_CPUS_ON_NODE \
-db /data/shared/ncbi-nr/nr \
-query Past_proteins_names_v1.0.faa.prot4nr \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_ncbi.xml

echo "STOP" $(date)
```


## Interproscan

`sbatch Past_InterProScan.sh`

```bash
#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="interproscan_out_error"
#SBATCH --output="interproscan_out"
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1/InterProScan

echo "START $(date)"

# Load module
module load InterProScan/5.52-86.0-foss-2021a
module load Java/11.0.2
java -version

interproscan.sh --cpu $SLURM_CPUS_ON_NODE ...
interproscan.sh -version
interproscan.sh -f XML -i /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_proteins_v1.fasta -b Past.interpro.20220113  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i Past.interpro.20220113.xml -b Past.interpro.20220113

# -i is the input data
# -b is the output file base
# -f is formats
# -iprlookup enables mapping
# -goterms is GO Term
# -pa is pathway mapping
# -version displays version number

echo "DONE $(date)"
```

## Exporting xml Files

SwissProt
```
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1/PastGeneModels_maxhit.xml /Users/kevinwong/Desktop/URI_PHD/Projects/Past_genome/Functional_Annotation_files/
```

nr
```
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1/PastGeneModels_ncbi.xml /Users/kevinwong/Desktop/URI_PHD/Projects/Past_genome/Functional_Annotation_files/
```

trembl
```
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1/Past_protein_blastp_trembl.xml /Users/kevinwong/Desktop/URI_PHD/Projects/Past_genome/Functional_Annotation_files/
```

Interproscan
```
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/functional_anno_v1/InterProScan/Past.interpro.20220113.xml /Users/kevinwong/Desktop/URI_PHD/Projects/Past_genome/Functional_Annotation_files/
```

## Statistics with AGAT

``` bash
[kevin_wong1@n065 Plutea]$ module load AGAT/0.8.1-foss-2020b
[kevin_wong1@n065 Plutea]$ agat_sp_statistics.pl --gff plut2v1.1.genes.gff3
Reading file plut2v1.1.genes.gff3
Parsing: 100% [======================================================]D 0h01m12sParsing Finished
Compute statistics
--------------------------------------------------------------------------------

Compute mrna with isoforms if any

Number of genes                              31126
Number of mrnas                              31126
Number of mrnas with utr both sides          9305
Number of mrnas with at least one utr        15997
Number of cdss                               31126
Number of exons                              204364
Number of five_prime_utrs                    12044
Number of three_prime_utrs                   13258
Number of exon in cds                        197258
Number of exon in five_prime_utr             16114
Number of exon in three_prime_utr            16195
Number of intron in cds                      166132
Number of intron in exon                     173238
Number of intron in five_prime_utr           4070
Number of intron in three_prime_utr          2937
Number gene overlapping                      1069
Number of single exon gene                   5296
Number of single exon mrna                   5296
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          6.6
mean five_prime_utrs per mrna                0.4
mean three_prime_utrs per mrna               0.4
mean exons per cds                           6.3
mean exons per five_prime_utr                1.3
mean exons per three_prime_utr               1.2
mean introns in cdss per mrna                5.3
mean introns in exons per mrna               5.6
mean introns in five_prime_utrs per mrna     0.1
mean introns in three_prime_utrs per mrna    0.1
Total gene length                            255612213
Total mrna length                            255612213
Total cds length                             43135478
Total exon length                            56249572
Total five_prime_utr length                  2866380
Total three_prime_utr length                 10247714
Total intron length per cds                  190476728
Total intron length per exon                 199362641
Total intron length per five_prime_utr       6152481
Total intron length per three_prime_utr      2502850
mean gene length                             8212
mean mrna length                             8212
mean cds length                              1385
mean exon length                             275
mean five_prime_utr length                   237
mean three_prime_utr length                  772
mean cds piece length                        218
mean five_prime_utr piece length             177
mean three_prime_utr piece length            632
mean intron in cds length                    1146
mean intron in exon length                   1150
mean intron in five_prime_utr length         1511
mean intron in three_prime_utr length        852
Longest gene                                 143011
Longest mrna                                 143011
Longest cds                                  52230
Longest exon                                 12544
Longest five_prime_utr                       12253
Longest three_prime_utr                      11162
Longest cds piece                            12180
Longest five_prime_utr piece                 12253
Longest three_prime_utr piece                6747
Longest intron into cds part                 37855
Longest intron into exon part                37855
Longest intron into five_prime_utr part      28545
Longest intron into three_prime_utr part     13676
Shortest gene                                159
Shortest mrna                                159
Shortest cds                                 9
Shortest exon                                3
Shortest five_prime_utr                      1
Shortest three_prime_utr                     1
Shortest cds piece                           1
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                4
Shortest intron into exon part               4
Shortest intron into five_prime_utr part     20
Shortest intron into three_prime_utr part    20

Re-compute mrna without isoforms asked. We remove shortest isoforms if any

Number of genes                              31126
Number of mrnas                              31126
Number of mrnas with utr both sides          9305
Number of mrnas with at least one utr        15997
Number of cdss                               31126
Number of exons                              204364
Number of five_prime_utrs                    12044
Number of three_prime_utrs                   13258
Number of exon in cds                        197258
Number of exon in five_prime_utr             16114
Number of exon in three_prime_utr            16195
Number of intron in cds                      166132
Number of intron in exon                     173238
Number of intron in five_prime_utr           4070
Number of intron in three_prime_utr          2937
Number gene overlapping                      1069
Number of single exon gene                   5296
Number of single exon mrna                   5296
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          6.6
mean five_prime_utrs per mrna                0.4
mean three_prime_utrs per mrna               0.4
mean exons per cds                           6.3
mean exons per five_prime_utr                1.3
mean exons per three_prime_utr               1.2
mean introns in cdss per mrna                5.3
mean introns in exons per mrna               5.6
mean introns in five_prime_utrs per mrna     0.1
mean introns in three_prime_utrs per mrna    0.1
Total gene length                            255612213
Total mrna length                            255612213
Total cds length                             43135478
Total exon length                            56249572
Total five_prime_utr length                  2866380
Total three_prime_utr length                 10247714
Total intron length per cds                  190476728
Total intron length per exon                 199362641
Total intron length per five_prime_utr       6152481
Total intron length per three_prime_utr      2502850
mean gene length                             8212
mean mrna length                             8212
mean cds length                              1385
mean exon length                             275
mean five_prime_utr length                   237
mean three_prime_utr length                  772
mean cds piece length                        218
mean five_prime_utr piece length             177
mean three_prime_utr piece length            632
mean intron in cds length                    1146
mean intron in exon length                   1150
mean intron in five_prime_utr length         1511
mean intron in three_prime_utr length        852
Longest gene                                 143011
Longest mrna                                 143011
Longest cds                                  52230
Longest exon                                 12544
Longest five_prime_utr                       12253
Longest three_prime_utr                      11162
Longest cds piece                            12180
Longest five_prime_utr piece                 12253
Longest three_prime_utr piece                6747
Longest intron into cds part                 37855
Longest intron into exon part                37855
Longest intron into five_prime_utr part      28545
Longest intron into three_prime_utr part     13676
Shortest gene                                159
Shortest mrna                                159
Shortest cds                                 9
Shortest exon                                3
Shortest five_prime_utr                      1
Shortest three_prime_utr                     1
Shortest cds piece                           1
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                4
Shortest intron into exon part               4
Shortest intron into five_prime_utr part     20
Shortest intron into three_prime_utr part    20

--------------------------------------------------------------------------------

Bye Bye.
```

```bash
[kevin_wong1@n065 past_struc_annotations_v1]$ agat_sp_statistics.pl --gff Pastreoides_all_v1.gff
Reading file Pastreoides_all_v1.gff
Parsing: 100% [======================================================]D 0h10m33sParsing Finished
Compute statistics
--------------------------------------------------------------------------------

Compute match_part with isoforms if any

Number of matchs                             905337
Number of protein_matchs                     300777
Number of introns                            1062102
Number of match_parts                        2268216
Number gene overlapping                      0
Number of single exon match                  905337
Number of single exon protein_match          300777
mean introns per match                       1.2
mean introns per protein_match               3.5
mean match_parts per match                   2.5
mean match_parts per protein_match           7.5
Total match length                           939607611
Total protein_match length                   676247819
Total intron length                          1166327186
Total match_part length                      449528244
mean match length                            1037
mean protein_match length                    2248
mean intron length                           1098
mean match_part length                       198
Longest match                                104218
Longest protein_match                        140218
Longest intron                               46623
Longest match_part                           29944
Shortest match                               6
Shortest protein_match                       66
Shortest intron                              4
Shortest match_part                          2

Re-compute match_part without isoforms asked. We remove shortest isoforms if any

Number of matchs                             905337
Number of protein_matchs                     300777
Number of introns                            1062102
Number of match_parts                        2268216
Number gene overlapping                      0
Number of single exon match                  905337
Number of single exon protein_match          300777
mean introns per match                       1.2
mean introns per protein_match               3.5
mean match_parts per match                   2.5
mean match_parts per protein_match           7.5
Total match length                           939607611
Total protein_match length                   676247819
Total intron length                          1166327186
Total match_part length                      449528244
mean match length                            1037
mean protein_match length                    2248
mean intron length                           1098
mean match_part length                       198
Longest match                                104218
Longest protein_match                        140218
Longest intron                               46623
Longest match_part                           29944
Shortest match                               6
Shortest protein_match                       66
Shortest intron                              4
Shortest match_part                          2

--------------------------------------------------------------------------------

Compute mrna with isoforms if any

Number of genes                              64636
Number of mrnas                              64636
Number of mrnas with utr both sides          1623
Number of mrnas with at least one utr        8185
Number of cdss                               64636
Number of exons                              317304
Number of five_prime_utrs                    5227
Number of three_prime_utrs                   4581
Number of exon in cds                        315399
Number of exon in five_prime_utr             6439
Number of exon in three_prime_utr            5240
Number of intron in cds                      250763
Number of intron in exon                     252668
Number of intron in five_prime_utr           1212
Number of intron in three_prime_utr          659
Number gene overlapping                      1149
Number of single exon gene                   6221
Number of single exon mrna                   6221
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          4.9
mean five_prime_utrs per mrna                0.1
mean three_prime_utrs per mrna               0.1
mean exons per cds                           4.9
mean exons per five_prime_utr                1.2
mean exons per three_prime_utr               1.1
mean introns in cdss per mrna                3.9
mean introns in exons per mrna               3.9
mean introns in five_prime_utrs per mrna     0.0
mean introns in three_prime_utrs per mrna    0.0
Total gene length                            279311982
Total mrna length                            279311982
Total cds length                             58377165
Total exon length                            60411427
Total five_prime_utr length                  564947
Total three_prime_utr length                 1469315
Total intron length per cds                  216586326
Total intron length per exon                 218900555
Total intron length per five_prime_utr       1735255
Total intron length per three_prime_utr      537735
mean gene length                             4321
mean mrna length                             4321
mean cds length                              903
mean exon length                             190
mean five_prime_utr length                   108
mean three_prime_utr length                  320
mean cds piece length                        185
mean five_prime_utr piece length             87
mean three_prime_utr piece length            280
mean intron in cds length                    863
mean intron in exon length                   866
mean intron in five_prime_utr length         1431
mean intron in three_prime_utr length        815
Longest gene                                 105022
Longest mrna                                 105022
Longest cds                                  26259
Longest exon                                 10349
Longest five_prime_utr                       1687
Longest three_prime_utr                      2206
Longest cds piece                            10349
Longest five_prime_utr piece                 1501
Longest three_prime_utr piece                2206
Longest intron into cds part                 32611
Longest intron into exon part                32611
Longest intron into five_prime_utr part      11040
Longest intron into three_prime_utr part     7173
Shortest gene                                39
Shortest mrna                                39
Shortest cds                                 9
Shortest exon                                1
Shortest five_prime_utr                      1
Shortest three_prime_utr                     1
Shortest cds piece                           1
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                4
Shortest intron into exon part               4
Shortest intron into five_prime_utr part     5
Shortest intron into three_prime_utr part    5

Re-compute mrna without isoforms asked. We remove shortest isoforms if any

Number of genes                              64636
Number of mrnas                              64636
Number of mrnas with utr both sides          1623
Number of mrnas with at least one utr        8185
Number of cdss                               64636
Number of exons                              317304
Number of five_prime_utrs                    5227
Number of three_prime_utrs                   4581
Number of exon in cds                        315399
Number of exon in five_prime_utr             6439
Number of exon in three_prime_utr            5240
Number of intron in cds                      250763
Number of intron in exon                     252668
Number of intron in five_prime_utr           1212
Number of intron in three_prime_utr          659
Number gene overlapping                      1149
Number of single exon gene                   6221
Number of single exon mrna                   6221
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          4.9
mean five_prime_utrs per mrna                0.1
mean three_prime_utrs per mrna               0.1
mean exons per cds                           4.9
mean exons per five_prime_utr                1.2
mean exons per three_prime_utr               1.1
mean introns in cdss per mrna                3.9
mean introns in exons per mrna               3.9
mean introns in five_prime_utrs per mrna     0.0
mean introns in three_prime_utrs per mrna    0.0
Total gene length                            279311982
Total mrna length                            279311982
Total cds length                             58377165
Total exon length                            60411427
Total five_prime_utr length                  564947
Total three_prime_utr length                 1469315
Total intron length per cds                  216586326
Total intron length per exon                 218900555
Total intron length per five_prime_utr       1735255
Total intron length per three_prime_utr      537735
mean gene length                             4321
mean mrna length                             4321
mean cds length                              903
mean exon length                             190
mean five_prime_utr length                   108
mean three_prime_utr length                  320
mean cds piece length                        185
mean five_prime_utr piece length             87
mean three_prime_utr piece length            280
mean intron in cds length                    863
mean intron in exon length                   866
mean intron in five_prime_utr length         1431
mean intron in three_prime_utr length        815
Longest gene                                 105022
Longest mrna                                 105022
Longest cds                                  26259
Longest exon                                 10349
Longest five_prime_utr                       1687
Longest three_prime_utr                      2206
Longest cds piece                            10349
Longest five_prime_utr piece                 1501
Longest three_prime_utr piece                2206
Longest intron into cds part                 32611
Longest intron into exon part                32611
Longest intron into five_prime_utr part      11040
Longest intron into three_prime_utr part     7173
Shortest gene                                39
Shortest mrna                                39
Shortest cds                                 9
Shortest exon                                1
Shortest five_prime_utr                      1
Shortest three_prime_utr                     1
Shortest cds piece                           1
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                4
Shortest intron into exon part               4
Shortest intron into five_prime_utr part     5
Shortest intron into three_prime_utr part    5

--------------------------------------------------------------------------------

Bye Bye.
```

```
[kevin_wong1@n063 Prus]$ module load AGAT/0.8.1-foss-2020b
[kevin_wong1@n063 Prus]$ agat_sp_statistics.pl --gff Prus_gene_prediction.gff
Reading file Prus_gene_prediction.gff
Parsing: 100% [======================================================]D 0h01m13sParsing Finished
Compute statistics
--------------------------------------------------------------------------------

Compute mrna with isoforms if any

Number of genes                              39453
Number of mrnas                              39453
Number of cdss                               39453
Number of exons                              195533
Number of exon in cds                        195533
Number of intron in cds                      156080
Number of intron in exon                     156080
Number gene overlapping                      0
Number of single exon gene                   10066
Number of single exon mrna                   10066
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          5.0
mean exons per cds                           5.0
mean introns in cdss per mrna                4.0
mean introns in exons per mrna               4.0
Total gene length                            237912951
Total mrna length                            237912951
Total cds length                             46736451
Total exon length                            46736451
Total intron length per cds                  191176500
Total intron length per exon                 191176500
mean gene length                             6030
mean mrna length                             6030
mean cds length                              1184
mean exon length                             239
mean cds piece length                        239
mean intron in cds length                    1224
mean intron in exon length                   1224
Longest gene                                 149486
Longest mrna                                 149486
Longest cds                                  36834
Longest exon                                 9978
Longest cds piece                            9978
Longest intron into cds part                 10000
Longest intron into exon part                10000
Shortest gene                                300
Shortest mrna                                300
Shortest cds                                 300
Shortest exon                                3
Shortest cds piece                           3
Shortest intron into cds part                20
Shortest intron into exon part               20

Re-compute mrna without isoforms asked. We remove shortest isoforms if any

Number of genes                              39453
Number of mrnas                              39453
Number of cdss                               39453
Number of exons                              195533
Number of exon in cds                        195533
Number of intron in cds                      156080
Number of intron in exon                     156080
Number gene overlapping                      0
Number of single exon gene                   10066
Number of single exon mrna                   10066
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          5.0
mean exons per cds                           5.0
mean introns in cdss per mrna                4.0
mean introns in exons per mrna               4.0
Total gene length                            237912951
Total mrna length                            237912951
Total cds length                             46736451
Total exon length                            46736451
Total intron length per cds                  191176500
Total intron length per exon                 191176500
mean gene length                             6030
mean mrna length                             6030
mean cds length                              1184
mean exon length                             239
mean cds piece length                        239
mean intron in cds length                    1224
mean intron in exon length                   1224
Longest gene                                 149486
Longest mrna                                 149486
Longest cds                                  36834
Longest exon                                 9978
Longest cds piece                            9978
Longest intron into cds part                 10000
Longest intron into exon part                10000
Shortest gene                                300
Shortest mrna                                300
Shortest cds                                 300
Shortest exon                                3
Shortest cds piece                           3
Shortest intron into cds part                20
Shortest intron into exon part               20

--------------------------------------------------------------------------------
```
```
[kevin_wong1@n064 Paus]$ module load AGAT/0.8.1-foss-2020b
[kevin_wong1@n064 Paus]$ agat_sp_statistics.pl --gff paus_mRNA.gff
Reading file paus_mRNA.gff
Parsing: 100% [======================================================]D 0h03m03sParsing Finished
Compute statistics
--------------------------------------------------------------------------------

Compute transcript with isoforms if any

Number of genes                              30603
Number of transcripts                        35910
Number of mrnas with utr both sides          35517
Number of mrnas with at least one utr        35907
Number of cdss                               35910
Number of exons                              314893
Number of five_prime_utrs                    35681
Number of introns                            252017
Number of start_codons                       34552
Number of stop_codons                        34675
Number of three_prime_utrs                   35743
Number of transcription_end_sites            34542
Number of exon in cds                        285689
Number of exon in five_prime_utr             56405
Number of exon in three_prime_utr            43838
Number of intron in cds                      249779
Number of intron in exon                     278983
Number of intron in five_prime_utr           20724
Number of intron in intron                   222367
Number of intron in three_prime_utr          8095
Number gene overlapping                      593
Number of single exon gene                   3970
Number of single exon transcript             4021
mean transcripts per gene                    1.2
mean cdss per transcript                     1.0
mean exons per transcript                    8.8
mean five_prime_utrs per transcript          1.0
mean introns per transcript                  7.0
mean start_codons per transcript             1.0
mean stop_codons per transcript              1.0
mean three_prime_utrs per transcript         1.0
mean transcription_end_sites per transcript  1.0
mean exons per cds                           8.0
mean exons per five_prime_utr                1.6
mean exons per three_prime_utr               1.2
mean introns in cdss per transcript          7.0
mean introns in exons per transcript         7.8
mean introns in five_prime_utrs per transcript0.6
mean introns in introns per transcript       6.2
mean introns in three_prime_utrs per transcript0.2
Total gene length                            348974335
Total transcript length                      463448447
Total cds length                             60898458
Total exon length                            99450105
Total five_prime_utr length                  10676138
Total intron length                          301789183
Total start_codon length                     103656
Total stop_codon length                      104025
Total three_prime_utr length                 27875509
Total transcription_end_site length          34542
Total intron length per cds                  297498344
Total intron length per exon                 363998342
Total intron length per five_prime_utr       49561106
Total intron length per intron               37321291
Total intron length per three_prime_utr      15710850
mean gene length                             11403
mean transcript length                       12905
mean cds length                              1695
mean exon length                             315
mean five_prime_utr length                   299
mean intron length                           1197
mean start_codon length                      3
mean stop_codon length                       3
mean three_prime_utr length                  779
mean transcription_end_site length           1
mean cds piece length                        213
mean five_prime_utr piece length             189
mean three_prime_utr piece length            635
mean intron in cds length                    1191
mean intron in exon length                   1304
mean intron in five_prime_utr length         2391
mean intron in intron length                 167
mean intron in three_prime_utr length        1940
Longest gene                                 232931
Longest transcript                           232931
Longest cds                                  63957
Longest exon                                 18417
Longest five_prime_utr                       16693
Longest intron                               87093
Longest start_codon                          3
Longest stop_codon                           3
Longest three_prime_utr                      18129
Longest transcription_end_site               1
Longest cds piece                            14759
Longest five_prime_utr piece                 16693
Longest three_prime_utr piece                18129
Longest intron into cds part                 87093
Longest intron into exon part                87093
Longest intron into five_prime_utr part      57732
Longest intron into intron part              14759
Longest intron into three_prime_utr part     42549
Shortest gene                                174
Shortest transcript                          174
Shortest cds                                 5
Shortest exon                                3
Shortest five_prime_utr                      1
Shortest intron                              6
Shortest start_codon                         3
Shortest stop_codon                          3
Shortest three_prime_utr                     1
Shortest transcription_end_site              1
Shortest cds piece                           3
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                50
Shortest intron into exon part               50
Shortest intron into five_prime_utr part     50
Shortest intron into intron part             4
Shortest intron into three_prime_utr part    50

Re-compute transcript without isoforms asked. We remove shortest isoforms if any

Number of genes                              30603
Number of transcripts                        30603
Number of mrnas with utr both sides          30242
Number of mrnas with at least one utr        30601
Number of cdss                               30603
Number of exons                              231831
Number of five_prime_utrs                    30396
Number of introns                            179972
Number of start_codons                       29381
Number of stop_codons                        29514
Number of three_prime_utrs                   30447
Number of transcription_end_sites            29393
Number of exon in cds                        208588
Number of exon in five_prime_utr             46979
Number of exon in three_prime_utr            36768
Number of intron in cds                      177985
Number of intron in exon                     201228
Number of intron in five_prime_utr           16583
Number of intron in intron                   155513
Number of intron in three_prime_utr          6321
Number gene overlapping                      489
Number of single exon gene                   3976
Number of single exon transcript             3976
mean transcripts per gene                    1.0
mean cdss per transcript                     1.0
mean exons per transcript                    7.6
mean five_prime_utrs per transcript          1.0
mean introns per transcript                  5.9
mean start_codons per transcript             1.0
mean stop_codons per transcript              1.0
mean three_prime_utrs per transcript         1.0
mean transcription_end_sites per transcript  1.0
mean exons per cds                           6.8
mean exons per five_prime_utr                1.5
mean exons per three_prime_utr               1.2
mean introns in cdss per transcript          5.8
mean introns in exons per transcript         6.6
mean introns in five_prime_utrs per transcript0.5
mean introns in introns per transcript       5.1
mean introns in three_prime_utrs per transcript0.2
Total gene length                            348974335
Total transcript length                      346221435
Total cds length                             48222964
Total exon length                            80303427
Total five_prime_utr length                  9038176
Total intron length                          215643389
Total start_codon length                     88143
Total stop_codon length                      88542
Total three_prime_utr length                 23042287
Total transcription_end_site length          29393
Total intron length per cds                  211904194
Total intron length per exon                 265918008
Total intron length per five_prime_utr       39962084
Total intron length per intron               27047457
Total intron length per three_prime_utr      13007009
mean gene length                             11403
mean transcript length                       11313
mean cds length                              1575
mean exon length                             346
mean five_prime_utr length                   297
mean intron length                           1198
mean start_codon length                      3
mean stop_codon length                       3
mean three_prime_utr length                  756
mean transcription_end_site length           1
mean cds piece length                        231
mean five_prime_utr piece length             192
mean three_prime_utr piece length            626
mean intron in cds length                    1190
mean intron in exon length                   1321
mean intron in five_prime_utr length         2409
mean intron in intron length                 173
mean intron in three_prime_utr length        2057
Longest gene                                 232931
Longest transcript                           232931
Longest cds                                  63957
Longest exon                                 18417
Longest five_prime_utr                       16693
Longest intron                               87093
Longest start_codon                          3
Longest stop_codon                           3
Longest three_prime_utr                      18129
Longest transcription_end_site               1
Longest cds piece                            14759
Longest five_prime_utr piece                 16693
Longest three_prime_utr piece                18129
Longest intron into cds part                 87093
Longest intron into exon part                87093
Longest intron into five_prime_utr part      57732
Longest intron into intron part              14759
Longest intron into three_prime_utr part     42549
Shortest gene                                174
Shortest transcript                          174
Shortest cds                                 5
Shortest exon                                3
Shortest five_prime_utr                      1
Shortest intron                              6
Shortest start_codon                         3
Shortest stop_codon                          3
Shortest three_prime_utr                     1
Shortest transcription_end_site              1
Shortest cds piece                           3
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                50
Shortest intron into exon part               50
Shortest intron into five_prime_utr part     50
Shortest intron into intron part             4
Shortest intron into three_prime_utr part    50

--------------------------------------------------------------------------------

Bye Bye.
```
