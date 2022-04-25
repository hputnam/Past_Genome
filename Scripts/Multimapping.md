# Testing Multimapping to the *P. astreoides genome*

Steps:
1. Map RNAseq data to PAST ab initio genome `.bam` (STAR)
2. Obtain mapping genome coordinates `genome.fasta` (Stringtie)
3. Convert to transcript coordinates `transcript.fasta` (gffread)
4. Map `transcript.fasta` to PAST ab initio genome

## 1. Map RNAseq data to PAST ab initio genome `.bam`

### Mapping with HISAT2

```
mkdir multimapping
cd multimapping
mkdir Past_genome_ref
```

#### Sequences:
* KW RNA seq reads from the same genome coral
* Link [here](https://github.com/hputnam/Past_Genome/blob/master/denovo_transcriptome_assembly_BUSCO.md) for the scripts for QC.

#### HISAT2 Script

`nano hisat2.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/multimapping
#SBATCH --cpus-per-task=3

#load packages
module load HISAT2/2.1.0-foss-2018b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# symbolically link 'clean' reads to hisat2 dir
ln -s /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/cleaned_reads/clean*.fasta ./

# index the reference genome for Porites astreoides output index to working directory
hisat2-build -f /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta ./Past_genome_ref # called the reference genome (scaffolds)
echo "Reference genome indexed. Starting alignment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fasta)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --dta --new-summary -x Past_genome_ref -U ${i} -S ${sample_name}.sam 2> "${sample_name}"_hisat2.err
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

echo "Mission Complete!" $(date)
```

#### HISAT2 Alignment Statistics

| Sequence Name | Total Reads |   Aligned 0 time  |   Aligned 1 time   |  Aligned >1 times | Overall alignment rate |
|:-------------:|:-----------:|:-----------------:|:------------------:|:-----------------:|:----------------------:|
|   Sample2_R1  |   8726675   |  2535524 (29.05%) |   4612342 (52.85%) |  1578809 (18.09%) |         70.95%         |
|   Sample2_R2  |   8726675   |  2706082 (31.01%) |   4470937 (51.23%) |  1549656 (17.76%) |         68.99%         |
|   Sample3_R1  |   7255348   |  1267524 (17.47%) |   4501001 (62.04%) |  1486823 (20.49%) |         82.53%         |
|   Sample3_R2  |   7255348   |  1403913 (19.35%) |   4385847 (60.45%) |  1465588 (20.20%) |         80.65%         |
|   Sample4_R1  |   7237699   |  1593178 (22.01%) |   4255093 (58.79%) |  1389428 (19.20%) |         77.99%         |
|   Sample4_R2  |   7237699   |  1710164 (23.63%) |   4153467 (57.39%) |  1374068 (18.98%) |         76.37%         |
|   Sample5_R1  |   7673755   |  1482145 (19.31%) |   4677711 (60.96%) |  1513899 (19.73%) |         80.69%         |
|   Sample5_R2  |   7673755   |  1597815 (20.82%) |   4577756 (59.65%) |  1498184 (19.52%) |         79.18%         |
| Sample_ALL_R1 |   30893477  |  6870864 (22.24%) |  18054756 (58.44%) |  5967857 (19.32%) |         77.76%         |
| Sample_ALL_R2 |   30893477  |  7410565 (23.99%) |  17596361 (56.96%) |  5886551 (19.05%) |         76.01%         |


```
scp -r kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/Past_Genome/multimapping/bam_files/ .

scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gff .
scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta .
scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_transcripts_v1.fasta .

```

## Mapping statistics with STAR

#### Sequences:

1.  5 KW RNA seq reads from the same genome coral
  * Link [here](https://github.com/hputnam/Past_Genome/blob/master/denovo_transcriptome_assembly_BUSCO.md) for the scripts for QC.
2. Kenkel 2013 transcriptome
3. Mansour 2016 transcriptome
4. Walker 2019 transcriptome
5. Wong 2022 transcriptome from MAKER


#### Generate Genome Index

```
mkdir STAR
cd STAR
mkdir GenomeIndex_Past
```

`nano genome_index.sh`

```
#!/bin/bash
#SBATCH --job-name="STAR_index"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="index_error"
#SBATCH --output="index_out"
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR

module load STAR/2.7.2b-GCC-8.3.0

# Indexing genome
STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/GenomeIndex_Past \
--genomeFastaFiles /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--sjdbGTFfile /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gff
```

#### Align to reference genome

Converting fastq files to fasta

```
cd /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/cleaned_reads

gunzip -c clean.Sample2_R1.fastq.gz  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > clean.Sample2_R1.fasta
gunzip -c clean.Sample2_R2.fastq.gz  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > clean.Sample2_R2.fasta
gunzip -c clean.Sample3_R1.fastq.gz  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > clean.Sample3_R1.fasta
gunzip -c clean.Sample3_R2.fastq.gz  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > clean.Sample3_R2.fasta
gunzip -c clean.Sample4_R1.fastq.gz  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > clean.Sample4_R1.fasta
gunzip -c clean.Sample4_R2.fastq.gz  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > clean.Sample4_R2.fasta
gunzip -c clean.Sample5_R1.fastq.gz  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > clean.Sample5_R1.fasta
gunzip -c clean.Sample5_R2.fastq.gz  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > clean.Sample5_R2.fasta
```

##### Aligning RNAseq reads

`mkdir align_rnaseq`

`nano alignrnaseq_star.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_STAR_out_error"
#SBATCH --output="Align_STAR_out"
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/align_rnaseq

# symbolically link to existing directory
ln -s /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/cleaned_reads/*fasta ./ #RNAseq reads

# Align reads to genome
module load STAR/2.7.2b-GCC-8.3.0

F=/data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/align_rnaseq

array1=($(ls $F/*_R1.fasta))
for i in ${array1[@]}
do
STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir ${i}_TMP \
--readFilesIn ${i} $(echo ${i}|sed s/_R1/_R2/) \
--genomeDir /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/GenomeIndex_Past \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $(echo ${i}|sed s/_R1.fasta//).
done
```


##### Aligning transcriptomes
`mkdir align`

`nano align_star.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_STAR_out_error"
#SBATCH --output="Align_STAR_out"
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/align

# symbolically link to existing directory
ln -s /data/putnamlab/kevin_wong1/Past_Genome/refs/Kenkel2013_past_transcriptome.fasta ./ #Kenkel 2013 transcriptome  
ln -s /data/putnamlab/kevin_wong1/REFS/Past_Mansour/p_ast2016.fasta ./ #Mansour 2016 transcriptome
ln -s /data/putnamlab/kevin_wong1/Past_Genome/refs/Walker_past_transcriptome_clean.fasta ./ #Walker 2019 transcriptome
ln -s /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_transcripts_v1.fasta  ./ #Wong 2022 transcriptome from MAKER

# Align reads to genome
module load STAR/2.7.2b-GCC-8.3.0

F=/data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/align

array1=($(ls $F/*fast*))
for i in ${array1[@]}
do
STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir ${i}_TMP \
--readFilesIn ${i} \
--genomeDir /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/GenomeIndex_Past \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${i}.
done
```

![ ](https://github.com/hputnam/Past_Genome/blob/master/images/STAR_mapping_stats.png)


## 2. Obtain mapping genome coordinates using Stringtie

### Stringtie

`mkdir stringtie`

`nano stringtie.sh`

```bash
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/multimapping/stringtie
#SBATCH --cpus-per-task=3

module load StringTie/2.1.4-GCC-9.3.0

# symbolically link to existing directory
ln -s /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/align_rnaseq/*sortedByCoord.out.bam ./

# Assemble
array1=($(ls *.bam))
for i in ${array1[@]}; do
stringtie -p 48 --rf -e -G /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gff  -o ./${i}.gtf ./${i}
done
```


## 3. Convert to transcript coordinates `transcript.fasta` (gffread)

### gffread

`nano gffread.sh`

```bash
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/multimapping/stringtie
#SBATCH --cpus-per-task=3

module load gffread/0.12.7-GCCcore-11.2.0

array1=($(ls *.Aligned.sortedByCoord.out.bam.gtf))
for i in ${array1[@]}; do
gffread -w ./${i}.transcripts.fa -g /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta ./${i}
done
```

## 4. Map `transcript.fasta` to PAST ab initio genome

### STAR

`cd ../STAR`

`mkdir align_transcripts`

`nano align_transcripts.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_STAR_out_error"
#SBATCH --output="Align_STAR_out"
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/align_transcripts

# symbolically link to existing directory
ln -s /data/putnamlab/kevin_wong1/Past_Genome/multimapping/stringtie/*.fa ./

# Align reads to genome
module load STAR/2.7.2b-GCC-8.3.0

F=/data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/align_transcripts

array1=($(ls $F/*.fa))
for i in ${array1[@]}
do
STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir ${i}_TMP \
--readFilesIn ${i} \
--genomeDir /data/putnamlab/kevin_wong1/Past_Genome/multimapping/STAR/GenomeIndex_Past \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${i}.
done
```

`less clean.Sample2.Aligned.sortedByCoord.out.bam.gtf.transcripts.fa.Log.final.out`

```bash
                                 Started job on |       Mar 09 11:30:44
                             Started mapping on |       Mar 09 11:31:55
                                    Finished on |       Mar 09 11:31:58
       Mapping speed, Million of reads per hour |       0.00

                          Number of input reads |       4
                      Average input read length |       532
                                    UNIQUE READS:
                   Uniquely mapped reads number |       1
                        Uniquely mapped reads % |       25.00%
                          Average mapped length |       649.00
                       Number of splices: Total |       3
            Number of splices: Annotated (sjdb) |       3
                       Number of splices: GT/AG |       3
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       0
                      Mismatch rate per base, % |       0.00%
                         Deletion rate per base |       0.00%
                        Deletion average length |       0.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       0.00
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       2
             % of reads mapped to multiple loci |       50.00%
        Number of reads mapped to too many loci |       1
             % of reads mapped to too many loci |       25.00%
                                  UNMAPPED READS:
  Number of reads unmapped: too many mismatches |       0
       % of reads unmapped: too many mismatches |       0.00%
            Number of reads unmapped: too short |       0
                 % of reads unmapped: too short |       0.00%
                Number of reads unmapped: other |       0
                     % of reads unmapped: other |       0.00%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
```
