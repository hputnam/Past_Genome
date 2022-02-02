# Testing Multimapping to the P.astreoides genome

```
mkdir multimapping
cd multimapping
mkdir Past_genome_ref
```

## Sequences:
* KW RNA seq reads from the same genome coral
* Link [here](https://github.com/hputnam/Past_Genome/blob/master/denovo_transcriptome_assembly_BUSCO.md) for the scripts for QC.

## HISAT2 Script

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
ln -s /data/putnamlab/kevin_wong1/20201221_P.astreoides_Ref_Transcriptome/cleaned_reads/clean*.fastq.gz ./

# index the reference genome for Porites astreoides output index to working directory
hisat2-build -f /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta ./Past_genome_ref # called the reference genome (scaffolds)
echo "Reference genome indexed. Starting alignment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --dta --new-summary -x Past_genome_ref -U ${i} -S ${sample_name}.sam 2> "${sample_name}"_hisat2.err
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

echo "Mission Complete!" $(date)
```

## HISAT2 Alignment Statistics

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
