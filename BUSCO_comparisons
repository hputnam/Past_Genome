# BUSCO Comparisons

### Genome Comparisons

Porites lutea
* Reference
* Link

```
cp plut_final_2.1.fasta.gz ../../kevin_wong1/Past_Genome/refs/
gunzip plut_final_2.1.fasta.gz
```

Porites rus
* Downloaded 20220131
* Reference
* Link

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/290/455/GCA_900290455.1_Prus/GCA_900290455.1_Prus_genomic.fna.gz
gunzip GCA_900290455.1_Prus_genomic.fna.gz
```

```
mkdir BUSCO Compare
cd BUSCO_Compare
mkdir genomes
mkdir transcriptomes
```

### Porites rus genome

`nano BUSCO_prus.sh`

```
#!/bin/bash
#SBATCH --job-name="Prus_BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/BUSCO_Compare/genomes
#SBATCH --mem=50GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#run BUSCO
busco --config config.ini \
-m genome \
-i /data/putnamlab/kevin_wong1/Past_Genome/refs/GCA_900290455.1_Prus_genomic.fna \
-o Prus_BUSCO \
-l /data/putnamlab/kevin_wong1/busco_downloads/metazoa_odb10 \
--offline

echo "BUSCO Mission complete!" $(date)
```

```
# BUSCO version is: 5.2.2
# The lineage dataset is: metazoa_odb10 (Creation date: 2020-09-10, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/Past_Genome/refs/GCA_900290455.1_Prus_genomic.fna
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

        ***** Results: *****

        C:67.8%[S:65.5%,D:2.3%],F:18.7%,M:13.5%,n:954      
        647     Complete BUSCOs (C)                        
        625     Complete and single-copy BUSCOs (S)        
        22      Complete and duplicated BUSCOs (D)         
        178     Fragmented BUSCOs (F)                      
        129     Missing BUSCOs (M)                         
        954     Total BUSCO groups searched   
```

### Porites lutea genome

`nano BUSCO_plut.sh`

```
#!/bin/bash
#SBATCH --job-name="Plut_BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/BUSCO_Compare/genomes
#SBATCH --mem=50GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#run BUSCO
busco --config config.ini \
-m genome \
-i /data/putnamlab/kevin_wong1/Past_Genome/refs/plut_final_2.1.fasta \
-o Plut_BUSCO \
-l /data/putnamlab/kevin_wong1/busco_downloads/metazoa_odb10 \
--offline

echo "BUSCO Mission complete!" $(date)
```

```
# BUSCO version is: 5.2.2
# The lineage dataset is: metazoa_odb10 (Creation date: 2020-09-10, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/Past_Genome/refs/plut_final_2.1.fasta
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

        ***** Results: *****

        C:93.7%[S:91.9%,D:1.8%],F:2.8%,M:3.5%,n:954        
        894     Complete BUSCOs (C)                        
        877     Complete and single-copy BUSCOs (S)        
        17      Complete and duplicated BUSCOs (D)         
        27      Fragmented BUSCOs (F)                      
        33      Missing BUSCOs (M)                         
        954     Total BUSCO groups searched  
```

### Porites astreoides genome

`nano BUSCO_past.sh`

```
#!/bin/bash
#SBATCH --job-name="Past_BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/BUSCO_Compare/genomes
#SBATCH --mem=50GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#run BUSCO
busco --config config.ini \
-m genome \
-i /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
-o Past_BUSCO \
-l /data/putnamlab/kevin_wong1/busco_downloads/metazoa_odb10 \
--offline

echo "BUSCO Mission complete!" $(date)
```

```
# BUSCO version is: 5.2.2
# The lineage dataset is: metazoa_odb10 (Creation date: 2020-09-10, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

        ***** Results: *****

        C:90.9%[S:77.5%,D:13.4%],F:4.1%,M:5.0%,n:954       
        867     Complete BUSCOs (C)                        
        739     Complete and single-copy BUSCOs (S)        
        128     Complete and duplicated BUSCOs (D)         
        39      Fragmented BUSCOs (F)                      
        48      Missing BUSCOs (M)                         
        954     Total BUSCO groups searched    
```
