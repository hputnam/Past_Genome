# BUSCO Comparisons

## Genome Comparisons

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


## Transcriptome Comparisons

### Porites astreoides (Wong and Putnam 2022)

`nano Wong_BUSCO.sh`

```
#!/bin/bash
#SBATCH --job-name="Wong_BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/BUSCO_Compare/transcriptomes
#SBATCH --mem=50GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#run BUSCO
busco --config config.ini \
-m transcriptome \
-i /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_transcripts_v1.fasta \
-o Wong_BUSCO \
-l /data/putnamlab/kevin_wong1/busco_downloads/metazoa_odb10 \
--offline

echo "BUSCO Mission complete!" $(date)
```
```
# BUSCO version is: 5.2.2
# The lineage dataset is: metazoa_odb10 (Creation date: 2020-09-10, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_transcripts_v1.fasta
# BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:79.6%[S:69.6%,D:10.0%],F:12.1%,M:8.3%,n:954      
        759     Complete BUSCOs (C)                        
        664     Complete and single-copy BUSCOs (S)        
        95      Complete and duplicated BUSCOs (D)         
        115     Fragmented BUSCOs (F)                      
        80      Missing BUSCOs (M)                         
        954     Total BUSCO groups searched                
```


### Porites astreoides (Kenkel et al 2013)

`nano Kenkel_BUSCO.sh`

```
#!/bin/bash
#SBATCH --job-name="Kenkel_BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/BUSCO_Compare/transcriptomes
#SBATCH --mem=50GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#run BUSCO
busco --config config.ini \
-m transcriptome \
-i /data/putnamlab/kevin_wong1/Past_Genome/refs/Kenkel2013_past_transcriptome.fasta \
-o Kenkel_BUSCO \
-l /data/putnamlab/kevin_wong1/busco_downloads/metazoa_odb10 \
--offline

echo "BUSCO Mission complete!" $(date)
```

```
# BUSCO version is: 5.2.2
# The lineage dataset is: metazoa_odb10 (Creation date: 2020-09-10, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/Past_Genome/refs/Kenkel2013_past_transcriptome.fasta
# BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:22.9%[S:22.2%,D:0.7%],F:31.1%,M:46.0%,n:954      
        219     Complete BUSCOs (C)                        
        212     Complete and single-copy BUSCOs (S)        
        7       Complete and duplicated BUSCOs (D)         
        297     Fragmented BUSCOs (F)                      
        438     Missing BUSCOs (M)                         
        954     Total BUSCO groups searched    
```

### Porites astreoides (Mansour et al 2016)

`nano Mansour_BUSCO.sh`

```
#!/bin/bash
#SBATCH --job-name="Mansour_BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/BUSCO_Compare/transcriptomes
#SBATCH --mem=50GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#run BUSCO
busco --config config.ini \
-m transcriptome \
-i /data/putnamlab/kevin_wong1/REFS/Past_Mansour/p_ast2016.fasta \
-o Mansour_BUSCO \
-l /data/putnamlab/kevin_wong1/busco_downloads/metazoa_odb10 \
--offline

echo "BUSCO Mission complete!" $(date)
```

```
BUSCO version is: 5.2.2
# The lineage dataset is: metazoa_odb10 (Creation date: 2020-09-10, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/REFS/Past_Mansour/p_ast2016.fasta
# BUSCO was run in mode: transcriptome

       ***** Results: *****

       C:30.3%[S:18.1%,D:12.2%],F:4.9%,M:64.8%,n:954      
       289     Complete BUSCOs (C)                        
       173     Complete and single-copy BUSCOs (S)        
       116     Complete and duplicated BUSCOs (D)         
       47      Fragmented BUSCOs (F)                      
       618     Missing BUSCOs (M)                         
       954     Total BUSCO groups searched   
```


### Porites astreoides (Walker et al 2019)

https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/CX9HWD

```
wget -O CoralTranscripts_Transcriptome_FINAL.fasta "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/CX9HWD/YDFYYK"

mv CoralTranscripts_Transcriptome_FINAL.fasta Walker_past_transcriptome.fasta
```

For some reason, this fasta file has duplicate of sequence >c84342_g1_i1 in it, so it will not run BUSCO. To remove the duplicate sequence I am doing the following:

`awk 'BEGIN {i = 1;} { if ($1 ~ /^>/) { tmp = h[i]; h[i] = $1; } else if (!a[$1]) { s[i] = $1; a[$1] = "1"; i++; } else { h[i] = tmp; } } END { for (j = 1; j < i; j++) { print h[j]; print s[j]; } }' < Walker_past_transcriptome.fasta > Walker_past_transcriptome_clean.fasta`


`nano Walker_BUSCO.sh`

```
#!/bin/bash
#SBATCH --job-name="Walker_BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/BUSCO_Compare/transcriptomes
#SBATCH --mem=50GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#run BUSCO
busco --config config.ini \
-m transcriptome \
-i /data/putnamlab/kevin_wong1/Past_Genome/refs/Walker_past_transcriptome_clean.fasta \
-o Walker_BUSCO \
-l /data/putnamlab/kevin_wong1/busco_downloads/metazoa_odb10 \
--offline

echo "BUSCO Mission complete!" $(date)
```
```
# BUSCO version is: 5.2.2
# The lineage dataset is: metazoa_odb10 (Creation date: 2020-09-10, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/Past_Genome/refs/Walker_past_transcriptome_clean.fasta
# BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:63.1%[S:26.5%,D:36.6%],F:21.6%,M:15.3%,n:954     
        602     Complete BUSCOs (C)                        
        253     Complete and single-copy BUSCOs (S)        
        349     Complete and duplicated BUSCOs (D)         
        206     Fragmented BUSCOs (F)                      
        146     Missing BUSCOs (M)                         
        954     Total BUSCO groups searched   
```
