# Functional Annotation - Combining blast2go outputs
# Author: Kevin Wong
# Last Modified: 20220310

# Load Packages
library("dplyr")
library("tidyverse")
library("ggplot2")

# Load Data
ncbi_blast <- read.csv("output/Functional_Annotation/ncbi_blast2go_table.csv", check.names = FALSE)
swissprot <- read.csv("output/Functional_Annotation/SwissProt_blast2go_table.csv", check.names = FALSE)
trembl <- read.csv("output/Functional_Annotation/trembl_blast2go_table.csv", check.names = FALSE)

# Replacing "spaces" in column headers with "_"
names(ncbi_blast) <- gsub(" ", "_", names(ncbi_blast))
names(swissprot) <- gsub(" ", "_", names(swissprot))
names(trembl) <- gsub(" ", "_", names(trembl))

# Replacing "#" in column headers with "num"
names(ncbi_blast) <- gsub("#", "num_", names(ncbi_blast))
names(swissprot) <- gsub("#", "num_", names(swissprot))
names(trembl) <- gsub("#", "num_", names(trembl))

# Selecting columns that will remain in the final data frame
const_df <- ncbi_blast %>%
  select(SeqName, Length, InterPro_IDs,	InterPro_GO_IDs,	InterPro_GO_Names)

# Selecting columns of interest and modifying header names in each data frame
ncbi_blast2 <- ncbi_blast %>%
  select(SeqName, Description, num_Hits,	`e-Value`,	sim_mean, num_GO, GO_IDs, GO_Names, Enzyme_Codes, Enzyme_Names)

colnames(ncbi_blast2)[2:10] <- paste("BLAST", colnames(ncbi_blast2)[2:10], sep = "_")

swissprot2 <- swissprot %>%
  select(SeqName, Description, num_Hits,	`e-Value`,	sim_mean, num_GO, GO_IDs, GO_Names, Enzyme_Codes, Enzyme_Names)

colnames(swissprot2)[2:10] <- paste("SwissProt", colnames(swissprot2)[2:10], sep = "_")

trembl2 <- trembl %>%
  select(SeqName, Description, num_Hits,	`e-Value`,	sim_mean, num_GO, GO_IDs, GO_Names, Enzyme_Codes, Enzyme_Names)

colnames(trembl2)[2:10] <- paste("TrEMBL", colnames(trembl2)[2:10], sep = "_")

# Merging all data frames together
df_list <- list(const_df, ncbi_blast2, swissprot2, trembl2)
Past_annotation_full <- df_list %>% reduce(full_join, by='SeqName')
rownames(Past_annotation_full)<-NULL

write_csv(Past_annotation_full, "output/Functional_Annotation/Past_annotation_20220310.csv")

## Finding # of annotations from each database

Past_annotation_full <- read.csv("output/Functional_Annotation/Past_annotation_20220310.csv", check.names = FALSE)

ncbi_blast3 <- Past_annotation_full %>%
  select(SeqName, BLAST_num_Hits) %>%
  filter(BLAST_num_Hits != "NA") # 25513

swissprot3 <- Past_annotation_full %>%
  select(SeqName, SwissProt_num_Hits) %>%
  filter(SwissProt_num_Hits != "NA") # 30444

trembl3 <- Past_annotation_full %>%
  select(SeqName, TrEMBL_num_Hits) %>%
  filter(TrEMBL_num_Hits != "NA") # 24359 

all <- 25513 + 30444 + 24359 

test1 <- semi_join(swissprot3, trembl3) # all swissprot are unique (30444)

test2 <- merge(trembl3, ncbi_blast3, by = "SeqName") # 24359 from trembl

test3 <- 25513 - 24359 # 1154 from NCBI 

all <- 30444 + 24359 + 1154
none <- 64636 - all


allper <- (all/64636)*100 #86.6%

nonper <- 100 - (all/64636)*100 #13.4%, 

SPper <- (30444 / 64636)*100 #47.1%
TRper <- (24359 / 64636)*100 #37.7%
NCper <- (1154 / 64636)*100 #1.8%



