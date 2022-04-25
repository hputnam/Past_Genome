---
layout: page
title: Genome annotation with MAKER
---

### What is MAKER?

MAKER is a genome annotation pipeline from the [Yandell lab](http://www.yandell-lab.org/software/maker.html), [Campbell et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4286374/). 
As described on the page "Its purpose is to allow smaller eukaryotic and prokaryotic genome projects to independently annotate their genomes and to create genome databases. MAKER identifies repeats, aligns ESTs and proteins to a genome, produces ab-initio gene predictions and automatically synthesizes these data into gene annotations having evidence-based quality values."  Since it makes ab initio gene predictions, training ab-initio predictors is an important step to produce good results.

A tutorial is available with details and explanations at [MAKER Tutorial](http://gmod.org/wiki/MAKER_Tutorial).

To annotate a genome using Maker, you need the following files:

- The genome sequence in fasta format
- Assembled RNA-seq transcripts of the species in fasta format (TO DO) 
- Protein sequences in fasta format, usually from closely related species or from a curated sequence database like UniProt/SwissProt. (TO DO)
Maker will align the transcript and protein sequences on the genome sequence to determine gene positions.

### Assessing genome quality 
Before running the full annotation process, we need first to evaluate the quality of the genome assembly. This can be done using [BUSCO](https://busco.ezlab.org/) which gives a completeness score. (TO DO)

### MAKER control files
The first thing is to create new control files for your MAKER run. The control files tell MAKER how to run and where to find additional software such as different gene callers.

In the desired folder for your MAKER run type:

`maker -CTL`

This will create three files:

- **maker_bopts.ctl** containing settings for BLAST and Exonerate.
- **masker_exe.ctl** with all the paths to different executables used by MAKER on your system.
- **maker_opts.ctl** is the file controlling MAKERs running behavior. 

These files contain paths to files that are used by `maker` and our choices for analysis. Edit the maker_opts.ctl file to specify the genome assembly sequence, experimental alignment evidence and which gene finding method to use.  
