---
title: "R Notebook"
output: html_notebook
---

# Introduction

The purpose of this Markdown is to document variant processing and QC of the 470k exome release of the UK Biobank data.

# Table of Contents

- [Setup](#setup)
- [Utility Functions](#utility-functions)
  * [Cloud Workstation](#cloud-workstation)
  * [Querying jobs by status](#querying-jobs-by-status)
  * [Checking Finished Jobs](#checking-finished-jobs)
- [Variant Processing and Filtering](#variant-processing-and-filtering)
  * [1. Splitting BCFs](#1-splitting-bcfs)
  * [2. Filtering and VEP Annotation](#2-filtering-and-vep-annotation)
  * [3. Create BGEN from BCF](#3-create-bgen-from-bcf)
    + [3a. Make BGEN list file](#3a-make-bgen-list-file)
    + [3b. Run make-bgen](#3b-run-make-bgen)
  * [4. Collapse Variants](#4-collapse-variants)
    + [4a. Generating Genetic File List](#4a-generating-genetic-file-list)
- [4b. Running Collapse Variants](#4b-running-collapse-variants)
  * [5. Run Association Testing](#5-run-association-testing)
    + [5a. Setting Genetic Data:](#5a-setting-genetic-data-)
      - [WES Samples List](#wes-samples-list)
      - [Define WBA grouping](#define-wba-grouping)
      - [Running mrcepid-buildgrms](#running-mrcepid-buildgrms)
    + [5b. Setting Covariates](#5b-setting-covariates)
    + [5c. Post-processing of statistics](#5c-post-processing-of-statistics)
    + [5d. Generate Transcript List](#5d-generate-transcript-list)
      - [Download and Input Initial Transcript Lists](#download-and-input-initial-transcript-lists)
      - [Setting Relative Position of Genes in the Genome](#setting-relative-position-of-genes-in-the-genome)
      - [Get per-gene Counts](#get-per-gene-counts)
      - [Identifying Genes with 0 Variants](#identifying-genes-with-0-variants)
      - [Sample Counts](#sample-counts)
      - [Writing Transcript File for DNANexus](#writing-transcript-file-for-dnanexus)
    + [5f. Running Associations](#5f-running-associations)
      - [Burden Tests](#burden-tests)
      - [Individual Gene(s)](#individual-gene-s-)
      - [PheWAS](#phewas)
- [Loading Association Tests](#loading-association-tests)
  * [Functions for Loading Data](#functions-for-loading-data)
  * [Functions for Plotting](#functions-for-plotting)
  * [Example Code For Loading and Plotting](#example-code-for-loading-and-plotting)



# Setup

```{r setup}

library(data.table)
library(tidyverse)
library(patchwork)
library(lemon)
library(readxl)
library(broom)

theme <- theme(panel.background=element_rect(fill="white"),line=element_line(size=1,colour="black",lineend="round"),axis.line=element_line(size=1),text=element_text(size=16,face="bold",colour="black"),axis.text=element_text(colour="black"),axis.ticks=element_line(size=1,colour="black"),axis.ticks.length=unit(.1,"cm"),strip.background=element_rect(fill="white"),axis.text.x=element_text(angle=45,hjust=1),legend.position="blank",panel.grid.major=element_line(colour="grey",size=0.5),legend.key=element_blank())

## Default theme w/legend
theme.legend <- theme + theme(legend.position="right")

```

# Utility Functions

This section contains simple examples of utility functions for working on the UK Biobank Research Access Platform (RAP). You may need to refine them yourself.

## Cloud Workstation

This will open a cloud work station with 2 CPUs, 4 Gb of memory and 100 Gb of storage. If more is needed, modify the `instance-type` parameter.

The -isnapshot parameter can also be used to load an image with some basic tools installed:

-   bcftools
-   samtools
-   tabix/bgzip
-   plink
-   plink2
-   qctool
-   bgenix

I have generated such an image within several projects, but it will be up to you to do the same.

```{bash, eval = F}

dx run app-cloud_workstation --ssh --instance-type mem1_ssd1_v2_x4 -imax_session_length=8h -isnapshot=project-GFPBQv8J0zVvBX9XJyqyqYz1:July_28_2022_15_06.snapshot

```

## Querying jobs by status

```{bash, eval = F}
dx find jobs --name 'BCFSplitter' --created-after "2022-08-04" --state 'failed' -n 100 --json --verbose > failed_splitter.json
```

## Checking Finished Jobs

Generalised workflow for checking finished VCF/BCFs and launching missed jobs

```{bash, eval = F}

# Gets a list of all bcfs that need to be processed and gets the root name
# Should be 4685 files long
dx ls -l filtered_vcfs/ukb23148_c*_b*_v1_chunk?.norm.filtered.tagged.missingness_filtered.annotated.cadd.bcf | perl -ane 'chomp $_; if ($F[5] =~ /(ukb23148_c[0-9XY]{1,2}_b\d+_v1_chunk\d)\./) {$root = $1; if ($F[6] =~ /^\((\S+)\)$/) {$id = $1; print "$root\t$id\n"}}' > vcfs.txt

# Gets a list of all finished bcfs (and file IDs) -- this will need to  be modified:
dx ls -l collapsed_vcfs/ukb23148_c*_b*_v1_chunk*.norm.filtered.tagged.missingness_filtered.annotated.cadd.MISS.tar.gz | perl -ane 'chomp $_; if ($F[5] =~ /(ukb23148_c[0-9XY]{1,2}_b\d+_v1_chunk\d)\./) {print "$1\n"}' > finished.txt

# Use matcher.pl to extract files that still need to be finished
~/Documents/Current\ Projects/Fertility/UKBBFertility/scripts/matcher.pl -file1 finished.txt -file2 vcfs.txt | perl -ane 'chomp $_; print "$F[1]\n";' > incomplete.txt

# Then just redo split / dx upload as below

```

# Variant Processing and Filtering

To perform variant processing and qualite control, I will be using the pipeline documented here: <https://github.com/mrcepid-rap#vcf-filtering-and-rare-variant-burden-testing>. Please see this GitHub repository for more details as to how this pipeline works. Here I am just documenting commands used to run this pipeline for the purposes of reproducibility. All work for this project is being completed on the UK Biobank Research Access Platform (RAP) project MRC - Variant Filtering 450k (project-G6BJF50JJv8p4PjGB9yy7YQ2). Analysis will procede in 6 steps (after '-' is the applet being used on the RAP):

1.  Splitting BCFs - [mrcepid-bcfsplitter](http://github.com/mrcepid-rap/mrcepid-bcfsplitter)
2.  Filtering and VEP annotation - [mrcepid-filterbcf](https://github.com/mrcepid-rap/mrcepid-filterbcf)
3.  Make combined .bgen format files and annotations for each chromosome - [mrcepid-makebgen](https://github.com/mrcepid-rap/mrcepid-makebgen)
4.  Collapse variants into resource bundles for association testing - [mrcepid-collapsevariants](https://github.com/mrcepid-rap/mrcepid-collapsevariants)
5.  Association testing - [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting)

There is also one supplementary step that prepares genetic resource files and sample lists for association testing:

-   Building of GRMs and sample inclusion/exclusion lists - [mrcepid-buildgrms](https://github.com/mrcepid-rap/mrcepid-buildgrms)

## 1. Splitting BCFs

The RAP stores data in 977 separate vcf.gz files in the folder `/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - interim 450k release/`. All files come pre-indexed with a .tbi index file.

Due to the size of these files, the first analysis step requires splitting these files into separate chunks. In brief, each vcf file is split into chunks equal to:

$$ y = sgn(x)[|\frac{n.var_{vcf}}{5000}| + 0.5] $$

Thus, a split bcf file can have no fewer than 2500 sites and no more than 7500 sites. In practice, this results in most original vcf.gz files being split into between 4-5 smaller .bcfs. Note that this is the *ONLY* part of the pipeline that I have not parallelized. This was because I did this analysis prior to deciding to parallelize all other parts of the codebase.

**BIG Note** There is a small bug in this pipeline that will duplicate a site if the following are true:

1.  The site has an identical position value to another site in the vcf (i.e. split multiallelics)
2.  The site and it's position duplicate SPAN the junction of a bcf split
    -   e.g. if variant 1 of a multiallelic pair is variant no. 5000 and variant 2 is variant no. 5001

This bug will result in both variants in the pair being present in both split bcf files. In theory, this could also happen to a 3 variant multiallelic, but I have not encountered such a situation. I account for this error later in the pipeline when I extract variants according to filtering criteria (mrcepid-collapsevariants). This case is also rare. Across the 47 vcf.gz files that cover chromosome 7, this situation happened once.

```{bash, eval = F}

# 1. get a list of the initial VCF
dx ls -l 'Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23157*.vcf.gz' | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {$id = $1; $F[5] =~ s/.vcf.gz//; print "$F[5]\t$id\n";}' > bcf_list.txt 

# 2. Split to a manageable per-size job:
split -l 15 -d ../bcfsplitter_input.lst splitjob_

# 3. Upload to DNAnexus
dx upload --brief -r --destination batch_lists/ batch_lists/*

# 4. Launch jobs:
dx ls -l batch_lists/ | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-bcfsplitter --priority low --yes --brief --destination split_vcfs/ -iinput_vcfs=$1\n";}' | less
```

## 2. Filtering and VEP Annotation

This step performs filtering and VEP annotation according to parameters decided by the MRC Epidemiology RAP working group. Please see the [github page](https://github.com/mrcepid-rap/mrcepid-filterbcf) for more information.

**BIG NOTE:** There are default file locations for all of the annotations. If you do not have access to the internal MRC project(s) you will need to provide locations to all annotations files as inputs. The below command points to a .json that does not exist byt default with this repo, but gives an indication of how these files can be provided in bulk.

All split vcf files are stored in `/filtered_vcfs/` and match a regex like `ukb23148_c*_b*_v1_chunk?.bcf`. There are a total of 4,722 split bcfs.

```{bash, eval = F}

# Get a list of all split bcfs.
# The total number of files being processed here is 4722.
dx ls -l filtered_vcfs/ukb23148_c*_b*_v1_chunk?.bcf | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "$1\n";}' > bcf_list.txt 

# Split this file list into individual "jobs":
# This will generate 78 files with 60 BCFs per file. This means that each core will be ~used 2x in each job
split -a 3 -d -l 14 bcf_list.txt bcf_list_

# Upload these lists onto DNA Nexus:
dx upload --brief --no-progress bcf_list_* --destination batch_lists/

# Generate a list of jobs to run on DNA Nexus, and submit
dx ls -l batch_lists/ | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {$file = $1; if ($F[5] =~ /bcf_list_(\d+)/) {$coords = "coordinates_$1" . ".tsv"; print "dx run mrcepid-filterbcf --destination filtered_annotated_vcfs/ --priority normal -f jsons/filterbcf.json --yes --brief -iinput_vcfs=$file -icoordinates_name=$coords\n";}}' | less

```

## 3. Create BGEN from BCF

### 3a. Make BGEN list file

Need to login to a cloud workstation (see above) and download all of the coordinates files from each original job in step (2). I then make a single .tsv.gz file listing the file location of all filtered & annotated bcf/vep in the following format:

    #chrom    start   stop    fileprefix    bcf_dxpy    vep_dxpy
    chr1    1   1000000   ukb23157_c1_b1_v1_chunk1    file-1234567890ABCDEFGHIJK    file-0987654321KJIHGFEDCBA
    chr1    1000001   2000000   ukb23157_c1_b1_v1_chunk2    file-ABCDEFGHIJK1234567890    file-KJIHGFEDCBA0987654321

Commands below make this file on the cloud workstation:

```{bash make bgen list, eval = F}

dx-su-contrib
mkdir coords & cd coords
dx download -a --no-progress --lightweight filtered_annotated_vcfs/coordinates_*.tsv

## Do some quick QC to make sure we have 4722 UNIQUE files:
# Should give 4722 as output
cat *.tsv | awk '{print $4}' | sort | uniq | wc -l 

# Make a header file...
echo -e '#chrom\tstart\tstop\tfileprefix\tbcf_dxpy\tvep_dxpy' > coord_header.tsv
# Concatenate everything together:
cat coord_header.tsv > bcf_coordinates.tsv; cat *.tsv | sort -k 1,1 -k 2,2n coordinates_*.tsv >> bcf_coordinates.tsv
bgzip 

# Upload
dx upload --destination filtered_annotated_vcfs/ bcf_coordinates.tsv.gz

```

### 3b. Run make-bgen

```{bash make bgen, eval = F}

perl -e 'for my $chr (1..22,"X","Y") {print "dx run mrcepid-makebgen --destination filtered_bgen/ --priority normal -ichromosome=\"chr$chr\" -icoordinate_file=file-GFx0z70J0zVpJzVF0QG6Fjfk;\n";}' | bash

```

## 4. Collapse Variants

### 4a. Generating Genetic File List

The input `bgen_index` to mrcepid-runassociationtesting is a tab-delimited file that looks like the following:

```
chrom   vep   bgen    index   sample
chr1    file-1234567890ABCDEFGH   file-0987654321ABCDEFGH   file-1234567890HGFEDCBA   file-0987654321HGFEDCBA
```

The following script should automate the production of this file:

`scripts/build_bgen_list.py`

Remember to upload the file to DNANexus like:

`dx upload --destination filtered_bgen/ bgen_locs.tsv`

# 4b. Running Collapse Variants

Here we generate a list of variant types that we want to collapse on to generate sets of variants for associationtesting. Obviously these can be tweaked, but here I am generating a set of initial variants useful for most scenarios.

**WARNING:** DO NOT enable CADD below unless you are an academic or have a licence!

```{r collapse variants, eval = F}

MAFs <- list('MAF_01' = "AF<0.001",
             'MAF_1' = "AF<0.01",
             'AC_1' = "AC==1")
CSQs <- list('HC_PTV' = 'PARSED_CSQ=="PTV" & LOFTEE=="HC"',
            'PTV' = 'PARSED_CSQ=="PTV"',
            'MISS' = 'PARSED_CSQ=="MISSENSE"',
            #'MISS_CADD25' = 'PARSED_CSQ=="MISSENSE" & CADD>=25',
            'MISS_REVEL0_5' = 'PARSED_CSQ=="MISSENSE" & REVEL>=0.5',
            'MISS_REVEL0_7' = 'PARSED_CSQ=="MISSENSE" & REVEL>=0.7',
            #'MISS_CADD25_REVEL0_7' = 'PARSED_CSQ=="MISSENSE" & CADD>=25 & REVEL>=0.7',
            'SYN' = 'PARSED_CSQ=="SYN"',
            'DMG' = '((PARSED_CSQ=="PTV" & LOFTEE=="HC") | (PARSED_CSQ=="MISSENSE" & REVEL>=0.5))')

for (maf in names(MAFs)) {
  
  for (csq in names(CSQs)) {
    
    base_command <- 'dx run mrcepid-collapsevariants --priority normal --destination collapsed_variants --brief --yes -ibgen_index=file-GFzfyyQJ0zVqXZQBK4592QQg'
    expression <- paste0("'",paste(MAFs[maf], CSQs[csq], 'FILTER=="PASS"', sep = " & "), "'")
    filename <- paste(csq, maf, sep = '-')
    cat(paste(base_command, # standard input for all runs
              paste0('-ifiltering_expression=',expression), # filtering expression
              paste0('-ifile_prefix=', filename), # file name
              sep = ' '))
    cat('\n')
    
    
  }
  
  cat('\n')
  
}

# Note -- It's normally easier to add this to a shell script and submit rather than copy-pasting
```

## 5. Run Association Testing

mrcepid-runassociation testing requires several pre-built files (other than variant definitions) to be able to run:

A.  Genetic Data - all produced by mrcepid-buildgrms
    -   Exclusion / Inclusion list - this can be modified according to user need, but default lists are made here
    -   QCd array data - in plink2 bed/bim/bam format
    -   List of low MAC variants in the QCd array data - for filtering with some associationtesting tools
    -   sparse GRM - generated from the prebuilt GRM from Bycroft et al.

B. Base Covariates - A tab-delimited file of base covariates that all models will include
    - PC1..10
    - age
    - sex
    - wes_batch

C. Generating statistics files for quality control
    - Per-individual stats
    - All VEP annotations

D. Transcripts List - A list of valid ENST/MANE transcripts with annotations

E. File listing locations of genetic data

### 5a. Setting Genetic Data:

Building the set of files in item (1) above, also requires some pre-built data. We document how to generate these files below.

#### WES Samples List

This needs to be done on a cloud workstation

```{bash, eval = F}

# Download (Quotes ["] are required due to spaces)
dx mkdir project_resources/wes_resources/
dx download "Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23157_cY_b0_v1.vcf.gz"
# bcftools query -l to get a list of samples. There SHOULD be 469,835 individuals
bcftools query -l ukb23157_cY_b0_v1.vcf.gz > wes_samples.470k.txt 
dx upload --destination project_resources/wes_resources/ wes_samples.470k.txt

```

#### Define WBA grouping

The file `pca_participant.tsv` is generated by using the phenotype browser on the DNANexus website. Relevant fields were added manually and include:

-   PC1-4 (field 22009-0.[1-4])
-   Genetic ancestry grouping (field 22006-0.0)
-   Self-provided ancestry grouping (field 21000-0.0)

This file was then queried locally. All we have done is calculate a kmeans clustering into 4 groups (general human super-population structure), and taken the 'mode' group as European ancestry (because most participants are of European genetic ancestry). We then query this against the self-reported ancestry grouping (field 21000) to ensure consistency between ours and participant assessment

Individuals that fall within this group are included in the 'European' ancestry definition for the purposes of association testing. These resulting file (`wba.txt`) is then uploaded onto DNANexus into the `project_resources/genetics/` folder for use with the `mrcepid-makegrms` applet.

```{r pca}

pca <- fread("data_files/ancestry/pca_participant.tsv", colClasses = c("eid"="character"))

# Load ancestry codings for easier labelling
ancestry_codings <- fread("data_files/ancestry/coding1001.tsv")
ancestry_codings[,parent_id:=if_else(coding < 10, coding, parent_id)]
setkey(ancestry_codings,"parent_id")
ancestry_codings[,meaning:=factor(meaning,levels=ancestry_codings[,meaning])]

cols <- c('eid', paste('22009-0',seq(1,4), sep = '.'), '22006-0.0', '21000-0.0')
setnames(pca, cols, c('n_eid', paste0('PC', seq(1,4)), 'wba', 'ancestry'))

pca <- merge(pca,ancestry_codings[,c("coding","meaning")],by.x="ancestry",by.y="coding",all.x=T)
setcolorder(pca,c("n_eid","ancestry","meaning","wba","PC1","PC2","PC3","PC4"))

pca <- pca[!is.na(PC1)]
pca[,wba:=if_else(is.na(wba), 0, 1)]

cluster_data <- pca[,c("PC1","PC2","PC3","PC4")]
cluster_data <- as.matrix(cluster_data)
rownames(cluster_data) <- pca[,n_eid]

set.seed(6547)

pca_cluster <- kmeans(cluster_data, 4, iter.max = 10, nstart = 1, algorithm = c("Hartigan-Wong"))
cluster_table <- data.table("n_eid" = names(pca_cluster$cluster), "cluster" = as.integer(pca_cluster$cluster))

pca <- merge(pca, cluster_table, by = "n_eid")

european_ancestries <- c(1, 1001, 1002, 1003, 6, -3, -1)
south_asian_ancestries <- c(3, 3001, 3002, 3003, 3004, 6, -3, -1)
african_ancestries <- c(4, 4001, 4002, 4003, 6, -3, -1)

# Note that you may have to set these cluster numbers yourself based on clustering results!
# This prioritises unknowns to the first cluster checked. I.e. if an individual with unknown (6,-3,-1) is in cluster 2, 
# they will go to the EUR category, then SAS, then AFR
pca[,"ancestry":=if_else(cluster == 2 & (ancestry %in% european_ancestries | is.na(ancestry)), "eur", 
                         if_else(cluster == 4 & (ancestry %in% south_asian_ancestries | is.na(ancestry)), "sas",
                                 if_else(cluster == 1 & (ancestry %in% african_ancestries | is.na(ancestry)), "afr", 
                                 as.character(NA))))]

pca[,table(meaning,cluster)]
pca[,table(meaning,ancestry,useNA = "always")]

ggplot(pca, aes(PC1, PC2, colour = as.factor(cluster))) + geom_point(size = 0.5) + theme.legend
ggplot(pca, aes(PC3, PC4, colour = as.factor(cluster))) + geom_point(size = 0.5) + theme.legend

fwrite(pca[,c("n_eid","ancestry")], 'data_files/ancestry/ancestry.txt', quote = F, sep = ' ', row.names = F, col.names = T, na = "NA")
```

#### Running mrcepid-buildgrms

This will generate the required files for running burden tests. Replace the inputs with the DNANexus file IDs generated above. Also, make sure to create a folder for this data to be deposited in if you haven't already (hard-coded to project_resource/genetics/, below).

```{bash}

dx run mrcepid-buildgrms --destination project_resources/genetics_v2/ -isample_ids_file=file-GFjJ0yjJ0zVb9BK82ZQFkx0y -iancestry_file=file-GGbPpqQJ0zVf9F0y4fXbyqxJ --priority=normal

```

### 5b. Setting Covariates

Base covariates were processed on in a RAP python notebook environment (see [here](https://documentation.dnanexus.com/user/jupyter-notebooks)). A python notebook to repeat this analysis is available as part of this repository:

`python_notebooks/base_covariates.ipynb`

### 5c. Post-processing of statistics

Following this step I launched a [cloud workstation](#cloud-workstation) to generate final per-gene and per-individual statistics and generate a list file of variant coordinates within each split and annotated bcf-file for future users. The [cloud workstation](#cloud-workstation) was launched as described above. 

Note that this code requires the use of the 'mash.pl' script, which is included in the `scripts/` folder. This will need to be uploaded to the RAP in order to use as directed below when downloading (`dx download mash.pl`)

Following launch, I used the following commands:

```{bash post-process stats, eval = F}

# Download all statistics files from the 'annotate-cadd' step into two separate directories:
mkdir per_indv/
cd per_indv/ \
  && dx download --no-progress filtered_annotated_vcfs/ukb*_c*_b*_v1_chunk*.per_indv.tsv

# Download required processing files
dx download project_resources/wes_resources/wes_samples.470k.txt
dx download mash.pl

# Merge per-individual statistics files together:
ls per_indv/*.per_indv.tsv > per_indv_list.txt
./mash.pl > 470k_indv_counts.tsv
bgzip 470k_indv_counts.tsv
dx upload --destination project_resources/wes_resources/ 470k_indv_counts.tsv.gz

# Merge per-variant statistics files together:
mkdir vep/ && cd vep/ && dx download --no-progress filtered_bgen/*.vep.tsv.gz
zgrep -hv CHROM *.tsv.gz > 470k_vep.tsv
zgrep CHROM chrY.filtered.vep.tsv.gz > header.tsv
sort -k 1,1 -k2,2n 470k_vep.tsv > 470k_vep.sorted.tsv
cat header.tsv 470k_vep.sorted.tsv > 470k_vep.sorted.header.tsv
mv 470k_vep.sorted.header.tsv 470k_vep.sorted.tsv
bgzip 470k_vep.sorted.tsv
tabix -c 'C' -s 1 -b 2 -e 2 470k_vep.sorted.tsv.gz
```

### 5d. Generate Transcript List

#### Download and Input Initial Transcript Lists

To get a reasonable set of protein-coding transcripts for Hg38 I generated two separate files:

1.  `data_files/hgTables.txt.gz` - File downloaded from the UCSC Genome Browser - Table Browser. The purpose of this file is to calculate total coding sequence length, which cannot be derived from BioMart (file #2). I selected the latest GENCODE (v41 currently) and selected the following rows:

    -   Transcript name (name)
    -   Gene Start (chromStart)
    -   Coding sequence start (thickStart)
    -   Coding sequence end (thickEnd)
    -   Number of exons (blockCount)
    -   Exon start coordinates (chromStarts)
    -   Exon lengths (blockSizes)
    -   transcript type (e.g. protein_coding; transcriptType)

This file is then processed like: `python3 scripts/get_cds_len.py`

2.  `data_files/mart_export.txt` - File downloaded from ENSEMBL BioMart. The purpose of this file is to get valid ENSEMBL transcripts for all possible genes in Hg38, with canonical and MANE transcript for each gene. I selected the latest version and GRCh38 and selected the following rows:

    -   Gene stable ID (form of ENSG)
    -   Transcript stable ID (form of ENST)
    -   MANE Select ID (form of NM)
    -   Transcript Length - Note that this is **NOT** coding sequence length, as above
    -   HGNC Symbol
    -   Is ENSEMBL canonical transcript (1 for yes, blank for no)
    -   Chromosome
    -   Start coordinate
    -   End coordinate
    -   Transcript type (e.g. protein_coding)

I then use these files to perform the following to get set of valid transcripts for all downstream statistics:

```{r import transcripts}

transcripts <- fread("data_files/mart_export.txt")
setnames(transcripts,names(transcripts),c("ENSG","ENST","MANE","transcript_length","SYMBOL","CANONICAL","chrom","start","end","BIOTYPE"))
transcripts[,CANONICAL:=if_else(is.na(CANONICAL),F,T)]
transcripts[,chrom:=as.character(chrom)]
transcripts <- transcripts[(MANE != "" | CANONICAL == T) & grepl("^[\\dXY]{1,2}$", chrom, perl = T) & BIOTYPE == "protein_coding"]
transcripts[,dummy:=1]

# Get real CDS length
cds_table <- fread("data_files/cds_lengths.txt")
setnames(cds_table,names(cds_table),c("ENST","cds_length"))
# Remove duplicated PAR transcripts (same cds_len on X and Y)
cds_table <- unique(cds_table)

# merge
transcripts <- merge(transcripts, cds_table, by = "ENST", all.x = T)

# set a relative location in the human genome for all transcripts for manhattan plots:
get_manh_pos <- function(start, c) {
  
  c <- paste0('chr',c)
  cuml_length <- hg38_index[chrom == c, cuml_length]
  total_length <- hg38_index[chrom == 'chrY', length + cuml_length]
  manh_pos <- (start + cuml_length) / total_length
  return(manh_pos)
  
}
hg38_index <- fread("data_files/hs38DH.fa.fai", col.names = c("chrom","length","cuml_length","linecount","linewidth"))

transcripts[,manh.pos:=get_manh_pos(start, chrom), by = 1:nrow(transcripts)]

```

#### Setting Relative Position of Genes in the Genome

This is just to enable easy plotting on Manhattan plots.

```{r mean gene position}

# Get chromosome locations for the plot
mean.chr.pos <- transcripts[,mean(manh.pos), by = chrom]

```

#### Get per-gene Counts

Note that the vep.tsv.gz file cannot be included in this repository, so you will have to place the file (once generated) into the correct folder as indicated below

```{r fig.height=6, fig.width=8}

vep.stats <- fread("data_files/annotations/470k_vep.sorted.tsv.gz")
vep.stats[,dummy:=1]
counts <- vep.stats[PARSED_CSQ=="SYN" & FILTER == "PASS" & MAF < 1e-3,sum(dummy),by=c("ENST")]

transcripts.annotated <- merge(counts,transcripts,by="ENST",all.y=T)
transcripts.annotated[,syn.count:=if_else(is.na(V1),0,V1)]
transcripts.annotated[,V1:=NULL]
transcripts.annotated[,coord:=paste0("chr", chrom,":",start,"-",end)]

# Highlight genes with 0 syn variants on the plot
transcripts.annotated[,highlight:=if_else(syn.count == 0, T, F)]

ggplot(transcripts.annotated,aes(cds_length,syn.count,colour=highlight)) + geom_point(size = 0.5) + xlim(0,30000) + ylim(0,3100) + geom_abline(slope = 0.08062,intercept=13.25662) + theme
```

#### Identifying Genes with 0 Variants

This will pull out all genes with "0" synonymous variants and annotate them according to preset rules based on manual inspection.

```{r fig.height=5, fig.width=10}

missing <- transcripts.annotated[syn.count == 0]
fwrite(transcripts.annotated[syn.count == 0], "data_files/transcripts/missing_genes.tsv", col.names = F, row.names = F, sep = "\t", quote = F)

# Checking number of baits for missing genes
# Need to bgzip/tabix xgen_plus_spikein.GRCh38.bed
# perl -ane 'chomp $_; print "echo $F[0]\ntabix xgen_plus_spikein.GRCh38.bed.gz \"$F[14]\" | wc -l;\n";' missing_genes.tsv | bash > bait_counts.txt
# perl -ane 'chomp $_; if ($_ =~ /ENST/) {print "$_\t"} elsif ($_ =~ /(\d+)/) {print "$1\n"}' bait_counts.txt > bait_counts.form.txt
baits <- fread("data_files/transcripts/bait_counts.form.txt", col.names = c("ENST","bait.count"))
missing <- merge(missing,baits,by="ENST")

# Checking number of IDd variants for missing genes
# perl -ane 'chomp $_; print "echo $F[0]\ntabix ../470k_vep.sorted.tsv.gz \"$F[14]\" | wc -l;\n";' missing_genes.tsv | bash > variant_counts.txt
# perl -ane 'chomp $_; if ($_ =~ /ENST/) {print "$_\t"} elsif ($_ =~ /(\d+)/) {print "$1\n"}' variant_counts.txt > variant_counts.form.txt
variants <- fread("data_files/transcripts/variant_counts.form.txt", col.names=c("ENST","variant.count"))
missing <- merge(missing,variants, by = "ENST")

# Checking proportion of variants in a given gene with FAIL for FILTER
fail.prop <- data.table(pivot_wider(vep.stats[,sum(dummy),by = c("ENST","FILTER")],id_cols = "ENST", names_from = "FILTER", values_from = "V1", values_fill = 0))
fail.prop[,fail.prop:=FAIL / (FAIL + PASS)]
missing <- merge(missing, fail.prop[,c("ENST", "fail.prop")], all.x = T, by = "ENST")

# Setting fail categories
missing[,fail.cat:=if_else(bait.count == 0, "NOT_SEQ", "")]
missing[,fail.cat:=if_else(grepl("\\-\\S{2}",SYMBOL,perl = T) & fail.cat == "", "READTHROUGH", fail.cat)]
missing[,fail.cat:=if_else(chrom == "Y" & fail.cat == "", "CHR_Y", fail.cat)]
missing[,fail.cat:=if_else(SYMBOL == "" & fail.cat == "", "NO_NAME", fail.cat)]
missing[,fail.cat:=if_else(fail.prop > 0.75 & fail.cat == "" & !is.na(fail.prop), "NO_PASS_VARS", fail.cat)]
missing[,fail.cat:=if_else(variant.count == 0 & fail.cat == "", "NO_VARS_IN_VCF",fail.cat)]
missing[,fail.cat:=if_else(fail.cat == "", "OTHER",fail.cat)]

# Merge back into the main plot
transcripts.annotated <- merge(transcripts.annotated,missing[,c("ENST","fail.cat")],by="ENST",all.x=T)
transcripts.annotated[,dummy:=1]
transcripts.annotated[,fail:=if_else(is.na(fail.cat),F,T)]
per.chr.fail <- data.table(pivot_wider(transcripts.annotated[,sum(dummy),by=c("fail","chrom")], names_from = fail, values_from = V1))
per.chr.fail[,chrom:=factor(chrom,levels=c(1:22,"X","Y"))]
setkey(per.chr.fail,chrom)
per.chr.fail[,prop.fail:=(`TRUE`/(`TRUE` + `FALSE`))*100]

ggplot(per.chr.fail,aes(chrom,prop.fail)) + geom_col() + xlab("Chromosome") + ylab("Proportion Failed Genes") + theme + theme(panel.grid.major.x=element_blank())

per.cat.fail <- transcripts.annotated[!is.na(fail.cat),sum(dummy),by=c("fail.cat","chrom")]
per.cat.fail[,chrom:=factor(chrom,levels=c(1:22,"X","Y"))]
ggplot(per.cat.fail,aes(chrom, V1, fill=fail.cat, group=fail.cat)) + geom_col() + theme.legend + theme(panel.grid.major.x=element_blank())
```

#### Sample Counts

```{r sample counts}

wba <- fread("data_files/ancestry/wba.txt")
setnames(wba,"n_eid","eid")
wba[,eid:=as.character(eid)]

indv.counts <- fread("data_files/470k_indv_counts.tsv.gz")
indv.counts[,eid:=as.character(eid)]

indv.counts <- merge(indv.counts,wba[,c("eid","European_ancestry")],by="eid")

ggplot(indv.counts,aes(SYN)) + geom_histogram(binwidth=1) + geom_vline(xintercept=80) + theme

ggplot(indv.counts[European_ancestry == 1],aes(PTV)) + geom_histogram(binwidth=1) + theme
ggplot(indv.counts[European_ancestry == 1],aes(MISSENSE)) + geom_histogram(binwidth=1) + theme
ggplot(indv.counts[European_ancestry == 1],aes(SYN)) + geom_histogram(binwidth=1) + xlim(0,160) + theme

quantile(rpois(n = 469430, lambda = 80), probs = c(0.997))
indv.counts[,high.SYN:=if_else(SYN>quantile(rpois(n = 453342, lambda = 80), probs = c(0.995))[[1]],1,0)]
table(indv.counts[,c("high.SYN","European_ancestry")])
```

#### Writing Transcript File for DNANexus

Finally, I write a file that is uploaded to DNANexus for annotation during Association Testing.

```{r set transcripts}

# Set specific column order to enable tab-indexing
transcripts.annotated <- transcripts.annotated[,c("chrom","start","end","ENST","ENSG","MANE","transcript_length","SYMBOL","CANONICAL","BIOTYPE","cds_length","syn.count","coord","fail.cat","fail","manh.pos")]

# Sort:
setkeyv(transcripts.annotated,cols=c("chrom","start","end"))

# And write a file
fwrite(transcripts.annotated,"transcripts.tsv",quote=F,row.names=F,sep="\t")

```

Make sure to bzip, and tabix after...

```
bgzip transcripts.tsv
tabix -c c -s 1 -b 2 -e 3 transcripts.tsv.gz
dx upload --destination project_resources/wes_resources/ transcripts.tsv.gz*
```

### 5f. Running Associations

The following sections document different use-cases for the mrcepid-runassociationtesting app/applet.

#### Burden Tests

This uses the `launch.sh` script located in `./scripts/`. See that file for more information on how association testing is run for individual tools. This script will launch a given tool on all masks for a given phenotype. The input parameters are:

1. tarballs
2. phenofile(s)
3. Phenoname (can be 'null')
4. is binary? (must be true/false)
5. Output prefix
6. Sex (0/1/2)
7. Additional covariates file (can be 'null')
8. Additional quantitative covariates in additional covariates file (can be 'null')
9. Additional categorical covariates in additional covariates file (can be 'null')
10. Tool to use (bolt/glm/saige/regenie/staar). Can be 'null' to run all 5 tools or a comma-separated list (e.g. bolt,glm)

**NOTE**: To use this script, you will need to replace the default file settings provided as part of the "DEFAULT_COVARS" variable in the `./launch.sh` script.

```{bash, eval = F}

# Make a list of all association test tar files for masks with MAF < 0.1%:
dx ls -l collapsed_variants/*.tar.gz |  grep 'MAF_01' | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "$1\n";}' > variant_mask_list.txt
dx upload variant_mask_list.txt --destination tarball_lists/

# This will launch a single burden test using BOLT. Additional tools can be run by modifying the final parameter (e.g. bolt,saige,staar)
./launch.sh file-GGJZB4jJ0zVbJv8J5K8z2KF9 file-GGf1gF8JJy85qqGk4x4GQ22F null false menopause_34_eur 0 null null null bolt file-GGbZFKQJP7J5223J4v6k49P3 null # EUR Related
./launch.sh file-GGJZB4jJ0zVbJv8J5K8z2KF9 file-GGf1gF8JJy85qqGk4x4GQ22F null false menopause_34_sas 0 null null null bolt file-GGbZFPjJP7J0vZ6f4v65G4Zj null # SAS Related
./launch.sh file-GGJZB4jJ0zVbJv8J5K8z2KF9 file-GGf1gF8JJy85qqGk4x4GQ22F null false menopause_34_afr 0 null null null bolt file-GGbZFQjJP7J92Xq94x9bXJXG null # AFR Related
./launch.sh file-GGJZB4jJ0zVbJv8J5K8z2KF9 file-GGf1gF8JJy85qqGk4x4GQ22F null false menopause_34_new 0 null null null bolt file-GGbZFJ8JP7JKbG8v4xkQVGY7 file-GGbQ3yQJ0zVxXg8y4yY2ByXZ # New (470k samples) only
./launch.sh file-GGJZB4jJ0zVbJv8J5K8z2KF9 file-GGf1gF8JJy85qqGk4x4GQ22F null false menopause_34_all 0 null null null bolt file-GGbZFJ8JP7JKbG8v4xkQVGY7 null # All (Eur/Sas/Afr/Admixed) Related

```

#### Individual Gene(s)

```{bash}

dx run mrcepid-runassociationtesting --instance-type mem3_ssd1_v2_x16 --priority normal --destination results/ -imode='extract' -ioutput_prefix="menopause" -iinput_args='--gene_ids BRCA2 CHEK2 ZNF518A --association_tarballs file-GGJZ2F8J0zVpJfYj3K3095Q9 --phenofile  --sex 0 --inclusion_list file-GGJb118J80QpZGk5Fy4q6bJZ --transcript_index file-GFzk5gjJ0zVQXPQX8p4jj2BJ  --base_covariates file-GGPKqYQJJy8Gv2XKJ3v3KK85 --bgen_index file-GFzfyyQJ0zVqXZQBK4592QQg'

```


#### PheWAS

See below for how the input list of cancer phenotypes was created

```{bash}

dx run mrcepid-runassociationtesting --instance-type mem3_ssd1_v2_x16 --priority normal --destination results/ -imode='phewas' -ioutput_prefix="menopause_phewas" -iinput_args='--gene_ids BRCA2 CHEK2 ZNF518A --association_tarballs file-GGJZ2F8J0zVpJfYj3K3095Q9 --phenofile file-GGjgvfQJJy8F3Pvb1K8Qqx3x --sex 0 --inclusion_list file-GGbZFKQJP7J5223J4v6k49P3 --transcript_index file-GFzk5gjJ0zVQXPQX8p4jj2BJ  --base_covariates file-GGZkYk8JJy8GFjjF6kYG01g8 --sparse_grm file-GGbZPz8JP7J4vVXp4vj5JGyx --sparse_grm_sample file-GGbZPzjJP7JF98v34vxg8QkV --is_binary'

```

# Generating Medical Phenotypes

## Downloading And Processing Phenotypes

Medical phenotypes were generated in several steps:

1. Download first occurence variables using a JupyterLab notebook on DNANexus: `python_notebooks/download_first_occurence.ipynb`
2. Download summary diagnoses using a JupyterLab notebook on DNANexus: `python_notebooks/download_icd_data.ipynb`
    * Cancer ICD-10
    * Cancer ICD-9
    * HES ICD-10
    * HES ICD-9
    * Primary Death ICD-10
    * Secondary Death ICD_10 
3. Merge above information, add age of onset information, and convert from ICD-9 to ICD-10: `python_notebooks/icd9_to10.Rmd`
4. Remove duplicate codes and create final table of all conditions for all participants: `scripts/parse_ICD10_search_table.py`

The last step was run via the applet that I wrote called 'script_runner':

```{bash process phenos}

dx run script_runner -iscript=parse_ICD10_search_table.py -irequired_files=file-GG3jKp8JJy89gGbQ3Zb5561Q

```

## Creating Phenofiles

The script `scripts/parse_ICD10.py` can then be used to generate phenofiles for associationtesting:

```{bash}

# Binary phenotype for an ICD-10 grouping (note lack of 'Block' in -d flag):
# C81-C96 = Blood Cancers
# D46 = Myelodysplastic Syndrome
# B20-B24 = HIV/AIDS
./scripts/parse_ICD10.py -d C81-C96 D46 B20-B24 K743 K744 K745 K746 -i data_files/medical_data/coding19.tsv -p data_files/medical_data/vital_stats.tsv -c data_files/medical_data/processed_icd10.tsv

# Time-to-event phenotype for a single ICD-10 code:
# Recommended to use the --na flag for ANY tte analyses due to ambiguity in individuals with a code but without an incidence age
./scripts/parse_ICD10.py -d 'F29' -i data_files/medical_data/coding19.tsv -p data_files/medical_data/vital_stats.tsv -c data_files/medical_data/processed_icd10.tsv --tte --na

# Binary phenotype for multiple codes at once:
# We use NA here to set individuals with F42 but NOT F429 to NA as we do not know the subtype 
./scripts/parse_ICD10.py -d F20 F30 F429 -i data_files/medical_data/coding19.tsv -p data_files/medical_data/vital_stats.tsv -c data_files/medical_data/processed_icd10.tsv --na --output mhd.pheno

# All cancer codes:
./scripts/parse_ICD10.py -d "C00 C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C30 C31 C32 C33 C34 C37 C38 C39 C40 C41 C42 C43 C44 C45 C46 C47 C48 C49 C50 C51 C52 C53 C54 C55 C56 C57 C58 C60 C61 C62 C63 C64 C65 C66 C67 C68 C69 C70 C71 C72 C73 C74 C75 C76 C77 C78 C79 C80 C81 C82 C83 C84 C85 C86 C88 C90 C91 C92 C93 C94 C95 C96 C97" -i data_files/medical_data/coding19.tsv -p data_files/medical_data/vital_stats.tsv -c data_files/medical_data/processed_icd10.tsv -o cancer_codes.pheno

```

# Loading Association Tests

## Functions for Loading Data

```{r loading function}

load.and.plot.data <- function(file.names = c(), p.val.col, AC.col, tool.name, marker.file = NULL, ymin = 0, ymax = 15) {
  
  # Read all data files and cat
  result.table <- data.table()
  for (file in file.names) {
    curr.tab <- fread(file)
    result.table <- rbind(result.table,curr.tab)
  }
  
  # Set gene and p. value column to something standard:
  setnames(result.table, p.val.col, "p.value.selected")
  result.table[,log.p:=-log10(p.value.selected)] # Add -log10 p value
  result.table[,chrom:=factor(chrom, levels = mean.chr.pos[,chrom])] 
  
  masks <- result.table[MASK != "",unique(MASK)]
  mafs <- result.table[MAF != "",unique(MAF)]
  
  # Make standard plots
  plots = list()
  
  for (mask in masks) {
    for (maf in mafs) {
      name = paste(mask,maf,sep="-")
      manh.plot <- plot.manh(result.table[get(AC.col) > 10 & MASK == mask & MAF == maf], "log.p",ymin=ymin,ymax=ymax)
      qq.plot <- plot.qq(result.table[get(AC.col) > 10 & MASK == mask & MAF == maf], "p.value.selected",ymin=ymin,ymax=ymax)
      comb.plot <- manh.plot + qq.plot + plot_layout(ncol = 2, nrow = 1, widths = c(3,1.3)) + plot_annotation(title = paste(tool.name, mask, maf, sep = " - "), theme = theme(plot.title=element_text(size=18,face="bold",colour="black")))
      
      plots[[name]] = list('manh.plot' = manh.plot,
                           'qq.plot' = qq.plot,
                           'comb.plot' = comb.plot)
    }
  }
  
  if (!is.null(marker.file)) {
    return(list('gene.table' = result.table, 
                'marker.file' = marker.file,
                'plots' = plots,
                'masks' = masks))
  } else {
    return(list('gene.table' = result.table, 
                'marker.file' = NULL,
                'plots' = plots,
                'masks' = masks))
  }
  
}

```

## Functions for Plotting

```{r plotting functions}

plot.manh <- function(stats, p.var, ymin = 0, ymax = 15) {
 
  manh.plot <- ggplot(stats, aes(manh.pos, get(p.var), colour = chrom)) + 
    geom_point() +
    geom_hline(yintercept = -log10(1.6e-6), colour = "red", linetype = 2) +
    geom_text(inherit.aes = F, data = stats[get(p.var)>-log10(1.4e-6)], aes(manh.pos, get(p.var), label = SYMBOL), position = position_nudge(0.01,0.2),hjust=0,angle=45) +
    scale_x_continuous(name = "Chromosome", label = mean.chr.pos[,chrom], breaks = mean.chr.pos[,V1]) +
    scale_y_continuous(name = expression(bold(-log[10](italic(p)))), limits = c(ymin,ymax)) +
    scale_colour_manual(values = rep(c("#53878D","#7AC6CC"),12), breaks = c(as.character(1:22),"X","Y")) +
    coord_capped_cart(bottom="both",left="both") + # Comes from the "lemon" package
    theme + theme(panel.grid.major = element_blank())

  return(manh.plot)
}

plot.qq <- function(stats, p.var, ymin = 0, ymax = 15) {
  
  ## QQplot
  qqplot.data <- data.table(observed = stats[,get(p.var)],
                            SYMBOL = stats[,SYMBOL])
  setkey(qqplot.data,observed)
  qqplot.data <- qqplot.data[!is.na(observed)]
  qqplot.data[,observed:=-log10(observed)]
  qqplot.data[,expected:=-log10(ppoints(nrow(qqplot.data)))]
  
  qq.plot <- ggplot(qqplot.data, aes(expected, observed)) +
    geom_point() +
    geom_text(inherit.aes = F, data = qqplot.data[observed > -log10(1.4e-6)], aes(expected, observed, label = SYMBOL), position = position_nudge(-0.04,0),hjust=1) +
    geom_abline(slope = 1, intercept = 0,colour="red") +
    scale_x_continuous(name = expression(bold(Expected~-log[10](italic(p)))), limits = c(0,5)) +
    scale_y_continuous(name = "", limits = c(ymin,ymax)) +
    theme + theme(panel.grid.major = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  
  return(qq.plot)
  
}

```

## Example Code For Loading and Plotting

The following code block shows a simple example of how we load burden tests. Will need to modify p.val.col / tool.name / AC.col to adjust for the specific tool, but otherwise should be compatible with any burden tool that mrcepid-runassociationtesting is compatible with.

```{r load data, fig.height=6, fig.width=15}

# Actually load data and generate plots
bolt.ret <- load.and.plot.data(file.names = c("../scratch/menopause_34.bolt.genes.BOLT.stats.tsv.gz"),
                               p.val.col="P_BOLT_LMM",
                               tool.name = "BOLT",
                               AC.col = "AC",
                               marker.file = "../scratch/menopause_34.bolt.markers.BOLT.stats.tsv.gz",
                               ymax = 25)

# Show all the plots for MAF_01 (can change to AC_1/MAF_1/etc. for additional MAF cutoffs)
for (mask in names(bolt.ret$plots)[grepl("MAF_01", names(bolt.ret$plots))]) {
  print(bolt.ret$plots[[mask]]$comb.plot)
}
```
