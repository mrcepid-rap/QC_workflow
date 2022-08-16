---
title: "R Notebook"
output: html_notebook
---

# Introduction

The purpose of this Markdown is to document variant processing and QC of the 470k exome release of the UK Biobank data.

# Table of Contents

- [Setup](#setup)
- [Utility Functions](#utility-functions)
  * [Cloud Workstation {#cloud-workstation}](#cloud-workstation---cloud-workstation-)
  * [Querying jobs by status](#querying-jobs-by-status)
  * [Checking Finished Jobs](#checking-finished-jobs)
- [Variant Processing and Filtering](#variant-processing-and-filtering)
  * [1. Splitting BCFs](#1-splitting-bcfs)
  * [2. Filtering and VEP Annotation](#2-filtering-and-vep-annotation)
    + [2a. Post-processing of statistics](#2a-post-processing-of-statistics)
  * [3. Create BGEN from BCF](#3-create-bgen-from-bcf)
    + [3a. Make BGEN list file](#3a-make-bgen-list-file)
    + [3b. Run make-bgen](#3b-run-make-bgen)
  * [4. Collapse Variants](#4-collapse-variants)
  * [5. Run Association Testing](#5-run-association-testing)
    + [5a. Setting Genetic Data:](#5a-setting-genetic-data-)
      - [WES Samples List](#wes-samples-list)
      - [Define WBA grouping](#define-wba-grouping)
    + [5b. Setting Covariates](#5b-setting-covariates)
      - [Defining WES Batch](#defining-wes-batch)
    + [5c. Generate Transcript List](#5c-generate-transcript-list)
      - [Setting Relative Position of Genes in the Genome](#setting-relative-position-of-genes-in-the-genome)
      - [Gene Counts](#gene-counts)
      - [Identifying Genes with 0 Variants](#identifying-genes-with-0-variants)
      - [Sample Counts](#sample-counts)
      - [Writing Transcript File for DNANexus](#writing-transcript-file-for-dnanexus)
    + [5d. Running Associations](#5d-running-associations)
      - [Run the Tool](#run-the-tool)
      - [Running Individual Genes/Phenotypes](#running-individual-genes-phenotypes)
      - [Manual](#manual)
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

## Cloud Workstation {#cloud-workstation}

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

### 2a. Post-processing of statistics

Following this step I launched a [cloud workstation](#cloud-workstation) to generate final per-gene and per-individual statistics and generate a list file of variant coordinates within each split and annotated bcf-file for future users. The [cloud workstation](#cloud-workstation) was launched as described above. Following launch, I used the following commands:

```{bash, eval = F}

# Download all statistics files from the 'annotate-cadd' step into two separate directories:
mkdir per_indv/
mkdir per_gene/
cd per_indv/ \
  && dx download --no-progress filtered_vcfs/ukb23148_c*_b*_v1_chunk*.norm.filtered.tagged.missingness_filtered.annotated.cadd.per_indv.tsv \
  && cd ../per_gene/ \
  && dx download --no-progress filtered_vcfs/ukb23148_c*_b*_v1_chunk*.norm.filtered.tagged.missingness_filtered.annotated.vep.tsv.gz

# Download required processing files
dx download project_resources/genetics/ukb_450k_samples.txt
dx download mash.pl

# Merge per-individual statistics files together:
ls per_indv/*.per_indv.tsv > per_indv_list.txt
./mash.pl > 450k_indv_counts.tsv
bgzip 450k_indv_counts.tsv
dx upload --destination results/ 450k_indv_counts.tsv.gz

# Merge per-variant statistics files together:
zcat per_gene/*.tsv.gz > 450k_vep.tsv
sort -k1,1 -k2,2n 450k_vep.tsv > 450k_vep.sorted.tsv
bgzip 450k_vep.sorted.tsv
tabix -s 1 -b 2 -e 2 450k_vep.sorted.tsv.gz
 dx upload --destination results/ 450k_vep.sorted.tsv.gz
 
# Generate coordinate file:
ls per_gene/*.tsv.gz > variants_list.txt
./coordinate_maker.pl > coordinates.tsv
dx upload coordinates.tsv

# Then locally...
dx ls -l filtered_vcfs/*.norm.filtered.tagged.missingness_filtered.annotated.vep.tsv.gz.tbi | perl -ane 'chomp $_; if ($F[5] =~ /(ukb23148_c[0-9XY]{1,2}_b\d+_v1_chunk\d)\./) {$root = $1; if ($F[6] =~ /^\((\S+)\)$/) {$id = $1; print "$root\t$id\t$F[5]\n"}}' > file_locs.txt
~/Documents/Current\ Projects/Fertility/UKBBFertility/scripts/matcher.pl -file2 coordinates.files.tsv -file1 file_locs.txt -r -col2 3 | perl -ane 'chomp $_; print "$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[5]\t$F[6]\n";' | sort -k1,1 -k2,2n -k3,3n > coordinates.files.tsv
bgzip coordinates.files.tsv
tabix -s 1 -b 2 -e 3 coordinates.files.tsv.gz
dx upload coordinates.files.tsv.gz
```

## 3. Create BGEN from BCF

### 3a. Make BGEN list file

Need to login to a cloud workstation (see above) and download all of the coordinates files from each original job in step (2). I then make a single .tsv.gz file listing the file location of all filtered & annotated bcf/vep in the following format:

    #chrom    start   stop    fileprefix    bcf_dxpy    vep_dxpy
    chr1    1   1000000   ukb23157_c1_b1_v1_chunk1    file-1234567890ABCDEFGHIJK    file-0987654321KJIHGFEDCBA
    chr1    1000001   2000000   ukb23157_c1_b1_v1_chunk2    file-ABCDEFGHIJK1234567890    file-KJIHGFEDCBA0987654321

Commands below make this file on the cloud workstation:

```{bash, eval = F}

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

```{bash, eval = F}

perl -e 'for my $chr (1..22,"X","Y") {print "dx run mrcepid-makebgen --destination filtered_bgen/ --priority normal -ichromosome=\"chr$chr\" -icoordinate_file=file-GFx0z70J0zVpJzVF0QG6Fjfk;\n";}' | bash

```

```{r}

bcf_coords <- fread("../scratch/bcf_coordinates.tsv.gz")
bcf_coords[,dummy:=1]
bcf_coords[,row:=.I]

bcf_coords[,root_name:=str_split(fileprefix, "_chunk",simplify = T)[,1]]
bcf_coords[,sum(dummy),by=root_name]

starts = bcf_coords[,start]
stops = bcf_coords[,stop]

distances <- starts[2:length(starts)] - stops[1:(length(stops)-1)]
distances <- c(NA, distances)
bcf_coords[,distance:=distances]

vep <- fread("../scratch/chr1.filtered.vep.tsv.gz")
```

## 4. Collapse Variants

```{bash, eval = F}
# MAF < 1e-3
# PTVs
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & PARSED_CSQ=="PTV" & LOFTEE=="HC"' -ifile_prefix=HC_PTV-MAF_01 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & PARSED_CSQ=="PTV"' -ifile_prefix=PTV-MAF_01 --destination collapsed_variants_new/ --brief --yes

# Missense
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE"' -ifile_prefix=MISS-MAF_01 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & CADD>=25' -ifile_prefix=MISS_CADD25-MAF_01 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & REVEL>=0.5' -ifile_prefix=MISS_REVEL0_5-MAF_01 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & REVEL>=0.7' -ifile_prefix=MISS_REVEL0_7-MAF_01 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & CADD>=25 & REVEL>=0.7' -ifile_prefix=MISS_CADD25_REVEL0_7-MAF_01 --destination collapsed_variants_new/ --brief --yes

# Synonymous
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & PARSED_CSQ=="SYN"' -ifile_prefix=SYN-MAF_01 --destination collapsed_variants_new/ --brief --yes

# Combined Deleterious
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AF<0.001 & FILTER=="PASS" & ((PARSED_CSQ=="PTV" & LOFTEE=="HC") | (PARSED_CSQ=="MISSENSE" & CADD>=25))' -ifile_prefix=DMG-MAF_01 --destination collapsed_variants_new/ --brief --yes

# AC = 1
# PTVs
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & PARSED_CSQ=="PTV" & LOFTEE=="HC"' -ifile_prefix=HC_PTV-AC_1 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & PARSED_CSQ=="PTV"' -ifile_prefix=PTV-AC_1 --destination collapsed_variants_new/ --brief --yes

# Missense
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE"' -ifile_prefix=MISS-AC_1 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & CADD>=25' -ifile_prefix=MISS_CADD25-AC_1 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & REVEL>=0.5' -ifile_prefix=MISS_REVEL0_5-AC_1 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & REVEL>=0.7' -ifile_prefix=MISS_REVEL0_7-AC_1 --destination collapsed_variants_new/ --brief --yes
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & CADD>=25 & REVEL>=0.7' -ifile_prefix=MISS_CADD25_REVEL0_7-AC_1 --destination collapsed_variants_new/ --brief --yes

# Synonymous
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & PARSED_CSQ=="SYN"' -ifile_prefix=SYN-AC_1 --destination collapsed_variants_new/ --brief --yes

# Combined Deleterious
dx run mrcepid-collapsevariants --priority low -ifiltering_expression='AC==1 & FILTER=="PASS" & ((PARSED_CSQ=="PTV" & LOFTEE=="HC") | (PARSED_CSQ=="MISSENSE" & CADD>=25))' -ifile_prefix=DMG-AC_1 --destination collapsed_variants_new/ --brief --yes
```

```{bash}

dx ls -l filtered_vcfs/ukb23148_c*_b*_v1_chunk?.norm.filtered.tagged.missingness_filtered.annotated.cadd.bcf | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "$1\n";}' > collapse_list.txt

# Split this file list into individual "jobs":
# This will generate 134 files with ~35 BCFs per file. This means that each core will be ~used 2x in each job
split -a 2 -l 35 collapse_list.txt collapse_list_

# Upload these lists onto DNA Nexus:
dx upload collapse_list_* --destination batch_lists/

# Generate a list of jobs to run on DNA Nexus, and submit
#MAF < 1e-3
# PTVs
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.001 & LOFTEE=\"HC\" & PARSED_CSQ=\"PTV\"'\'' -ifile_prefix='\''HC_PTV-MAF_01'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.001 & PARSED_CSQ=\"PTV\"'\'' -ifile_prefix='\''PTV-MAF_01'\'';\n";}' | bash

# Missense
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority normal --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.001 & PARSED_CSQ=\"MISSENSE\"'\'' -ifile_prefix='\''MISS-MAF_01'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority normal --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.001 & PARSED_CSQ=\"MISSENSE\" & CADD>=25'\'' -ifile_prefix='\''MISS_CADD25-MAF_01'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.001 & PARSED_CSQ=\"MISSENSE\" & REVEL>=0.5'\'' -ifile_prefix='\''MISS_REVEL0_5-MAF_01'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.001 & PARSED_CSQ=\"MISSENSE\" & REVEL>=0.7'\'' -ifile_prefix='\''MISS_REVEL0_7-MAF_01'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.001 & PARSED_CSQ=\"MISSENSE\" & REVEL>=0.7 & CADD>=25'\'' -ifile_prefix='\''MISS_CADD25_REVEL0_7-MAF_01'\'';\n";}' | bash

# Synonymous Variants
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority normal --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.001 & PARSED_CSQ=\"SYN\"'\'' -ifile_prefix='\''SYN-MAF_01'\'';\n";}' | bash

# AC = 1
# PTVs
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AC=1 & LOFTEE=\"HC\" & PARSED_CSQ=\"PTV\"'\'' -ifile_prefix='\''HC_PTV-AC_1'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AC=1 & PARSED_CSQ=\"PTV\"'\'' -ifile_prefix='\''PTV-AC_1'\'';\n";}' | bash

# Missense
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AC=1 & PARSED_CSQ=\"MISSENSE\"'\'' -ifile_prefix='\''MISS-AC_1'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AC=1 & PARSED_CSQ=\"MISSENSE\" & CADD>=25'\'' -ifile_prefix='\''MISS_CADD25-AC_1'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AC=1 & PARSED_CSQ=\"MISSENSE\" & REVEL>=0.5'\'' -ifile_prefix='\''MISS_REVEL0_5-AC_1'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AC=1 & PARSED_CSQ=\"MISSENSE\" & REVEL>=0.7'\'' -ifile_prefix='\''MISS_REVEL0_7-AC_1'\'';\n";}' | bash
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AC=1 & PARSED_CSQ=\"MISSENSE\" & REVEL>=0.7 & CADD>=25'\'' -ifile_prefix='\''MISS_CADD25_REVEL0_7-AC_1'\'';\n";}' | bash

# Synonymous Variants
dx ls -l batch_lists/collapse_list_* | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination collapsed_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AC=1 & PARSED_CSQ=\"SYN\"'\'' -ifile_prefix='\''SYN-AC_1'\'';\n";}' | bash
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
    
C. 

### 5a. Setting Genetic Data:

Building the set of files in item (1) above, also requires some pre-build data. We document how to generate these files below.

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

-   PC1-40 (field 22009-0.[1-40])
-   Genetic ancestry grouping (field 22006-0.0)

This file was then queried locally. All we have done is calculate an ellipse based on MAD\*5 for two given PCs in three combinations:

-   PC1 + PC2
-   PC3 + PC4
-   PC5 + PC6

Individuals that fall within this ellipse for all three cases above (as similarly done in Bycroft et al.) are included in the 'european' ancestry definition for the purposes of association testing. These resulting file (`wba.txt`) is then uploaded onto DNANexus into the `project_resources/genetics/` folder for use with the `mrcepid-makegrms` applet.

```{r pca}

pca <- fread("../scratch/pca_participant.tsv")

cols <- c('eid', paste('22009-0',seq(1,10), sep = '.'), '22006-0.0')
setnames(pca, cols, c('n_eid', paste0('PC', seq(1,10)), 'wba'))

pca <- pca[!is.na(PC1)]
pca[,wba:=if_else(is.na(wba), 0, 1)]

calc_ellipse <- function(one, two, mad.pc1, mad.pc2, median.pc1, median.pc2) {
  val <- (((one - median.pc1) ^ 2 / (mad.pc1*5)^2) + (((two - median.pc2)^2) / (mad.pc2*5)^2)) <= 1
  return(val)
}


calc_pc <- function(pc_one, pc_two) {
  
  mad.pc1 <- mad(pc_one)
  mad.pc2 <- mad(pc_two)
  
  median.pc1 <- median(pc_one)
  median.pc2 <- median(pc_two)
  
  ellipse_table <- data.table(one=pc_one, two=pc_two)
  ellipse_table[,result:=calc_ellipse(one, two, mad.pc1, mad.pc2, median.pc1, median.pc2)]
  
  return(ellipse_table[,result])
  
}

pca[,pc1_pc2:=calc_pc(pca[,PC1], pca[,PC2])]
pca[,pc3_pc4:=calc_pc(pca[,PC3], pca[,PC4])]
pca[,pc5_pc6:=calc_pc(pca[,PC5], pca[,PC6])]

ggplot(pca,aes(PC1, PC2, colour = pc1_pc2)) + geom_point(size = 0.2) + theme
ggplot(pca,aes(PC3, PC4, colour = pc3_pc4)) + geom_point(size = 0.2) + theme
ggplot(pca,aes(PC5, PC6, colour = pc5_pc6)) + geom_point(size = 0.2) + theme

pca[,wba_v2:=pc1_pc2 & pc3_pc4 & pc5_pc6]
pca[,European_ancestry:=if_else(wba_v2==TRUE, 1, 0)]

fwrite(pca[,c("n_eid","European_ancestry")], '../scratch/wba.txt', quote = F, sep = ' ', row.names = F, col.names = T)
```

### 5b. Setting Covariates

#### Defining WES Batch

```{r WES Batch}

covariates <- fread("../scratch/covariates.tsv")
covariates[,eid:=as.character(eid)]

load.samples <- function(file) {
  
  samples <- fread(file,header = F)
  setnames(samples,"V1","eid")
  samples[,eid:=as.character(eid)]
  return(samples)
  
}

samples.50k <- load.samples("../scratch/samples_50k.txt")
samples.200k <- load.samples("../scratch/samples_200k.txt")
samples.450k <- load.samples("../scratch/samples_450k.txt")

covariates[,wes.batch:=if_else(eid %in% samples.50k[,eid], "50k",
                               if_else(eid %in% samples.200k[,eid], "200k",
                                       if_else(eid %in% samples.450k[,eid], "450k", as.character(NA))))]

fwrite(covariates, "../scratch/covariates.wes_added.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

```

### 5c. Generate Transcript List

To get a reasonable set of protein-coding transcripts for Hg38 I generated two separate files:

1.  `data_files/hgTables.txt` - File downloaded from the UCSC Genome Browser - Table Browser. The purpose of this file is to calculate total coding sequence length, which cannot be derived from BioMart (file #2). I selected GENCODE v38 and selected the following rows:

    -   Transcript name
    -   Coding sequence start
    -   Coding sequence end
    -   Number of exons
    -   Exon start coordinates
    -   Exon end coordinates
    -   transcript type (e.g. protein_coding)

This file is then processed like: `scripts/get_cds_len.pl > ./data_files/cds_lengths.txt`

2.  `data_files/mart_export` - File downloaded from ENSEMBL BioMart. The purpose of this file is to get valid ENSEMBL transcripts for all possible genes in Hg38, with canonical and MANE transcript for each gene. I selected GRCh38 and selected the following rows:

    -   Gene stable ID (form of ENSG)
    -   Transcript stable ID (form of ENST)
    -   MANE Select ID (form of NM)
    -   Transcript Length - Note that this is **NOT** coding sequence length, as above
    -   Gene Name (I believe HUGO name)
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

# Also load the raw exon models in order to plot later:
exon.models <- fread("data_files/hgTables.txt")
setnames(exon.models, names(exon.models), c("transcript","cdsStart","cdsEnd","numExons","exonStarts","exonEnds","transcriptType"))
exon.models[,transcript:=str_split(transcript,"\\.",simplify = T)[1],by=1:nrow(exon.models)]

```

#### Setting Relative Position of Genes in the Genome

This is just to enable easy plotting on Manhattan plots.

```{r mean gene position}

transcripts <- fread("data_files/transcripts.tsv.gz")
setnames(transcripts,"#chrom","chrom")

# Get chromosome locations for the plot
mean.chr.pos <- transcripts[,mean(manh.pos), by = chrom]

```

#### Gene Counts

Then, I use data calculated [above](#post-processing-of-statistics), to check expected number of rare synonymous variants per gene as a function of coding sequence length:

```{r fig.height=6, fig.width=8}

vep.stats <- fread("data_files/450k_vep.sorted.tsv.gz")
setnames(vep.stats,"#CHROM","CHROM")
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
fwrite(transcripts.annotated[syn.count == 0], "data_files/missing_genes.tsv", col.names = F, row.names = F, sep = "\t")

# Checking number of baits for missing genes
# perl -ane 'chomp $_; print "echo $F[0]\ntabix xgen_plus_spikein.GRCh38.bed.gz \"$F[9]\" | wc -l;\n";' missing_genes.tsv | bash > bait_counts.txt
# perl -ane 'chomp $_; if ($_ =~ /ENST/) {print "$_\t"} elsif ($_ =~ /(\d+)/) {print "$1\n"}' bait_counts.txt > bait_counts.form.txt
baits <- fread("data_files/bait_counts.form.txt")
setnames(baits, c("V1","V2"),c("ENST","bait.count"))
missing <- merge(missing,baits,by="ENST")

# Checking number of IDd variants for missing genes
# perl -ane 'chomp $_; print "echo $F[0]\ntabix 450k_vep.sorted.tsv.gz \"$F[9]\" | wc -l;\n";' missing_genes.tsv | bash > variant_counts.txt
# perl -ane 'chomp $_; if ($_ =~ /ENST/) {print "$_\t"} elsif ($_ =~ /(\d+)/) {print "$1\n"}' variant_counts.txt > variant_counts.form.txt
variants <- fread("data_files/variant_counts.form.txt")
setnames(variants, c("V1","V2"),c("ENST","variant.count"))
missing <- merge(missing,variants, by = "ENST")

# Checking proportion of variants in a given gene with FAIL for FILTER
fail.prop <- data.table(pivot_wider(vep.stats[,sum(dummy),by = c("ENST","FILTER")],id_cols = "ENST", names_from = "FILTER", values_from = "V1", values_fill = 0))
fail.prop[,fail.prop:=FAIL / (FAIL + PASS)]
missing <- merge(missing, fail.prop[,c("ENST", "fail.prop")], all.x = T, by = "ENST")

# Setting fail categories
missing[,fail.cat:=if_else(bait.count == 0, "NOT_SEQ", "")]
missing[,fail.cat:=if_else(grepl("-",SYMBOL) & fail.cat == "", "READTHROUGH", fail.cat)]
missing[,fail.cat:=if_else(chr == "Y" & fail.cat == "", "CHR_Y", fail.cat)]
missing[,fail.cat:=if_else(SYMBOL == "" & fail.cat == "", "NO_NAME", fail.cat)]
missing[,fail.cat:=if_else(fail.prop > 0.75 & fail.cat == "" & !is.na(fail.prop), "NO_PASS_VARS", fail.cat)]
missing[,fail.cat:=if_else(variant.count == 0 & fail.cat == "", "NO_VARS_IN_VCF",fail.cat)]

# Manually annotating the rest...
fwrite(missing[fail.cat == ""],"data_files/missing.tsv",col.names = T, row.names = F, sep = "\t", quote = F)

# And read back in and add to original file:
missing.annotated <- fread("data_files/missing.annotated.txt")
# Everything after here is stupid AF. Why did I code it this way?
missing.annotated <- merge(missing,missing.annotated[,c("ENST","fail.cat")],by="ENST")
missing.annotated <- missing.annotated[fail.cat.x==""]
setnames(missing.annotated,"fail.cat.y","fail.cat")
missing.annotated[,fail.cat.x:=NULL]
missing <- rbind(missing[fail.cat != ""], missing.annotated)

# Merge back into the main plot
transcripts.annotated <- merge(transcripts.annotated,missing[,c("ENST","fail.cat")],by="ENST",all.x=T)
transcripts.annotated[,dummy:=1]
transcripts.annotated[,fail:=if_else(is.na(fail.cat),F,T)]
per.chr.fail <- data.table(pivot_wider(transcripts.annotated[,sum(dummy),by=c("fail","chr")], names_from = fail, values_from = V1))
per.chr.fail[,chr:=factor(chr,levels=c(1:22,"X","Y"))]
setkey(per.chr.fail,chr)
per.chr.fail[,prop.fail:=(`TRUE`/(`TRUE` + `FALSE`))*100]

ggplot(per.chr.fail,aes(chr,prop.fail)) + geom_col() + xlab("Chromosome") + ylab("Proportion Failed Genes") + theme + theme(panel.grid.major.x=element_blank())

per.cat.fail <- transcripts.annotated[!is.na(fail.cat),sum(dummy),by=c("fail.cat","chr")]
per.cat.fail[,chr:=factor(chr,levels=c(1:22,"X","Y"))]
ggplot(per.cat.fail,aes(chr, V1, fill=fail.cat, group=fail.cat)) + geom_col() + theme.legend + theme(panel.grid.major.x=element_blank())
```

#### Sample Counts

```{r}

wba <- fread("data_files/UKB_European_IncList_15092021.txt")
setnames(wba,"n_eid","eid")
wba[,eid:=as.character(eid)]

indv.counts <- fread("data_files/450k_indv_counts.tsv.gz")
indv.counts[,eid:=as.character(eid)]

indv.counts <- merge(indv.counts,wba[,c("eid","European_ancestry")],by="eid")

pois.dist <- data.table(table(rpois(n = 453342, lambda = 77)))
pois.dist[,V1:=as.integer(V1)]

ggplot(indv.counts,aes(SYN)) + geom_histogram(binwidth=1) + geom_vline(xintercept=102) + theme

ggplot(indv.counts[European_ancestry == 1],aes(PTV)) + geom_histogram(binwidth=1) + theme
ggplot(indv.counts[European_ancestry == 1],aes(MISSENSE)) + geom_histogram(binwidth=1) + theme
ggplot(indv.counts[European_ancestry == 1],aes(SYN)) + geom_histogram(binwidth=1) + xlim(0,160) + theme

quantile(rpois(n = 453342, lambda = 77), probs = c(0.997))
indv.counts[,high.SYN:=if_else(SYN>quantile(rpois(n = 453342, lambda = 77), probs = c(0.995))[[1]],1,0)]
table(indv.counts[,c("high.SYN","European_ancestry")])
```

#### Writing Transcript File for DNANexus

Finally, I write a file that is uploaded to DNANexus for annotation during Association Testing. This file is currently stored in project `project-G6BJF50JJv8p4PjGB9yy7YQ2` as file `file-G7xyzF8JJv8kyV7q5z8VV3Vb`.

```{r}



```

### 5d. Running Associations

This uses the 'launch' script located in `./scripts/`. See that file for more information on how association testing is run for individual tools. This script will launch all tools on all masks for a given phenotype.

```{bash, eval = F}

# Make a list of all association test tar files:
dx ls -l collapsed_variants_new/*.tar.gz |  grep -v 'MISS-' | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "$1\n";}' > variant_mask_list.txt
dx upload variant_mask_list.txt --destination collapsed_variants_new/

## TEST PHENO:
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7jpjZ0JJv8jG8FKJGz7v1JZ false LOY_test 1 file-G7vPjPjJYVkB6qG4PbKz21gP loy_3way_GRS

# AAM
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7p0kZ0JYVk89zK628525Vkp false AAM_all 0 file-G7v19X8JYVk96gFpBKkpJ3pP menarche_GRS_696
# AAM - Capped 8/19
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7p0kjQJYVkG4Zg72Qg0bPBy false AAM_8_19 0 file-G7v19X8JYVk96gFpBKkpJ3pP menarche_GRS_696 

# LOY - 3Way
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7jpjZ0JJv8jG8FKJGz7v1JZ false LOY_3way 1 file-G7vPjPjJYVkB6qG4PbKz21gP loy_3way_GRS
# LOY - 2Way
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G529BXjJJjbkBxx1K3V86yVf false LOY_2way 1 file-G7vPjG0JYVkPJb35569p07Xk loy_2way_GRS
# LOX - 3Way
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G81x4ZjJYVk719K19YXV2FF5 false LOX_3way 0 file-XXXXXXXXXXXXXXXXXXXXXXXX lox_2way_GRS
# LOX - 2Way
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G857x90JYVk00Vj14V222px7 false LOX_2way 0 file-XXXXXXXXXXXXXXXXXXXXXXXX lox_2way_GRS

# T2D
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G5PFv00JXk80Qfx04p9X56X5 true T2D 2 file-G7vPj7QJYVk877BFP6F6VQKK t2d_unadjusted_GRS
# T2D Shrinkage GRS
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G5PFv00JXk80Qfx04p9X56X5 true T2D_shrinkage 2 file-G81j2x8JYVk9bX5qGy7V3z8V SCORE
# BMI
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7k6GK8JJv8vXgf9PZz3K0G9 false BMI 2 file-G7vPgyQJYVkPyFxv8Fg0gPZx bmi_GRS
# BMI_male
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7k6GK8JJv8vXgf9PZz3K0G9 false BMI 1 file-G7vPgyQJYVkPyFxv8Fg0gPZx bmi_GRS
# BMI_female
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7k6GK8JJv8vXgf9PZz3K0G9 false BMI 0 file-G7vPgyQJYVkPyFxv8Fg0gPZx bmi_GRS

# Height
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G80Y61QJYVk71gZgPgBG8K2G false height 2 file-XXXXXXXXXXXXXXXXXXXXXXXX height_GRS
# IGF1 
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G815vg0JYVk8F6bj5y9303ZP false IGF1 2 file-XXXXXXXXXXXXXXXXXXXXXXXX IGF1_GRS
# Birthweight
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G81vk70JYVk87bgq5Qq8z0BV false birthweight 2 file-XXXXXXXXXXXXXXXXXXXXXXXX birthweight
# birthweight_capped
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G81vkK8JYVk4q85J97vbbKxy false birthweight_capped 2 file-XXXXXXXXXXXXXXXXXXXXXXXX birthweight

# vb_all
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G876Z98J5K53Jjgj1kG5fJkk false vb_all 1 file-XXXXXXXXXXXXXXXXXXXXXXXX vb

# sac10_all
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G876Z98J5K56qvzK46JVBXX5 false sac10_all 2 file-XXXXXXXXXXXXXXXXXXXXXXXX sac

# neb_males
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G85ggVQJYVk7jGk0GBfjF1j3 false neb_males 1 file-XXXXXXXXXXXXXXXXXXXXXXXX children_fathered
# neb_females
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G85ggVQJYVk7jGk0GBfjF1j3 false neb_females 0 file-XXXXXXXXXXXXXXXXXXXXXXXX live_births
#neb
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G85ggVQJYVk7jGk0GBfjF1j3 false neb 2 file-XXXXXXXXXXXXXXXXXXXXXXXX neb

# childlessness_males
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G85ggZQJYVk88352Fqg3FpGF true childlessness_males 1 file-XXXXXXXXXXXXXXXXXXXXXXXX childlessness
# childlessness_females
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G85ggZQJYVk88352Fqg3FpGF true childlessness_females 0 file-XXXXXXXXXXXXXXXXXXXXXXXX childlessness
# childlessness
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G85ggZQJYVk88352Fqg3FpGF true childlessness 2 file-XXXXXXXXXXXXXXXXXXXXXXXX childlessness

# menopause_all
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G82b3v0JYVkJ7GXfJ3y54qkB false menopause_all 0 file-XXXXXXXXXXXXXXXXXXXXXXXX menopause_grs
# menopause_34min
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G82b4xjJYVk24bfJ47QP4F87 false menopause_34min 0 file-XXXXXXXXXXXXXXXXXXXXXXXX menopause_grs
# menopause_40to60
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G82b558JYVk567fy5V8Z01fZ false menopause_40to60 0 file-XXXXXXXXXXXXXXXXXXXXXXXX menopause_grs
# menopause_POI
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G82b58QJYVkJb1J03Bq2Xpj3 true menopause_POI 0 file-XXXXXXXXXXXXXXXXXXXXXXXX menopause_grs
# menopause_invnorm
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G839B8jJJv8X2yFvKK5xKBQv false menopause_invnorm 0 file-XXXXXXXXXXXXXXXXXXXXXXXX menopause_grs

#WHR adjBMI
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G80Z29QJYVk3Z8kzPg0PvpV1 false whr_bmi 2 file-XXXXXXXXXXXXXXXXXXXXXXXX whr_grs
# Depression
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7vPy4QJYVk907JP2Q0fkxFZ true depression 2 file-XXXXXXXXXXXXXXXXXXXXXXXX deppression_grs
# probable_MDD
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7vPy9QJYVk7QFx33B7z7gpb true probably_MDD 2 file-XXXXXXXXXXXXXXXXXXXXXXXX deppression_grs

# BMD
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7y1y6QJYVk44v1j6FxkZ946 false heelBMD 2 file-XXXXXXXXXXXXXXXXXXXXXXXX bmd_grs
# RBC
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G7y1yJjJYVk8Ygpg6JBVgFX7 false rbc 2 file-XXXXXXXXXXXXXXXXXXXXXXXX rbc_grs

# Lives with Partner
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G87q0f0JYVkGVv64BkgQk09x true lives_with_partner 2 file-XXXXXXXXXXXXXXXXXXXXXXXX partner_grs

# Adjusted_TL
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G81x9kjJYVkKPQ2F3gq028QV false adjusted_tl 2 file-XXXXXXXXXXXXXXXXXXXXXXXX adjusted_tl

# Baby birthweight
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G87xjg8JYVk56qfYBqG2Yg1K false babys_bw_1to7 0 file-XXXXXXXXXXXXXXXXXXXXXXXX bw_grs
./launch.sh file-G7zPvZ0JJv8v06j8Gv2ppxpJ file-G87xkX8JYVk8y19F3ky9ZJ3P false babys_bw_2.5to4.5 0 file-XXXXXXXXXXXXXXXXXXXXXXXX bw_grs
```

#### Run the Tool

```{bash}

dx run mrcepid-buildgrms --priority low --destination project_resources/genetics/ -isample_ids_file=file-GFjJ0yjJ0zVb9BK82ZQFkx0y -iwba_file=file-GFjX7y8J0zVgB7JbPq7K23qP

```

#### Running Individual Genes/Phenotypes

#### Manual

```{bash}

dx run mrcepid-extractvariants --destination results/gene_specific/ -iassociation_tarball=collapsed_variants_new/MISS_CADD25-MAF_01.tar.gz -iphenofile=file-G82b4xjJYVk24bfJ47QP4F87 -iis_binary=false -ioutput_prefix=MISS_CADD25.VARS1.T2D -igene_id=ENST00000375663 -iinclusion_list=file-G6qXvjjJ2vfQGPp04ZGf6ygj --name VARS1.anm_all --yes --brief

dx run mrcepid-extractvariants --destination results/gene_specific/ -iassociation_tarball=collapsed_variants_new/HC_PTV-MAF_01.tar.gz -iphenofile=file-G7jpjZ0JJv8jG8FKJGz7v1JZ -iis_binary=false -ioutput_prefix=HC_PTV.GIGYF1.T2D -igene_id=GIGYF1 -iinclusion_list=file-G6qXvjjJ2vfQGPp04ZGf6ygj --name GIGYF1.anm_all --yes --brief

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
bolt.ret <- load.and.plot.data(file.names = c("ukbb_data/T2D_ExWAS/T2D.bolt.genes.BOLT.stats.tsv.gz"),
                               p.val.col="P_BOLT_LMM_INF",
                               tool.name = "BOLT",
                               AC.col = "AC",
                               marker.file = "ukbb_data/T2D_ExWAS/T2D.bolt.markers.BOLT.stats.tsv.gz",
                               ymax = 60)

# Show all the plots for MAF_01 (can change to AC_1/MAF_1/etc. for additional MAF cutoffs)
for (mask in names(bolt.ret$plots)[grepl("MAF_01", names(bolt.ret$plots))]) {
  print(bolt.ret$plots[[mask]]$comb.plot)
}

```
